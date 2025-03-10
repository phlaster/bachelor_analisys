#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)

using ArgParse
using Dates
using Logging
using DelimitedFiles
using ProgressMeter
using .Threads: nthreads, @threads, threadid

function parse_commandline()
    s = ArgParseSettings(prog="download_genomes.jl",
        description="""Download genomes for all accessions in the chosen file as separate zip archives""")
    @add_arg_table! s begin
        "--input-file", "-i"
            help = "Path to accessions .tsv file"
            required = true
        "--output-dir", "-o"
            help = "Output directory (will be created if absent)"
            required = true
        "--max-retries", "-r"
            help = "Maximum download attempts per accession"
            arg_type = Int
            default = 3
        "--extract", "-e"
            help = "Extract downloaded archives"
            action = :store_true
        "--force", "-f"
            help = "Force archives extraction"
            action = :store_true
        "--parallel", "-p"
            help = "Number of simultaneous downloads \
            !!! to run unzip step with N parallel threads run the whole script as `julia --threads N script.jl [args]`"
            arg_type = Int
            default = 1
    end
    return parse_args(s)
end

# Download a single accession
function download_accession(accession::AbstractString, output_dir::AbstractString, max_retries::Int, min_duration::Ref{Float64})
    archive_dir = joinpath(output_dir, "archives")
    mkpath(archive_dir)
    archive_path = joinpath(archive_dir, "$accession.zip")
    archive_part = "$archive_path.part"

    attempts = 0
    success = false
    duration = 0.0

    while attempts < max_retries && !success
        attempts += 1
        rm(archive_part, force=true)

        @info "Attempt $attempts of $max_retries for accession $accession"
        start_time = now()
        try
            cmd = pipeline(`datasets download genome accession $accession --include genome,gff3 --filename $archive_part`, devnull)
            run(cmd)
            duration = (now() - start_time).value / 1000

            if min_duration[] != -1 && duration > 3 * min_duration[]
                @warn "Accession $accession: Download took $(duration)s (> $(3 * min_duration[])s), retrying..."
                continue
            end

            if min_duration[] == -1 || duration < min_duration[]
                min_duration[] = duration
            end

            mv(archive_part, archive_path, force=true)
            success = true
            @info "Accession $accession: Downloaded successfully in $(duration)s"
        catch e
            @warn "Accession $accession: Attempt $attempts failed: $e"
            if attempts == max_retries
                @error "Accession $accession: Failed after $max_retries attempts"
            end
        end
    end

    return success
end

function extract_archives(output_dir::String; force=false)
    archive_dir = joinpath(output_dir, "archives")
    genome_dir = joinpath(output_dir, "genomes")
    mkpath(genome_dir)
    
    zips = filter(f -> endswith(f, ".zip"), readdir(archive_dir, join=true))
    prog = Progress(length(zips))
    @threads for zip in zips
        accession = splitext(basename(zip))[1]
        extraction_dir = joinpath(genome_dir, accession)
        mkpath(extraction_dir)
        target_fna = joinpath(extraction_dir, "$accession.fna")
        target_gff = joinpath(extraction_dir, "$accession.gff")
        target_json = joinpath(extraction_dir, "$accession.json")

        if any(isfile, [target_fna, target_gff, target_json]) && !force
            error("$extraction_dir has name collisions")
            exit()
        end
        
        # Step 1: Extract the archive into a temporary directory
        temp_dir = mktempdir()
        try
            run(pipeline(`unzip -o $zip -d $temp_dir`, devnull))
            
            # Step 2: Locate and move the required files
            data_dir = joinpath(temp_dir, "ncbi_dataset", "data", accession)
            if isdir(data_dir)
                # Find and move the FASTA file
                fna_files = filter(f -> endswith(f, "_genomic.fna"), readdir(data_dir))
                if length(fna_files) == 1
                    original_fna = joinpath(data_dir, fna_files[1])
                    mv(original_fna, target_fna, force=true)
                else
                    @warn "Expected one _genomic.fna file in $data_dir, found $(length(fna_files))"
                end
                
                # Move and rename the GFF file
                original_gff = joinpath(data_dir, "genomic.gff")
                if isfile(original_gff)
                    mv(original_gff, target_gff, force=true)
                else
                    @warn "GFF file missing in $data_dir"
                end
            else
                @warn "Data directory not found for $accession"
            end
            
            # Move and rename the JSON assembly report
            json_file = joinpath(temp_dir, "ncbi_dataset", "data", "assembly_data_report.jsonl")
            if isfile(json_file)
                mv(json_file, target_json, force=true)
            else
                @warn "Assembly report JSON missing for $accession"
            end
            
            # Step 3: Clean up temporary directory
            rm(temp_dir, recursive=true)
        catch e
            @error "Failed to process $zip for $accession: $e (Thread $(threadid()))"
            rm(temp_dir, recursive=true, force=true)
        end
        next!(prog)
    end
    @info "Extraction and reorganization completed"
end

# Main function
function main()
    global_logger(ConsoleLogger(stdout, Logging.Info))
    args = parse_commandline()
    check_dependencies(["datasets", "unzip"])

    input_file = abspath(args["input-file"])
    output_dir = abspath(mkpath(args["output-dir"]))
    max_retries = max(1, args["max-retries"])
    extract = args["extract"]
    parallel = max(1, args["parallel"])
    force = args["force"]

    if !isfile(input_file)
        @error "Input file '$input_file' does not exist"
        exit(1)
    end
    accessions = unique(readdlm(input_file, '\t'; skipstart=1)[:, 1])
    if isempty(accessions)
        @error "Input file '$input_file' is empty"
        exit(1)
    end

    # Filter accessions to download only those that are missing
    archive_dir = joinpath(output_dir, "archives")
    accessions_to_download = [acc for acc in accessions if !isfile(joinpath(archive_dir, "$acc.zip"))]
    skipped = length(accessions) - length(accessions_to_download)
    if skipped > 0
        @info "Skipping $skipped accessions that are already downloaded"
    end
    @info "Processing $(length(accessions_to_download)) accessions with $parallel parallel downloads"

    min_duration = Ref{Float64}(-1) # Tracks fastest download time

    active_tasks = Task[]
    failed_accessions = String[]

    @showprogress for accession in accessions_to_download
        while length(active_tasks) >= parallel
            finished_task = wait_tasks(active_tasks)
            filter!(t -> t !== finished_task, active_tasks)
        end

        task = @async begin
            success = download_accession(accession, output_dir, max_retries, min_duration)
            if !success
                push!(failed_accessions, accession)
            end
        end
        push!(active_tasks, task)
    end

    for task in active_tasks
        wait(task)
    end

    if !isempty(failed_accessions)
        @error "Failed to download $(length(failed_accessions)) accessions: $(join(failed_accessions, ", "))"
    else
        @info "All accessions downloaded successfully"
    end

    if extract
        extract_archives(output_dir; force=force)
    end

    @info "Download process completed"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end