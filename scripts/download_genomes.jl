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

function download_accession(accession::AbstractString, output_dir::AbstractString, max_retries::Int)
    archive_dir = joinpath(output_dir, "archives")
    mkpath(archive_dir)
    archive_path = joinpath(archive_dir, "$accession.zip")
    archive_part = "$archive_path.part"

    attempts = 0
    success = false

    while attempts < max_retries && !success
        attempts += 1
        rm(archive_part, force=true)

        try
            cmd = pipeline(`datasets download genome accession $accession --include genome,gff3 --filename $archive_part`, stdout=devnull, stderr=devnull)
            run(cmd)
            mv(archive_part, archive_path, force=true)
            success = true
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
    prog = Progress(length(zips); desc="Extracting archives")
    @threads for zip in zips
        accession = splitext(basename(zip))[1]
        extraction_dir = joinpath(genome_dir, accession)
        mkpath(extraction_dir)
        target_fna = joinpath(extraction_dir, "$accession.fna")
        target_gff = joinpath(extraction_dir, "$accession.gff")
        target_json = joinpath(extraction_dir, "$accession.json")

        if all(isfile, [target_fna, target_gff, target_json]) && !force
            next!(prog)
            continue
        end
        
        temp_dir = mktempdir()
        try
            run(pipeline(`unzip -o $zip -d $temp_dir`, devnull))
            
            data_dir = joinpath(temp_dir, "ncbi_dataset", "data", accession)
            if isdir(data_dir)
                # Handle .fna file
                fna_files = filter(f -> endswith(f, "_genomic.fna"), readdir(data_dir, join=true))
                if !isempty(fna_files)
                    mv(first(fna_files), target_fna, force=true)
                else
                    @warn "No _genomic.fna file found in $data_dir for $accession"
                end
                
                # Handle .gff file
                gff_files = filter(f -> endswith(f, "genomic.gff"), readdir(data_dir, join=true))
                if !isempty(gff_files)
                    mv(first(gff_files), target_gff, force=true)
                else
                    @warn "GFF file missing in $data_dir for $accession"
                end
            else
                @warn "Data directory not found for $accession"
            end
            
            # Handle .json file
            json_files = filter(f -> endswith(f, "assembly_data_report.jsonl"), readdir(joinpath(temp_dir, "ncbi_dataset", "data"), join=true))
            if !isempty(json_files)
                mv(first(json_files), target_json, force=true)
            else
                @warn "Assembly report JSON missing for $accession"
            end

            # Check if all required files are present; remove directory if not
            if !all(isfile, [target_fna, target_gff, target_json])
                @warn "Incomplete extraction for $accession: missing one or more files. Removing directory."
                rm(extraction_dir, recursive=true, force=true)
            end
        catch e
            @error "Failed to process $zip for $accession: $e (Thread $(threadid()))"
            rm(extraction_dir, recursive=true, force=true)  # Clean up on error
        finally
            rm(temp_dir, recursive=true, force=true)
        end
        next!(prog)
    end
    clear_last_lines(1)
    @info "Extraction and reorganization completed"
end

function verify_extracted_directories(output_dir::String)
    genome_dir = joinpath(output_dir, "genomes")
    
    if !isdir(genome_dir)
        @warn "Genomes directory not found: $genome_dir. Skipping verification."
        return
    end
    
    subdirs = filter(d -> isdir(joinpath(genome_dir, d)), readdir(genome_dir))
    
    if isempty(subdirs)
        @info "No extracted directories found in $genome_dir. Verification complete."
        return
    end
    
    @showprogress desc="Final verification" for subdir in subdirs
        extraction_dir = joinpath(genome_dir, subdir)
        target_fna = joinpath(extraction_dir, "$subdir.fna")
        target_gff = joinpath(extraction_dir, "$subdir.gff")
        target_json = joinpath(extraction_dir, "$subdir.json")
        
        files = readdir(extraction_dir)
        
        expected_files = [target_fna, target_gff, target_json]
        if !all(isfile, expected_files) || length(files) != 3
            @warn "Directory $subdir is invalid: expected exactly 3 files (.fna, .gff, .json). Found: $(length(files)). Removing directory."
            rm(extraction_dir, recursive=true, force=true)
        end
    end
    clear_last_lines(1)
    @info "Final verification of extracted directories completed."
end

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

    archive_dir = joinpath(output_dir, "archives")
    accessions_to_download = [acc for acc in accessions if !isfile(joinpath(archive_dir, "$acc.zip"))]
    skipped = length(accessions) - length(accessions_to_download)
    if skipped > 0
        @info "Skipping $skipped accessions that are already downloaded"
    end
    @info "Processing $(length(accessions_to_download)) accessions with $parallel parallel downloads"

    sem = Channel{Bool}(parallel)
    for i in 1:parallel
        put!(sem, true)
    end

    tasks = Task[]
    failed_accessions = String[]
    failed_lock = ReentrantLock()

    @showprogress desc="Downloading accessions" barlen=40 showspeed=true for accession in accessions_to_download
        take!(sem)
        task = @async begin
            try
                success = download_accession(accession, output_dir, max_retries)
                if !success
                    lock(failed_lock)
                    push!(failed_accessions, accession)
                    unlock(failed_lock)
                end
            finally
                put!(sem, true)
            end
        end
        push!(tasks, task)
    end

    for task in tasks
        wait(task)
    end

    if !isempty(failed_accessions)
        @error "Failed to download $(length(failed_accessions)) accessions: $(join(failed_accessions, ", "))"
    else
        @info "All accessions downloaded successfully"
    end

    if extract
        extract_archives(output_dir; force=force)
        verify_extracted_directories(output_dir)
    end

    @info "Download process completed"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end