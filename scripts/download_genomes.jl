#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
# const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
# include(UTILS_FILE)

using ArgParse
using Dates
using Logging

# Parse command-line arguments
function parse_commandline()
    s = ArgParseSettings(prog="Genome Downloader", 
                         description="Download genomic data for each accession separately.")
    @add_arg_table! s begin
        "--input-file", "-i"
            help = "Path to accessions file"
            required = true
        "--output-dir", "-o"
            help = "Output directory"
            required = true
        "--max-retries", "-r"
            help = "Maximum download attempts per accession"
            arg_type = Int
            default = 3
        "--extract", "-e"
            help = "Extract downloaded archives"
            action = :store_true
        "--parallel", "-p"
            help = "Number of simultaneous downloads"
            arg_type = Int
            default = 1
    end
    return parse_args(s)
end

# Check external dependencies
function check_dependencies()
    for cmd in ["datasets", "unzip"]
        try
            run(pipeline(`which $cmd`, "/dev/null"))
        catch
            @error "Required command '$cmd' not found in PATH"
            exit(1)
        end
    end
    @info "All dependencies found"
end

# Download a single accession
function download_accession(accession::AbstractString, output_dir::AbstractString, max_retries::Int, min_duration::Ref{Float64})
    archive_dir = joinpath(output_dir, "archives")
    mkpath(archive_dir)
    archive_path = joinpath(archive_dir, "$accession.zip")
    archive_part = "$archive_path.part"
    
    # Check if the archive already exists
    if isfile(archive_path)
        @info "Accession $accession: Archive already exists, skipping download"
        return true
    end
    
    attempts = 0
    success = false
    duration = 0.0
    
    while attempts < max_retries && !success
        attempts += 1
        rm(archive_part, force=true)
        
        @info "Attempt $attempts of $max_retries for accession $accession"
        start_time = now()
        try
            cmd = pipeline(`datasets download genome accession $accession --include genome,gff3 --filename $archive_part`, "/dev/null")
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

function extract_archives(output_dir::String)
    archive_dir = joinpath(output_dir, "archives")
    genome_dir = joinpath(output_dir, "genomes")
    mkpath(genome_dir)
    
    zips = filter(f -> endswith(f, ".zip"), readdir(archive_dir, join=true))
    
    Threads.@threads for zip in zips
        @info "Extracting $zip (Thread $(Threads.threadid()))"
        try
            run(pipeline(`unzip -o $zip -d $genome_dir`, "/dev/null"))
        catch e
            @error "Failed to extract $zip: $e (Thread $(Threads.threadid()))"
        end
    end
    @info "Extraction completed"
end

# Wait for any task to finish
function wait_any(tasks::Vector{Task})
    while true
        for t in tasks
            if istaskdone(t)
                return t
            end
        end
        sleep(0.1)  # Avoid busy-waiting
    end
end

# Main function
function main()
    global_logger(ConsoleLogger(stdout, Logging.Info))
    args = parse_commandline()
    check_dependencies()
    
    input_file = abspath(args["input-file"])
    output_dir = abspath(args["output-dir"])
    max_retries = max(1, args["max-retries"])
    extract = args["extract"]
    parallel = max(1, args["parallel"])
    
    if !isfile(input_file)
        @error "Input file '$input_file' does not exist"
        exit(1)
    end
    accessions = unique!(filter(!isempty, strip.(readlines(input_file))))
    if isempty(accessions)
        @error "Input file '$input_file' is empty"
        exit(1)
    end
    
    @info "Processing $(length(accessions)) accessions with $parallel parallel downloads"
    min_duration = Ref{Float64}(-1)  # Tracks fastest download time
    
    active_tasks = Task[]
    failed_accessions = String[]
    
    for accession in accessions
        while length(active_tasks) >= parallel
            finished_task = wait_any(active_tasks)
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
        extract_archives(output_dir)
    end
    
    @info "Download process completed"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end