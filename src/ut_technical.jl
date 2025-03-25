function clear_last_lines(n::Int)
    for _ in 1:n
        print("\e[1A")
        print("\e[2K")
    end
    flush(stdout)
end

function check_dependencies(execs)
    for cmd in execs
        try
            run(pipeline(`which $cmd`, devnull))
        catch
            @error "Required command '$cmd' not found in PATH"
            exit(1)
        end
    end
end

function wait_tasks(tasks::Vector{Task})
    while true
        for t in tasks
            if istaskdone(t)
                return t
            end
        end
        sleep(0.5)
    end
end