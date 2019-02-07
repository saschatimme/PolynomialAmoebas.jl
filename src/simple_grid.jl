export simple_grid

function simple_grid(F::AbstractFiber, grid, generator::StartValueGenerator, avoid_check=_always_check, B=empty(grid); options=MembershipTestOptions())
    B = empty(grid)
    Fs = [F]
    for _=2:Threads.nthreads()
        push!(Fs, deepcopy(F))
    end

    partition = partition_work(length(grid))
    Threads.@threads for i=1:Threads.nthreads()
        simple_grid_kernel!(B, Fs[i], partition[i], grid, generator, avoid_check, options)
    end
    B
end

function simple_grid_kernel!(B, F::AbstractFiber, range, grid, generator::StartValueGenerator, avoid_check, options)
    for k in range
        if !B[k]
            w = grid[k]
            avoid, val = avoid_check(B, k)
            if !avoid
                B[k] = membershiptest(F, w, generator, options).successfull
            else
                B[k] = val
            end
        end
    end
end

_always_check(B, x) = (false, false)
