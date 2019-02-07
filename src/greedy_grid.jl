function greedy_grid(F::AbstractFiber, grid, generator::StartValueGenerator,
    initial_queue,
    bitmap=empty(grid), queued=copy(bitmap.data); callback=_do_nothing, options=MembershipTestOptions())
    # queued = copy(bitmap.data)
    queue = copy(initial_queue)

    if Threads.nthreads() > 1
        Fs = [F]
        for _=2:Threads.nthreads()
            push!(Fs, deepcopy(F))
        end
        partition = partition_work(length(queue))
        Threads.@threads for tid = 1:Threads.nthreads()
            t_queue = queue[partition[tid]]
            G = Fs[tid]
            greedy_grid_kernel!(bitmap, queued, t_queue, G, grid, generator, options, callback)
        end
    else
        greedy_grid_kernel!(bitmap, queued, queue, F, grid, generator, options, callback)
    end

    bitmap
end

function greedy_grid_kernel!(bitmap, queued, queue, F::AbstractFiber, grid, generator::StartValueGenerator, options, callback)
    while !isempty(queue)
        k = pop!(queue)
        if !bitmap[k]
            w = grid[k]
            result = membershiptest(F, w, generator, options)
            if result.successfull
                bitmap[k] = true
                callback(result, k, bitmap)
                addneighbours!(queue, queued, k)
            end
        end
    end
end

_do_nothing(result, k, bitmap) = nothing

function addneighbours!(queue, queued::BitMatrix, ij)
    i, j = Tuple(ij)
    i0 = max(i - 1, 1)
    i1 = min(i + 1, size(queued, 1))

    j0 = max(j - 1, 1)
    j1 = min(j + 1, size(queued, 2))
    for jj = j0:j1, ii = i0:i1
        if !queued[ii, jj]
            queued[ii, jj] = true
            push!(queue,  CartesianIndex(ii, jj))
        end
    end
end

function addneighbours!(queue, queued::BitArray{3}, ijk)
    i, j, k = Tuple(ijk)
    i0 = max(i - 1, 1)
    i1 = min(i + 1, size(queued, 1))

    j0 = max(j - 1, 1)
    j1 = min(j + 1, size(queued, 2))

    k0 = max(k - 1, 1)
    k1 = min(k + 1, size(queued, 3))
    for kk=k0:k1, jj = j0:j1, ii = i0:i1
        if !queued[ii, jj, kk]
            queued[ii, jj, kk] = true
            push!(queue, CartesianIndex(ii, jj, kk))
        end
    end
end


function initial_neighbour_queue(B::BitArray{N}) where N
    queued = copy(B)
    initial_queue = CartesianIndex{N}[]
    for (ind, val) in pairs(B)
        if val
            addneighbours!(initial_queue, queued, ind)
        end
    end
    initial_queue, queued
end

# compared to `greedy_grid ` does this uses previous neighboring startvalues.
function greedy_grid_memorized(F::AbstractFiber, grid, generator::StartValueGenerator,
    initial_queue::Vector,
    startvalues::Array,
    bitmap=empty(grid), queued=copy(bitmap.data); options=MembershipTestOptions())
    # queued = copy(bitmap.data)

    Fs = [F]
    for _=2:Threads.nthreads()
        push!(Fs, deepcopy(F))
    end

    partition = partition_work(length(initial_queue))
    Threads.@threads for tid = 1:Threads.nthreads()
        t_queue = initial_queue[partition[tid]]
        G = Fs[tid]
        greedy_grid_memorized_kernel!(bitmap, queued, t_queue, startvalues, G, grid, generator, options)
    end

    bitmap
end
function greedy_grid_memorized_kernel!(bitmap, queued, queue, startvalues, F::AbstractFiber, grid, generator::StartValueGenerator, options)
    tmp_neighbor_startvalues = Vector{eltype(startvalues)}()
    while !isempty(queue)
        k = pop!(queue)
        # @show length(queue)
        if !bitmap[k]
            w = grid[k]
            neighbor_startvalues!(tmp_neighbor_startvalues, F, startvalues, bitmap, k)
            result = membershiptest(F, w, generator, tmp_neighbor_startvalues, options)
            if result.successfull
                startvalues[k] = result.solution
                bitmap[k] = true
                addneighbours!(queue, queued, bitmap, k)
            end
        end
    end
end

function addneighbours!(queue, queued::BitMatrix, bitmap::Bitmap2D, k)
    i, j = Tuple(k)
    i0 = max(i - 1, 1)
    i1 = min(i + 1, size(queued, 1))

    j0 = max(j - 1, 1)
    j1 = min(j + 1, size(queued, 2))
    for ii = i0:i1, jj = j0:j1
        if !queued[ii, jj] || !bitmap[ii, jj]
            queued[ii, jj] = true
            push!(queue, CartesianIndex(ii, jj))
        end
    end
end

function addneighbours!(queue, queued::BitArray{3}, bitmap::Bitmap3D, l)
    i, j, k = Tuple(l)
    i0 = max(i - 1, 1)
    i1 = min(i + 1, size(queued, 1))

    j0 = max(j - 1, 1)
    j1 = min(j + 1, size(queued, 2))

    k0 = max(k - 1, 1)
    k1 = min(k + 1, size(queued, 3))
    for ii = i0:i1, jj = j0:j1, kk=k0:k1
        if !queued[ii, jj, kk] || !bitmap[ii, jj, kk]
            queued[ii, jj, kk] = true
            push!(queue, CartesianIndex(ii, jj, kk))
        end
    end
end
function neighbor_startvalues!(out, F, startvalues, B::Bitmap2D, k)
    empty!(out)
    i, j = Tuple(k)
    offset = 1
    for jj in j-offset:j+offset, ii in i-offset:i+offset
        if safe_getindex(B, ii, jj) && !isnan(startvalues[ii, jj])
            predict_startvalue!(out, F, startvalues[ii, jj], B.grid[ii, jj], B.grid[k])
        end
    end
end

function neighbor_startvalues!(out, F, startvalues, B::Bitmap3D, index)
    empty!(out)
    i, j, k = Tuple(index)
    offset = 1
    for kk in k-offset:k+offset, jj in j-offset:j+offset, ii in i-offset:i+offset
        if safe_getindex(B, ii, jj, kk) && !isnan(startvalues[ii, jj, kk])
            predict_startvalue!(out, F, startvalues[ii, jj, kk], B.grid[ii, jj, kk], B.grid[index])
        end
    end
    out
end

function predict_startvalue!(out, F, x_known, y_known, y_target)
    push!(out, x_known)
end
function predict_startvalue!(out, F::ImaginaryFiber2D, x_known, y_known, y_target)
    xy = SVector(x_known[1], x_known[2], y_known[1], y_known[2])
    U = F.U
    @inbounds begin
        ∇f_re = SP.gradient(F.f_re, xy)
        for i=1:4
            U[1, i] = ∇f_re[i]
        end
        U[2, 3] = U[1, 1]
        U[2, 4] = U[1, 2]
        U[2, 1] = -U[1, 3]
        U[2, 2] = -U[1, 4]
    end

    ∂xy = SMatrix{2,2}(U[1,1], U[2,1], U[1, 2], U[2, 2])
    ∂t = SMatrix{2,2}(U[1,3], U[2,3], U[1, 4], U[2, 4]) *
        SVector(y_target[1] - y_known[1], y_target[2] - y_known[2])

    Δx, ok = mysolve(∂xy, ∂t)
    if ok
        push!(out, x_known - Δx)
        # push!(out, x_known - 0.5Δx)
        push!(out, x_known)

    else
        @show "not ok"
        push!(out, x_known)
    end
end
