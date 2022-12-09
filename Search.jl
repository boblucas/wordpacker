module Search
using StaticArrays

export search_multithreaded, search

# ~10% of CPU time is on this function
@inline function get_child(nodes::Vector{UInt32}, node::UInt32, letter::UInt8)
    # absolute refs
    # return nodes[node+1] + UInt32(2) * (count_ones(nodes[node] & ((UInt32(1) << letter) - UInt32(1))) % UInt32)
    # relative refs
    return node + nodes[node+1] + UInt32(2) * (count_ones(nodes[node] & ((UInt32(1) << letter) - UInt32(1))) % UInt32)
end

# ~20% of CPU time, just on the lookup on the &= line.
@inline function get_mask(indices::SMatrix{MaxCross, Symbols, UInt8}, index, nodes::Vector{UInt32}, roots::MVector{Paths, UInt32}) :: UInt32 where {Paths, Symbols, MaxCross}
	result::UInt32 = 0b111111111111111111111111111
    for j in 1:MaxCross
        d = indices[j,index]
        if j == 1 || d > UInt8(0)
            result &= nodes[roots[d]]
        end
        # ? if d == 0 break is something we can do but is slower unless there is a lot of variance in count
    end
    return result
end

const print_lock = Threads.SpinLock()
function print_result(result::MVector{Symbols, UInt8}) where {Symbols}
    lock(print_lock) do
        Core.println(join([Char(x+'a') for x in result]))
    end
    # TODO: reverse engineer actual words based on path deformation and topology optimilisations
end

function search(nodes::Vector{UInt32}, dawgs::MVector{Paths, UInt32}, paths::SMatrix{MaxCross, Symbols, UInt8}, starts::MVector{Starts, UInt8}, callback) where { Symbols, Paths, Starts, MaxCross }
    @inbounds begin
        stack = MVector{Symbols, UInt8}(zeros(Symbols))
        maskStack = MVector{Symbols, UInt32}(zeros(Symbols))
        maskStack[1] = get_mask(paths, 1,  nodes, dawgs)
        parents = MMatrix{Paths, Symbols, UInt32}(zeros((Paths, Symbols)))

        i = 1

        if maskStack[1] == 0 return end

        for start in starts
            # select letter
            stack[i] = start

            # verify that the letter in question is an option
            if ((maskStack[i] >> start) & 1) == 0 return end

            # move 'down' to next symbol
            s = stack[i]
            for j in 1:MaxCross
                d = paths[j,i]
                if j == 1 || d > UInt8(0)
                    _d = dawgs[d]
                    parents[d,i] = _d 
                    dawgs[d] = get_child(nodes, _d, s)
                end
            end

            i += 1
            stack[i] = UInt32(0)
            maskStack[i] = get_mask(paths, i, nodes, dawgs)
        end

        n_starts = length(starts)
        while i-1 >= n_starts
            # Move to the next valid child
            stack[i] += trailing_zeros(maskStack[i] >> stack[i]) % UInt8

            # If we are not at the last node yet
            if i < Symbols
                # Move down
                s = stack[i]

                for j in 1:MaxCross
                    d = paths[j,i]
                    if j == 1 || d > UInt8(0)
                        _d = dawgs[d]
                        parents[d,i] = _d
                        dawgs[d] = get_child(nodes, _d, s)
                    end
                    # ? if d == 0 break is something we can do but is slower unless there is a lot of variance in count
                end

                i += 1
                stack[i] = UInt32(0)
                maskStack[i] = get_mask(paths, i, nodes, dawgs)

                # move up as long as there are no options
                while (maskStack[i] >> stack[i]) == UInt32(0) && i > 1
                    # Move up
                    i -= 1
    
                    for j in 1:MaxCross
                        d = paths[j,i]
                        if j == 1 || d > UInt8(0)
                            dawgs[d] = parents[d,i]
                        end
                    end
    
                    # And right
                    stack[i] += 1
                end
            else
                # Print result and move right
                callback(stack)
                #print_result(stack)
                stack[i] += 1
                while (maskStack[i] >> stack[i]) == 0 && i > 1
                    # Move up
                    i -= 1
    
                    for j in 1:MaxCross
                        d = paths[j,i]
                        if j == 1 || d > UInt8(0)
                            dawgs[d] = parents[d,i]
                        end
                    end
    
                    # And right
                    stack[i] += 1
                end
            end
        end
    end
end

function search_multithreaded(nodes::Vector{UInt32}, dawgs::MVector{Paths, UInt32}, paths::SMatrix{MaxCross, Symbols, UInt8}, callback) where { Symbols, Paths, MaxCross }
    letters = Base.product(0:26, 0:26, 0:26) |> collect
    Threads.@threads for (x,y,z) in letters
        search(nodes, MVector{Paths, UInt32}(dawgs), paths, MVector{3, UInt8}(x, y, z), callback)
    end
end

end

# example of how to call search directly without any topology file
# useful if you want to do many runs over a space of topologies, or for debugging/testing
if abspath(PROGRAM_FILE) == @__FILE__
    dict_file = open("/run/media/bob/a52a078b-47cd-4847-989d-9eaae30027f0/home/bob/programming/julia/8.bindict", "r")
    dict = read(dict_file)
    dict = reinterpret(UInt32, dict) |> collect
    paths = MVector{36, Vector{UInt8}}(Vector{UInt8}([1]),Vector{UInt8}([1,2]),Vector{UInt8}([1,3]),Vector{UInt8}([1,4]),Vector{UInt8}([1,5]),Vector{UInt8}([1,6]),Vector{UInt8}([1,7]),Vector{UInt8}([1,8]),Vector{UInt8}([2]),Vector{UInt8}([2,3]),Vector{UInt8}([2,4]),Vector{UInt8}([2,5]),Vector{UInt8}([2,6]),Vector{UInt8}([2,7]),Vector{UInt8}([2,8]),Vector{UInt8}([3]),Vector{UInt8}([3,4]),Vector{UInt8}([3,5]),Vector{UInt8}([3,6]),Vector{UInt8}([3,7]),Vector{UInt8}([3,8]),Vector{UInt8}([4]),Vector{UInt8}([4,5]),Vector{UInt8}([4,6]),Vector{UInt8}([4,7]),Vector{UInt8}([4,8]),Vector{UInt8}([5]),Vector{UInt8}([5,6]),Vector{UInt8}([5,7]),Vector{UInt8}([5,8]),Vector{UInt8}([6]),Vector{UInt8}([6,7]),Vector{UInt8}([6,8]),Vector{UInt8}([7]),Vector{UInt8}([7,8]),Vector{UInt8}([8]))
    paths = SMatrix{2, 36, UInt8}(1,0, 1,2, 1,3, 1,4, 1,5, 1,6, 1,7, 1,8, 2,0, 2,3, 2,4, 2,5, 2,6, 2,7, 2,8, 3,0, 3,4, 3,5, 3,6, 3,7, 3,8, 4,0, 4,5, 4,6, 4,7, 4,8, 5,0, 5,6, 5,7, 5,8, 6,0, 6,7, 6,8, 7,0, 7,8, 8,0)
    search_multithreaded(dict, MVector{8, UInt32}(1,1,1,1,1,1,1,1), paths)
    # Profile.print(format=:flat)
    # statprofilehtml()
end
