
include("Topology.jl")
using Main.Topology
include("Dawg.jl")
using Main.Dawg
include("Search.jl")
using Main.Search
using StaticArrays
using Serialization

function generate_dawgs_for_topology(topology)
    path_types = Dict()
    top_with_dicts = []
    for (path, dict) in topology
        p = (normalize([path])[1], dict)
        if !haskey(path_types, p)
            # TODO embed hash of file
            cache_filename = "cache/" * join([string(x) for x in p[1]], "_") * "_" * replace(dict, "/" => "_") * ".jls"
            if isfile(cache_filename)
                path_types[p] = open(deserialize, cache_filename)
            else
                path_types[p] = create_dawg(p[1], p[2])
                open(f -> serialize(f, path_types[p]), cache_filename, "w")
            end
        end
        push!(top_with_dicts, (path, dict, path_types[p]))
    end
    path_types = path_types|>collect|>sort

    full_dawg = cat([y for (x,y) in path_types]..., dims=1)
    indices = [1]
    for (_,y) in path_types[1:end-1]
        push!(indices, indices[end]+length(y))
    end

    path_types = Dict(zip([x for (x,y) in path_types], indices))

    top_with_dicts = [(path, dict, path_types[(normalize([path])[1], dict)]) for (path, dict, _) in top_with_dicts]

    return (top_with_dicts, full_dawg)
end

function invert_topology(paths)
    symbol_count = length(Set(cat(paths..., dims=1)))
    maximum_interaction = maximum([count(x->i∈x, paths) for i in 0:symbol_count])
    lookup = zeros(UInt8, (maximum_interaction, symbol_count))
    for s in 1:symbol_count
        for (p,path) in enumerate(paths)
            if s ∈ path
                lookup[(findall(x->x==0,lookup[:,s])[1]),s] = p
            end
        end
    end
    return lookup
end

const topology = load_topology(ARGS[1], ARGS[2])
println(topology)

const print_lock = Threads.SpinLock()
function print_result(symbols)
    raw = join([Char(x+'a') for x in symbols])

    formatted = [join([raw[s] for s in path]) for (path, dict) in topology]

    # enable if you don't want duplicate words
    if false && length(Set(formatted)) != length(formatted)
        return
    end

    lock(print_lock) do
        println(join(formatted, " "))
    end
end

# regular run
#topology2 = optimise_labeling(topology)
topology2 = compress_topology_labels(topology)

println("Generating DAWGS")
topology_dawgs,nodes = generate_dawgs_for_topology(topology2)
println("Inverting topology")
paths = invert_topology([path for (path,dict,dawg) in topology_dawgs])
dict_ref = Vector{UInt32}([dawg for (path,dict,dawg) in topology_dawgs])
println("Commencing search")

n = length(dict_ref)
_dict_ref = MVector{n, UInt32}(dict_ref)

n,m = size(paths)
_paths = SMatrix{n, m, UInt8}(paths)

search_multithreaded(nodes, _dict_ref , _paths, print_result)

