
# creates and stores all dawgs for a given topology (or just a single dictionary)
# in the most basic form we can just express a dict as a tree
# 2nd we can compress that tree's representation so that all children of a given node are contigious in memory and valueless, requiring only a bitset of the alphabet and a memory location.
# if in this scheme 2 nodes have identical child-trees we compress by removing one of them.

# TODO combine multiple dawgs into one array that is more compressed than just concat
# TODO cross rank compression
# TODO chained overlap compression

#
# load dictionary
#

module Dawg
using IterTools
using DataStructures

export create_dawg, load_words

# we should do a readlines and do all the various checks in parallel
function load_words(filename, path)
    words = []
    for line in filter(x -> x != "" && x[1] != "#", readlines(filename))
        word = split(line, ' ')[end]

        if length(word) == length(path) && all(c -> 0x61 <= UInt8(c) <= 0x7B, word)
            wt = sort(zip(path, word)|>collect)

            if all(x -> x[1][1] != x[2][1] || x[1][2] == x[2][2], partition(wt, 2, 1))
                w2 = join(map(x -> x[2], unique(x -> x[1], wt)))
                #println(word, " = ", w2)
                push!(words, w2)

            end

        end
    end
    sort!(words)
    unique!(words)
    return words
end

#
# Put words in tree
#
mutable struct Node
    children::Dict{Char, Union{Node, Nothing}}
end

function add(root, word)
    if word == "" return end
    for c in word[1:end-1]
        if !haskey(root.children, c)
            root.children[c] = Node(Dict())
        end
        root = root.children[c]
    end
    root.children[word[end]] = nothing
end

function get_children_mask(node::Node)
    if length(node.children) == 0
        return 0
    else
        return reduce(|, map(x -> 1 << (x - 'a'), keys(node.children)|>collect ))
    end
end

mutable struct CNode
    letters::UInt32
    children::UInt32
end

function treeToArray(wordlist::Vector{String})
    # first we store the unique values in the first column
    # then we iterate those values and store all unique values in the second column within each group
    # do so breath first until the last column
end


function treeToArray(in::Node, N)
    compacted = []
    for i in 1:N push!(compacted, []) end
    rank_count = zeros(Int64, N)

    q = Queue{Tuple{Node, Int64}}()
    enqueue!(q, (in, 0))
    while length(q) > 0
        next = dequeue!(q)
        node = next[1]
        rank = next[2]

        push!(compacted[rank+1], CNode(get_children_mask(node), length(node.children) > 0 && rank < N-1 ? rank_count[rank+1]+1 : 0))

        for child in map(x->x[2], sort(node.children|>collect, by=x->x[1]))
            if isnothing(child)
                continue
            end
            rank_count[rank+1] += 1
            enqueue!(q, (child, rank + 1))
        end
    end
    return compacted
end

using Combinatorics

function contains(x, y)
    for i in 1:length(x)-length(y)
        if x[i:i+length(y)-1] == y
            return i
        end
    end
    return nothing
end

function compress_rank(ranks, r)

    println("Collecting node-sequences")

    # all (indexes, node count) pairs that may be read from ranks[r]
    required_seqs = Set(map(node -> (node.children, count_ones(node.letters)), ranks[r-1]))
    # all sequences that are subsets of other sequences and so do not require idendendent representation
    # maps (index, node count) -> (index, offset, node count)
    seq_subsets = Dict()

    println("Creating flat representation")
    # flatten to completely simple array
    # - this can be done in parallel, since the exact location 
    simple = Vector{UInt32}()
    for node in ranks[r]
        push!(simple, node.letters)
        push!(simple, node.children)
    end

    # hash all relevant sub-sequences from that array and group by hash
    # perhaps we can split this into 2 tasks that are better parallelizable/optimised
    # - convert each item in required_seqs all of its k->v pairs, totally independent
    # - group that giant list by key for insertion into the hashmap (sorting is highly optimized)
    # this also means you can prevent hashing groups of size one from the get-go
    println("Finding equivalent sequences in ", length(required_seqs), " node sequences")
    subset_hash = Dict{Vector{UInt32}, Vector{NTuple{4, UInt32}}}()
    for (index, count) in required_seqs
        for index2 in index:index+count-1
            for count2 in 1:count-(index2-index)
                k = simple[(index2-1)*2+1:(index2-1)*2+count2*2]
                if !haskey(subset_hash, k) subset_hash[k] = [] end
                push!(subset_hash[k], (index, count, index2, count2))
            end
        end
    end
    for k in keys(subset_hash)
        if length(subset_hash[k]) <= 1
            delete!(subset_hash, k)
        end
    end
    
    println("Generating mappings for found sequences")
    # convert groups to useful format for further processing
    for group in values(subset_hash)
        whole = filter(x -> x[1:2] == x[3:4], group) |> collect
        if length(whole) > 0
            partial = filter(x -> x[1:2] != x[3:4], group) |> collect
            (ri,rn,ri2,_) = length(partial) > 0 ? partial[1] : whole[1]

            for (i1,n1,_,_) in whole[(length(partial) > 0 ? 1 : 2):end]
                seq_subsets[(i1,n1)] = ((ri, rn), ri2-ri, n1)
            end
        end
    end

    println("Removing equivalent sequences")
    # that means we only have to write these
    to_write_set = setdiff(required_seqs, keys(seq_subsets)|>Set) 
    to_write = to_write_set |> collect

    final_positions = Dict()

    println("Find partially overlapping sequences")
    # next we try to find partial overlap and resolve maximum overlap first
    # obviously each group can be evaluated independently, so parallelize
    overlapping = Vector{Tuple{Int64, Tuple{Int64, Int64}, Tuple{Int64, Int64}}}()
    for group in values(subset_hash)
        n = group[1][4]
        group = filter(x -> x[1:2] in to_write, group) |> collect
        endings   = filter(x -> sum(x[1:2]) == sum(x[3:4]) && x[1] != x[3], group) |> collect
        startings = filter(x -> x[1] == x[3] && x[2] != x[4], group) |> collect
        for e in endings
            for s in startings
                if e[1:2] != s[1:2]
                    push!(overlapping, (n, e[1:2], s[1:2]))
                end
            end
        end
    end

    sort!(overlapping, by=x->-x[1])
    unique!(x -> x[2:end], overlapping)

    result = []
    cuts = []
    while length(overlapping) > 0
        (o,a,b) = overlapping[1]
        deleteat!(overlapping, findall(x->x[2]==a||x[3]==a||x[2]==b||x[3]==b, overlapping))
        deleteat!(to_write, findall(x->x==a||x==b, to_write))
        final_positions[a] = length(result)+1
        final_positions[b] = length(result)+1+a[2]-o
        
        push!(cuts, length(result))
        result = [result; ranks[r][a[1]:a[1]+a[2]-1]; ranks[r][b[1]+o:b[1]+b[2]-1]]
    end
    
    println(length(to_write), " items length")
    for a in to_write
        final_positions[a] = length(result)+1
        result = [result; ranks[r][a[1]:a[1]+a[2]-1]]
    end

    println("rewrite graph")

    # convert seq_subsets to final positions
    while length(seq_subsets) > 0
        for (k,v) in seq_subsets|>collect
            if haskey(final_positions, v[1])
                delete!(seq_subsets, k)
                final_positions[k] = final_positions[v[1]] + v[2]
            end
        end
    end

    # and then rewrite results to reference final positions
    for node in ranks[r-1]
        node.children = 1+2(final_positions[(node.children, count_ones(node.letters))]-1)
    end

    # halve nodes that have no further children
    if r == length(ranks) && length(cuts) > 0
        best_cut = maximum(x -> (abs(length(result)÷2-x), x), cuts)[2]
        left = result[1:best_cut]
        right = result[best_cut+1:end]

        while length(left) < length(right) push!(left, CNode(0,0)) end
        while length(right) < length(left) push!(right, CNode(0,0)) end

        for node in ranks[r-1]
            if (node.children-1)÷2+1 > best_cut
                node.children = node.children - best_cut + 1
            end
        end

        result = [[CNode(left[i].letters, right[i].letters) for i in 1:length(left)]; [CNode(0,0)]]
    end

    return result
end

function recombine_ranks!(ranks)
    for rank in ranks[1:end-1]
        for j in 1:2:length(rank)
            rank[j+1] = length(rank) - j + rank[j+1]
        end
    end
    return cat(ranks..., dims=1)
end

function compress_all_ranks!(ranks)
    for i in length(ranks):-1:2
        precomp = length(ranks[i])
        ranks[i] = compress_rank(ranks, i)
        println(length(ranks[i]), "/", precomp)
    end
end

bits(x) = (i for i in 0:27 if (x>>i & 1) > 0)
function read_all(n, out, tree, node, word = "")
    if length(word) == n
        push!(out, word)
        return
    end
    
    for x in enumerate(bits(tree[node]))
        w2 = word * string('a' + x[2])
        child = node + tree[node+1] + 2(x[1]-1)
        read_all(n, out, tree, child, w2)
    end
end

function verify(words, tree, n)
    loaded_words = []
    read_all(n, loaded_words, tree, 1)

    words = Set(words)
    loaded_words = Set(loaded_words)

    too_many = setdiff(loaded_words, words)
    too_few = setdiff(words, loaded_words)

    if length(too_many) > 0
        println("TOO MANY ", too_many)
    end
    if length(too_few) > 0
        println("TOO FEW ", too_few)
    end
    return words == loaded_words
end

function create_dawg(path_form, dictionary_filename)

    # just a list of words (but mangled to fit path_form)
    println("Loading words from dict")
    words = load_words(dictionary_filename, path_form)

    println("Sorting...")
    sort!(words)

    println("Building trie")
    # NOTE: ideally we never build this trie and instead compile to a flat
    # representation in one go. treeToArray does a breath first search
    # so theoratically just recursively reading next column of the 'matrix' of words
    # (filtering on the all the previous columns being equal) is just fast
    # but doesn't require making a whole trie. and make flatting faster too to boot!

    # a regular trie of those words
    root = Node(Dict())
    for word in words
        if word == "" continue end

        active = root
        for c in word[1:end-1]
            if !haskey(active.children, c)
                active.children[c] = Node(Dict())
            end
            active = active.children[c]
        end
        active.children[word[end]] = nothing
    end

    println("Flatten to array format")

    # storing that trie in contigious arrays (1 for each rank) of nodes
    ranks = treeToArray(root, length(Set(path_form)))

    println("Compress")

    # removing duplicate nodes and rewriting references
    compress_all_ranks!(ranks)

    # flatting to a simple integer array (letter mask, child index)
    simple_ranks = Vector{Vector{UInt32}}()
    for rank in ranks
        simple_vector = Vector{UInt32}()
        sizehint!(simple_vector, 2length(rank))
        for node in rank
            push!(simple_vector, node.letters)
            push!(simple_vector, node.children)
        end
        push!(simple_ranks, simple_vector)
    end

    # combining all ranks into a single and rewriting references to be relative
    final_graph = recombine_ranks!(simple_ranks)
    #println(final_graph)
    #println(join([(string(x,base=2) * "." * string(Integer(y))) for (x,y) in partition(final_graph, 2, 2)], " "))

    # verify that after all that the words in the flat-tree are the same as in the loaded dict
    if verify(words, final_graph, length(Set(path_form)))
        println("ALL CORRECT")
    end
    
    return final_graph
end


#open("test.bindict", "w") do file
#    for child in final_graph
#        write(file, child.letters)
#        write(file, child.children)
#    end
#end

end

#using Main.Dawg
#PATH = map(x->parse(Int64, x), split(ARGS[2], ','))|>collect
#create_dawg(PATH, ARGS[1])