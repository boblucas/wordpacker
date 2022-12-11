# expresses topology in most effecient form, reducing symmetry and ordering labels

using Base.Iterators
using Combinatorics;

module Topology
export load_topology, save_topology, select_best_symmetry, unpack_symmetry, normalize, optimise_labeling, compress_topology_labels

# topologies are unabelled and paths unordered
# this function returns a canonical representation of a given topology.
# orders paths lexiographically (by their unlabelled representation) and
# normalizes labels of resulting topology
function compress_labels(paths)
	indices = sort(vcat(paths...) |> Set |> collect)
	return paths .|> path -> (path .|> (x -> findfirst(isequal(x), indices)));
end

function normalize(paths)
	return compress_labels(sort(paths, by=(p -> compress_labels([p])[1])));
end


function compress_topology_labels(topology)
	paths = compress_labels([x for (x,_) in topology])
	return [(p, t[2]) for (p, t) in zip(paths, topology)]
end


# while this optimization definitely reduces the size of the tree that needs to searched
# it often causes every path to be a unique topology which means we need more dictionaries
# this cache pressure often costs more performance than the additional tree searches that we save
function optimise_labeling(topology)
	function symbol_complexity(i)
		return count(x->x==i, cat([path for (path,_) in topology]..., dims=1))
	end

	symbol_count = maximum(cat([x[1] for x in topology]..., dims=1))

	for (to, from) in enumerate(sort(1:symbol_count|>collect, by=symbol_complexity, rev=true))
		topology = [([s == from ? to+1024 : s for s in path], dict) for (path, dict) in topology]
	end

	topology = [([s-1024 for s in path], dict) for (path, dict) in topology]

	return topology
end

# equivalent to normalize by stringifies a subset of your topology (convience function).
norm_hash(indices, paths) = join(normalize((x -> paths[x]).(indices[2:end])), ' ')

function adjacent_paths(evaluated_paths, paths, ids)
	if ids == [0]
		return 1:length(paths)
	end
	return filter(i -> (any(j -> j != 0 && length(intersect(paths[j], paths[i])) > 0, ids) && !(i in ids)), 1:length(paths))
end

function bfs(cb, root, children)
	opened, closed = [[root]], []
	while !isempty(opened) || !isempty(closed)
		# first callback on all open nodes
		while !isempty(opened)
			push!(closed, popfirst!(opened))
			cb(closed[end])
		end

		# then get the children of all those nodes
		while !isempty(closed)
			nodes = popfirst!(closed)
			# take each child of nodes and create a new path, one item longer, at the end of `opened`.
			push!(opened, ((c -> push!(copy(nodes), c)).(children(nodes)))...)
		end

		# filter out different paths that non-the-less visit the same nodes
		opened = opened .|> sort |> unique

		# and we have iterated to the next rank
	end
end

function find_repeating_substructure(t)
	groups = Dict{String, Array{Array{Int64}}}()
	evaluated_paths = Set{String}()

	# Converts a list of path indices into t into a key for the 'groups'
	norm_hash(indices) = join(normalize((x -> t[x]).(indices[2:end])), ' ')

	# BFS over all combinations of paths such that
	# - No new path is added when a subset of paths is alone in its group (has a unique key).
	# - A given subset is a single component where an edge is defined as sharing at least one label.
	group_get(subset) = get!(groups, norm_hash(subset), [])
	register(subset) = push!(group_get(subset), subset)
	do_recurse(bt) = length(bt) == 1 || length(group_get(bt)) > 1

	bfs(register, 0, x -> if do_recurse(x) adjacent_paths(evaluated_paths, t, x) else [] end)

	return collect(filter(x -> length(x) > 1, collect(values(groups))))
end

# Array<(Array<Int>, string)>
function load_topology(filename, fallbackdict)
	topology = []
	for line in filter(x -> x != "" && x[1] != '#', readlines(filename))
		line = split(line, ':')
		path = map(x -> parse(Int64, x), split(line[1], ','))
		if length(line) >= 2
			push!(topology, (path, line[2]))
		else
			push!(topology, (path, fallbackdict))
		end
	end

	if any(x -> 0 in x[1], topology)
		topology = [([x+1 for x in path], d) for (path, d) in topology]
	end

	return topology
end

function save_topology(filename, topology)
	println(topology)
	open(filename, "w") do f
		print(f, join(map(path -> join(map(string, path[1]), ",") * ":" * path[2], topology), "\n"))
	end
end

# filter out all non-path matching words and return Array<string>
function load_dictionary(path, filename)
	function matches_path(word)
		if length(word) != length(path) return false end
		return all(i -> word[i] == word[findfirst(x -> x == path[i], path[1:i])], 1:length(path))
	end
	return filter(matches_path, readlines(filename))
end

# (dictionary filename, letterposition, letter) -> letter frequency at letterposition
function precompute_frequence_table(topology)
	table = Dict()
	cache = Dict()
	for path in topology
		path_norm = normalize([path])[1]
		if haskey(cache, path_norm)
			table[path] = cache[path_norm]
			continue
		end

		d = load_dictionary(path[1], path[2])
		table[path] = []

		for i in 1:length(path[1])
			push!(table[path], zeros(26))
		end

		for word in d
			for i in 1:length(word)
				table[path][i][findfirst(x -> x == word[i], "abcdefghijklmnopqrstuvwxyz{")] += 1
			end
		end

		cache[path_norm] = table[path]
	end
	return table
end

# effectively solves each label independently, and combines.
# would be accurate if there were no relationships between parts of a word.
# but there are, and so tends to overestimate.
function solution_count(topology)
	freq = precompute_frequence_table(topology)
	N = maximum(map(x -> maximum(x[1]), topology))
	dictionarysizes = map(x -> sum(freq[x][1]), topology) |> collect

	aggregate = []
	# iterate all symbols in topology
	for l in 1:N-1
		total = ones(27)
		# consider all paths that contain this symbol
		for path in filter(path -> l in path[1], topology)
			# then find where in the path/word this symbol existis
			for i in filter(i -> path[1][i] == l, 1:length(path[1]))
				# and get the amount of words that follow that pattern
				# one value for each of the letters
				total = total .* freq[path][i]
			end
		end

		# for each letter that can be on symbol `l` we know have a frequency by simply taking the product
		# of all options of all paths that interact with that symbol
		# then we take a look at how many options we would have WITHOUT that symbol equivalance constraint
		# which is simply the product of the dictionaries of all relevant paths
		no_constraint_total = reduce(*, map(x -> dictionarysizes[x[1]], filter(path -> l in path[2][1], enumerate(topology)|>collect)))
		# the ratio gives us a measurement of how stringent this constraint is.
		println("Constraint of $l is $(sum(total)/no_constraint_total)")
		push!(aggregate, sum(total)/no_constraint_total)
	end

	# alternatively we could iterate all PAIRS of symbols
	# get all paths that contain at least one of those
	# check for each PAIR of letters what the frequency of fitting options is in each path when we fix those symbols
	# then again divide over the unconstrained case

	# numeric stability is bothersome, doing `return prod(dictionarysizes) * prod(aggregate)`
	result = 1.0
	while length(dictionarysizes) > 0 || length(aggregate) > 0
		if (result < 1.0 || length(aggregate) == 0) && length(dictionarysizes) > 0
			result *= pop!(dictionarysizes)
		end
		if (result >= 1.0 || length(dictionarysizes) == 0) && length(aggregate) > 0
			result *= pop!(aggregate)
		end
	end

	return result
end

function symmetry_score(symmetry, paths)
	# how often can a precomputed dictionary be applied
	folding = length(symmetry)
	return sqrt(folding) * sum(map(length, map(i -> paths[i], filter(x -> x != 0, symmetry[1]) )))
end

function select_best_symmetry(t)
	for subset in sort(find_repeating_substructure(map(x -> x[1], t)), by=x->symmetry_score(x, map(x -> x[1], t)), rev=true) 
		println(subset)
		example = map(i -> t[i], filter(i -> i > 0,subset[1]))
		if length(subset) >= 0 && length(example) > 1 && solution_count(example) < 1e10
			return subset
		end
	end
	return nothing
end

# dict filename => topology
function unpack_symmetry(t, symm)
	output = []
	t2 = deepcopy(t)
	blueprint = map(i -> t[i], filter(i -> i > 0, symm[1]))
	
	for x in symm
		example = map(i -> t[i], filter(i -> i > 0, x)) |> collect
		map(i -> t2[i] = [], filter(i -> i > 0, x)) |> collect
		push!(t2, (Iterators.flatten(map(x -> x[1], example)) |> collect, "symmetry.output"))
	end
	main = filter(x -> length(x) > 0, t2)|>collect

	return Dict("sym"=>optimise_labeling(blueprint), "main"=>optimise_labeling(main))
end



if abspath(PROGRAM_FILE) == @__FILE__
	# compute expected amount of solutions for given topology and dictionary
	t = load_topology(ARGS[1], ARGS[2])
	expected_total = solution_count(t)
	println(expected_total)

	# generate and store optimale topological representation
	if false
		s = select_best_symmetry(t)
		for topology in unpack_symmetry(t, s)
			save_topology(topology[1], topology[2])
		end
	end

	# show details of symmetric structure within topology
	if false
		tp = map(x -> x[1], t)
		for subset in sort(find_repeating_substructure(tp), by=x->symmetry_score(x, tp), rev=true) 
			if length(subset) >= 0
				example = map(i -> t[i], filter(i -> i > 0,subset[1]))
		
				println("--------", symmetry_score(subset, tp), " ", solution_count(example))
				for value in subset
					println(value[2:end])
				end
			end
		end
	end
end

end
