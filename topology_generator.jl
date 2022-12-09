# generates common topologies in all dimensions

using Base.Iterators

module Generators
export square, rect

function square(n)
    t(r) = r * (2n - r + 1) ÷ 2
    return (r -> vcat((i -> t(i)+r-i).(0:r-1), t(r):t(r)+n-r-1)).(0:n-1)
end

function rect(w, h)
    return vcat((i -> i*w:i*w+w-1).(0:h-1), (i -> i:w:w*h-1).(0:w-1))
end

end

using Main.Generators

function write_topologies()
    fnames = (s -> "$(s)").(names(Generators)[2:end])
    functions = eval.(names(Generators)[2:end])

    function write_topology_type(name, f)
        argn(f) = findfirst((i -> hasmethod(f, fill(Any, i))).(0:10)) - 1
        persist(x) = write(open("$name/$(join(x, 'x'))", "w"), join(join.(f(x...), ','), '\n'))
        mkpath(name)
        [product(fill(3:14-2argn(f), argn(f))...)...] .|> persist
    end

    (x -> write_topology_type.(x[1], x[2])).(zip(fnames, functions))
end

write_topologies()
