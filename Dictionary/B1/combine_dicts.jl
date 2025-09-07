using MAT, Glob, LinearAlgebra

function main()
    dict_folder = raw"C:\Users\LucyA\MSC_PROJECT\Dictionary\dict\b1\8mm_5001/"
    outfile = joinpath(dict_folder, "5001_combined.mat")

    files = sort(glob("*.mat", dict_folder))

    dicts = Matrix{ComplexF32}[]
    idxs  = Matrix{Float32}[]

    for fp in files
        S = matopen(fp) do f
            (; dict0 = read(f, "dict0"), idx = read(f, "idx"))
        end
        d0  = ComplexF32.(S.dict0)
        idx = Float32.(S.idx)

        if size(d0,1) < size(d0,2)
            d0 = permutedims(d0)
        end
        if size(idx,1) != 3
            idx = permutedims(idx)
        end

        push!(dicts, d0)
        push!(idxs,  idx)
    end

    dict0 = permutedims(vcat(dicts...))
    idx   = permutedims(hcat(idxs...))

    f = matopen(outfile, "w")
    write(f, "dict0", dict0)
    write(f, "idx",   idx)
    close(f)
end

main()
