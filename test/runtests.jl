using MovingBoundaryMinerals
using Test

push!(LOAD_PATH, "..")

function runtests()
    testdir = pwd()
    istest(f) = endswith(f, ".jl") && startswith(basename(f), "test_")
    testfiles = sort(
        filter(
            istest,
            vcat([joinpath.(root, files) for (root, dirs, files) in walkdir(testdir)]...),
        ),
    )
    nfail = 0
    printstyled("Testing package MovingBoundaryMinerals.jl\n"; bold=true, color=:white)

    for f in testfiles

        println("")
        println("Running tests from $f")
        try
            run(`$(Base.julia_cmd()) -O3 --startup-file=no $(joinpath(testdir, f))`)
        catch ex
            nfail += 1
        end
    end
    return nfail
end

exit(runtests())
