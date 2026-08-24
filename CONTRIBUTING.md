# Contributing

[MovingBoundaryMinerals.jl](https://github.com/AnStroh/MovingBoundaryMinerals.jl) is an open-source project and we are very happy to accept contributions from the community. Please feel free to [open issues](https://github.com/AnStroh/MovingBoundaryMinerals.jl/issues/new) or submit patches (preferably as [pull requests](https://github.com/AnStroh/MovingBoundaryMinerals.jl/pulls)) any time. For planned larger contributions, it is often beneficial to get in contact with one of the principal developers first (see [AUTHORS.md](AUTHORS.md)).

[MovingBoundaryMinerals.jl](https://github.com/AnStroh/MovingBoundaryMinerals.jl) and its contributions are licensed under the MIT license. As a contributor, you certify that all your contributions are in conformance with the *Developer Certificate of Origin (Version 1.1)*, which is reproduced below.

## Development setup

To work on the package itself, clone the repository and run everything from its root directory.

**Running the test suite:**
```julia-repl
julia --project=. -e 'using Pkg; Pkg.test()'
```
or, from within a Julia REPL started with `julia --project=.`:
```julia-repl
julia> ]
(MovingBoundaryMinerals) pkg> test
```
The test suite can take a while to run.

**Running the GUI test suite:**

The `GUI/` folder is its own isolated Julia project (see [Architecture](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/gui/#gui-architecture)), so its tests run separately, with multiple threads (a run is dispatched on a background thread - see `GUI/jobs.jl`):
```julia-repl
julia --project=GUI -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=GUI GUI/test/runtests.jl
```
It spins up every route in-process via `Oxygen.internalrequest` (no real server, no browser needed) against a throwaway results directory, so it never touches a real `GUI/results/` folder you might have runs saved in.

**Building the documentation locally:**
```julia-repl
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```
This writes the built site to `docs/build/`; open `docs/build/index.html` in a browser to preview it. The first run will be slower while dependencies precompile.

## Developer Certificate of Origin (Version 1.1)
The following text was taken from [https://developercertificate.org](https://developercertificate.org):

    Developer Certificate of Origin
    Version 1.1

    Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
    1 Letterman Drive
    Suite D4700
    San Francisco, CA, 94129

    Everyone is permitted to copy and distribute verbatim copies of this
    license document, but changing it is not allowed.


    Developer's Certificate of Origin 1.1

    By making a contribution to this project, I certify that:

    (a) The contribution was created in whole or in part by me and I
        have the right to submit it under the open source license
        indicated in the file; or

    (b) The contribution is based upon previous work that, to the best
        of my knowledge, is covered under an appropriate open source
        license and I have the right under that license to submit that
        work with modifications, whether created in whole or in part
        by me, under the same open source license (unless I am
        permitted to submit under a different license), as indicated
        in the file; or

    (c) The contribution was provided directly to me by some other
        person who certified (a), (b) or (c) and I have not modified
        it.

    (d) I understand and agree that this project and the contribution
        are public and that a record of the contribution (including all
        personal information I submit with it, including my sign-off) is
        maintained indefinitely and may be redistributed consistent with
        this project or the open source license(s) involved.
