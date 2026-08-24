# Changelog

All notable changes to `MovingBoundaryMinerals.jl` are documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- Example D2: diffusion-limited crystal growth of olivine driven by a user-defined, non-monotonic temperature-time path (temperature may increase and decrease between user-specified nodes, unlike D1's monotonic cooling path).
- A local, browser-based GUI (`GUI/`) for running single-crystal diffusion, moving-interface diffusion couple, and thermodynamically constrained growth/resorption models without writing Julia code: pre-filled forms with inline explanations of the less obvious parameters, cancellable background runs with an elapsed-time display and a live progress bar (`t / t_tot`), and downloadable results (plot at 300 dpi as PNG/PDF, raw profile data, and an `inputs.toml` recording every parameter for exact reproducibility) organized by a "Past runs" history page. See `GUI/README.md` and the [GUI documentation page](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/gui/).
- A package logo and favicon (`docs/src/assets/`), shown on the documentation site and in the README.
- Documentation: a [Getting Started](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/getting_started/) guide, an [Interpreting Output](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/interpreting_output/) page, and a [Quick Reference](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/quick_reference/) page; restored and expanded the General Remarks, Stefan Problem, and Digitization pages; added a proper bibliography via DocumenterCitations.jl; split the function reference into public API vs. internal helpers.
- `Stable` documentation badge, alongside the existing `Dev` badge.
- Documentation: [Boundary Conditions](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/boundary_conditions/) (the equations behind the outer Dirichlet/Neumann conditions and the three inner interface conditions — flux-balance, total mass-balance, and the thermodynamically constrained Stefan condition), [Mesh Refinement](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/mesh_refinement/) (the m-refinement and h-refinement grading equations behind each `RefineMethod`), and [Configuration Options](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/configuration_options/) (the switches that change which code path is used) pages; expanded Numerical Approach with the governing diffusion equation, the Arrhenius law, and the FEM/implicit-time discretization; expanded the Stefan Problem page with the interface-velocity equation and the `CompInt` lever-rule/domain-length relation.
- GUI: a **Delete all runs** button on the "Past runs" page, alongside the existing per-run **Delete** button, for clearing out `GUI/results/` without leaving the browser; both ask for confirmation first.
- GUI: every run now also saves `profiles.xlsx` and `profiles.jld2` (via XLSX.jl and JLD2.jl respectively) alongside the existing `.tab`/`inputs.toml` outputs — each bundles the initial/final profile data together with every input parameter and result into a single file, for users working in Excel or Julia. `inputs.toml` and both bundled exports now also carry the same small citation line already stamped on the plot and `.tab` files. (A `.mat` export via MAT.jl was tried and dropped: its only unconditional dependency is HDF5.jl, which pulls in four MPI implementations and AWS S3 client libraries for no reason a local single-array write needs — too heavy and failure-prone to install on the locked-down networks common in (geo)science/industry settings. `.xlsx` opens fine in MATLAB too.)
- GUI: an automated test suite (`GUI/test/`, run via `julia -t auto --project=GUI GUI/test/runtests.jl`), exercising every route in-process (no real server or browser) against a throwaway results directory — covers all three backends end to end, input validation, job cancellation, the "already running" conflict, and both delete endpoints. Wired into CI (`.github/workflows/GUI.yml`) on every change under `GUI/` or `src/`. `RESULTS_DIR` is now overridable via the `MBM_GUI_RESULTS_DIR` environment variable to make this possible; `app.jl`'s `main()` only runs when the file is launched directly, not when it's `include`d for testing.
- Documentation: [Benchmarks](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/benchmarks/) (every example tested against a closed-form analytical or semi-analytical solution — A1, A2, B1, B2, B4, B5, C1) and [Example Gallery](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/example_gallery/) (B6, B7, C2, D1 — realistic/illustrative output with no closed-form solution to compare against, kept separate from the benchmarks since they validate nothing) pages, each figure reproduced from [Stroh2025] with attribution.

### Fixed
- `CI.yml` no longer runs (and spuriously fails) on the `gh-pages` branch created by documentation deployment.
- A docstring-formatting bug (a stray blank line between the docstring and the `function` keyword) that silently dropped the docstrings of 13 exported functions, including `calculate_dt`, `save_figure`, and `h_refinement1`/`h_refinement2`.
- `MovingBoundaryMinerals.Benchmarks` functions (e.g. `analytical_sol_step_function`, `smith`, `rayleigh_fractionation`) are now properly reachable via `using MovingBoundaryMinerals` alone, instead of requiring a separate `using MovingBoundaryMinerals.Benchmarks`.
- Several thin/incomplete docstrings (missing arguments, missing units, incorrect type annotations) across `functions.jl`.
- Numerous spelling and grammar issues across the documentation and source comments.
- Documentation previously conflated the flux-balance interface condition with the (thermodynamically constrained) Stefan condition in several places (e.g. describing the interface velocity `V_ip` as solved-for in the flux-balance/mass-balance cases, where it is actually a prescribed input; and the `Residual` diagnostic as applying to all moving-interface examples rather than just the Stefan-condition family).

### Removed
- Dead package exports (`digitise_plot`, `ndgrid`, `lasaga`, `preallocations`) that referenced functions not defined anywhere in the package.

## [1.0.0] - 2025-12-02

### Added
- h-refinement grid strategies (`h_refinement1`, `h_refinement2`) as alternatives to m-refinement, selectable via `RefineMethod`.
- DOI badge and citation information in the README.

### Changed
- `nPoints` for h-refinement can now be given as a single number (both sides equal) or as `[left right]`.
- General code cleanup and example fixes.

## [0.1.0] - 2025-05-29

Initial public release: diffusion-couple and moving-boundary (Stefan) solvers, with benchmarked examples (A, B, C series) and the thermodynamically constrained crystal-growth example D1.

[Unreleased]: https://github.com/AnStroh/MovingBoundaryMinerals.jl/compare/v1.0.0...main
[1.0.0]: https://github.com/AnStroh/MovingBoundaryMinerals.jl/compare/v0.1.0...v1.0.0
[0.1.0]: https://github.com/AnStroh/MovingBoundaryMinerals.jl/releases/tag/v0.1.0
