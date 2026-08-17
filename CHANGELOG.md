# Changelog

All notable changes to `MovingBoundaryMinerals.jl` are documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- Example D2: diffusion-limited crystal growth of olivine driven by a user-defined, non-monotonic temperature-time path (temperature may increase and decrease between user-specified nodes, unlike D1's monotonic cooling path).
- A local, browser-based GUI (`GUI/`) for running single-crystal diffusion, moving-interface diffusion couple, and thermodynamically constrained growth/resorption models without writing Julia code: pre-filled forms with inline explanations of the less obvious parameters, cancellable background runs with an elapsed-time display, and downloadable results (plot at 300 dpi as PNG/PDF, raw profile data, and an `inputs.toml` recording every parameter for exact reproducibility) organized by a "Past runs" history page. See `GUI/README.md` and the [GUI documentation page](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/gui/).
- A package logo and favicon (`docs/src/assets/`), shown on the documentation site and in the README.
- Documentation: a [Getting Started](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/getting_started/) guide, an [Interpreting Output](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/interpreting_output/) page, and a [Quick Reference](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/quick_reference/) page; restored and expanded the General Remarks, Stefan Problem, and Digitization pages; added a proper bibliography via DocumenterCitations.jl; split the function reference into public API vs. internal helpers.
- `Stable` documentation badge, alongside the existing `Dev` badge.

### Fixed
- `CI.yml` no longer runs (and spuriously fails) on the `gh-pages` branch created by documentation deployment.
- A docstring-formatting bug (a stray blank line between the docstring and the `function` keyword) that silently dropped the docstrings of 13 exported functions, including `calculate_dt`, `save_figure`, and `h_refinement1`/`h_refinement2`.
- `MovingBoundaryMinerals.Benchmarks` functions (e.g. `analytical_sol_step_function`, `smith`, `rayleigh_fractionation`) are now properly reachable via `using MovingBoundaryMinerals` alone, instead of requiring a separate `using MovingBoundaryMinerals.Benchmarks`.
- Several thin/incomplete docstrings (missing arguments, missing units, incorrect type annotations) across `functions.jl`.
- Numerous spelling and grammar issues across the documentation and source comments.

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
