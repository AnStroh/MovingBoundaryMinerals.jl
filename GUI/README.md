# MovingBoundaryMinerals.jl GUI

A local, browser-based GUI for running three of the package's core models without writing any Julia code: single-crystal diffusion, a moving-interface diffusion couple, and thermodynamically constrained crystal growth/resorption.

Full documentation (with a step-by-step walkthrough): [GUI docs page](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/man/gui/).

## Starting it

Double-click one of these, from this `GUI/` folder:

- `run_gui.sh` (Linux)
- `run_gui.command` (macOS)
- `run_gui.bat` (Windows)

Or, from a terminal at the **repository root** (not this folder):

```bash
julia -t auto --project=GUI GUI/app.jl
```

This opens automatically in Google Chrome (Chromium/Edge also work) as a standalone app window. If none of those browsers are installed, or nothing opens automatically, the server URL (`http://127.0.0.1:8811`) is printed to the terminal — open it in any browser, including Firefox or Safari.

## Using it

1. Pick a mode from the navigation bar: Single crystal, Diffusion couple, or Thermodynamic growth.
2. The form comes pre-filled with working defaults — click **Run simulation** right away, or change any values first.
3. Once it finishes, download the result plot (PNG/PDF, 300 dpi) and the raw profile data (tab-delimited, or bundled with every input parameter into a single `.xlsx` or `.jld2` file), or find any past run again later via **Past runs** in the nav bar.

Every completed run is saved permanently under `GUI/results/<timestamp>[_<name>]/`, including an `inputs.toml` recording every parameter used, so it can be reproduced exactly. Delete a single run or every saved run from the **Past runs** page at any time.

## This is a simplified entry point, not the full package

The forms only expose the parameters most users would want to change; several settings (grid refinement strategy, crystallographic angles, non-monotonic temperature-time paths, alternative interface conditions) are fixed at sensible defaults and not shown here. For full control, edit the [example scripts](../examples/) directly, or adapt the backend functions in [`backends/`](backends/) — they're ordinary functions built on the package's public API, so nothing here is special or closed off.

## Contents

| Path | What it is |
|:--|:--|
| `app.jl` | Entry point: HTTP routes, job scheduling, browser launch |
| `backends/` | The three simulation functions the forms call (`run_single_crystal`, `run_diffusion_couple`, `run_thermo_growth`) |
| `templates/` | HTML page generation |
| `static/` | CSS |
| `jobs.jl`, `plotting.jl`, `results.jl` | Background-job handling, plotting, and saving results to disk |
| `results/` | Where completed runs are saved (created automatically, not tracked in git) |

See [CONTRIBUTING.md](../CONTRIBUTING.md) at the repository root for general development setup (running tests, building the docs).
