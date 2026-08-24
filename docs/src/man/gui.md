```@meta
CurrentModule = MovingBoundaryMinerals
```
# [Graphical user interface (GUI)](@id gui)

For users who don't want to write or edit Julia code, the package ships with a local, browser-based GUI covering three of its core model families. It lives in the [`GUI/`](https://github.com/AnStroh/MovingBoundaryMinerals.jl/tree/main/GUI) folder at the repository root.

## Launching it

Double-click `GUI/run_gui.sh` (Linux), `GUI/run_gui.command` (macOS), or `GUI/run_gui.bat` (Windows) — or, from a terminal at the repository root:

```julia-repl
julia -t auto --project=GUI GUI/app.jl
```

`-t auto` is required, not optional — it lets long simulations run on a background thread without freezing the interface (see [Architecture](@ref gui-architecture) below). This starts a local server at `http://127.0.0.1:8811` and opens it automatically in its own window.

### Which browser does it use?

**Google Chrome is directly supported** and is the browser the launcher looks for first — if Chrome is installed, the GUI opens in a dedicated Chrome "app window" (no address bar or tabs, so it looks and feels like a standalone program rather than a website). Chromium and Microsoft Edge are supported the same way, in case Chrome isn't installed.

If none of Chrome, Chromium, or Edge are found, the GUI falls back to opening the URL in your system's default browser instead — you'll just see the normal browser chrome (address bar, tabs) around it rather than the dedicated app window. And regardless of which browser gets auto-opened, the GUI is a completely ordinary local web page: **any modern browser can be pointed at it manually**, including Firefox or Safari. If a window doesn't open automatically, or you'd simply rather use a different browser, the URL is always printed to the terminal (`http://127.0.0.1:8811`) — copy it into whichever browser you prefer.

## Step-by-step: your first run

1. Launch the GUI as described above and wait for the window/tab to open.
2. Use the navigation bar at the top to pick a mode: **Single crystal**, **Diffusion couple**, or **Thermodynamic growth** (see [What it covers](@ref) below for what each one models).
3. The form is already filled in with working default values — you can click **Run simulation** immediately without changing anything, or adjust any of the numbers first. Fields with less obvious names (activation energy, CFL number, K_D, ...) have a short explanation printed just under their label.
4. Once running, the page shows an elapsed-time counter and a **Cancel** button. Only one simulation runs at a time; starting a second one while another is running isn't possible until the first finishes or is cancelled.
5. When it finishes, the resulting composition-profile plot appears on the page, along with **Download PNG**, **Download PDF**, and **Download data** links, and a link to open the full results folder.
6. To find a run again later, use the **Past runs** link in the navigation bar — every run you've ever completed is listed there, newest first, with its plot and downloads, whether or not the GUI has been restarted since.

## What it covers

Three modes, each a simplified, form-driven version of one example family:

| Mode | Based on | Physics |
|:--|:--|:--|
| Single crystal | [`Simple_Diff.jl`](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/Simple_Diff.jl) | Diffusion within a single homogeneous crystal, no moving boundary. |
| Diffusion couple | [`Diff_couple_Flux.jl`](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/Diff_couple_Flux.jl) | Two phases in contact with a moving interface, flux-balance condition. |
| Thermodynamic growth | [`D2.jl`](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/D2.jl) — a 2-point path reproduces [`D1.jl`](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/D1.jl)'s linear cooling, more points reproduce D2.jl's non-monotonic path | Thermodynamically constrained crystal growth/resorption from a digitized phase diagram. |

Fill in the form (every field is pre-filled with the same defaults as the underlying example, so you can just click "Run simulation" immediately) and the resulting composition profile is displayed once the run finishes. The more jargon-heavy fields (activation energy, CFL number, K_D, ...) have a short explanation shown right under their label.

!!! warning "The thermodynamic growth mode can take a long time"
    Depending on the chosen total time and resolution, a run can take from several minutes up to an hour or more with the default settings. The page shows a running elapsed-time counter and a **Cancel** button while a simulation is in progress; only one simulation runs at a time.

## Saved results and past runs

Every run that completes successfully is saved to its own timestamped folder under `GUI/results/` the moment it finishes, so results are never silently overwritten (cancelled or failed runs save nothing, since there's no valid result to keep):

- `plot.png` and `plot.pdf` — the result figure at 300 dpi
- `profile_initial.tab` and `profile_final.tab` — the raw (distance, composition) data, tab-delimited, for opening directly in Excel or any plotting/analysis tool
- `inputs.toml` — every input parameter used for the run, plus key scalar results (final time, mass-balance error), so the run can be found again and reproduced exactly
- `profiles.xlsx` and `profiles.jld2` — the initial/final profile data *and* everything from `inputs.toml`, bundled into a single file each for Excel and Julia (`profiles.jld2`, opened with `load(...)` from JLD2.jl) respectively — so you don't need to open `inputs.toml` separately just to see what a `.xlsx`/`.jld2` file was run with

Every one of these files — plots, profiles, and the bundled exports alike — carries a small citation line naming the package version and its Zenodo software archive ([StrohZenodo2025](@cite)), so a result stays traceable back to its source even if shared or moved on its own, without being visually intrusive (a thin grey footer on plots, a `#`-comment line in the `.tab` files, a `source` field/cell in `inputs.toml`/`.xlsx`/`.jld2`).

You can optionally name a run (e.g. "test-1") before clicking "Run simulation"; the name is appended to the folder's timestamp to make it easier to find later. Download links for all of these appear next to the result plot once a run finishes, and the **Past runs** page (linked in the navigation bar) lists every run you've ever done, newest first, with a thumbnail and download links for each — so results stay reachable even after closing and reopening the GUI, without needing to browse the file system directly.

Each run's row has its own **Delete** button, and a **Delete all runs** button appears above the table whenever at least one run is saved (it stays hidden on an empty history). Both ask for confirmation first and permanently remove the corresponding folder(s) under `GUI/results/` — there's no undo, so make sure you've downloaded anything you want to keep first.

## The GUI is a simplified entry point, not the full package

**The GUI intentionally exposes only a subset of what the package can do.** To keep the forms usable for someone with no programming background, each mode fixes several parameters at sensible defaults rather than exposing everything: the grid refinement strategy, crystallographic angles (mode 3), and the alternative interface conditions (total mass-balance) are not available through the form. (Thermodynamic growth is the exception — its temperature-time path field accepts any number of points, so both `D1.jl`'s linear path and `D2.jl`'s non-monotonic one are fully supported.) The GUI also only ever shows the final composition profile, not the fuller diagnostic plots (phase-diagram overlay, K_D evolution, mass-balance history, ...) the example scripts can produce.

**Directly working with the code always offers strictly more opportunities than the GUI does.** If you need a parameter, a boundary condition, a diagnostic, or a combination of settings the forms don't expose, you have two options, in increasing order of control:

1. Edit or copy one of the [example scripts](@ref "List of examples") directly (see [Getting Started](@ref getting-started)) — every parameter documented in the [Quick Reference](@ref quick-reference) is available there.
2. Adapt the GUI's own backend functions in `GUI/backends/*.jl` — these are ordinary Julia functions (`run_single_crystal`, `run_diffusion_couple`, `run_thermo_growth`) built directly on the package's public API (the same functions listed in [List of all functions](@ref)), so adding a field is a matter of exposing another keyword argument, not working around a restriction. The GUI is a deliberately narrowed view built on top of the same functions everyone else scripts against — nothing about it is special or closed off.

## [Architecture](@id gui-architecture)

For anyone maintaining or extending the GUI itself: it's a small [Oxygen.jl](https://github.com/OxygenFramework/Oxygen.jl) web app with its own isolated `GUI/Project.toml` (so its dependencies are never pulled in by a normal `] add MovingBoundaryMinerals`, mirroring how `docs/` has its own environment). Each "Run simulation" click starts the corresponding backend function on a background thread (`Threads.@spawn`, which is why `-t auto` matters) and the page polls for completion (`report_progress`, threaded through the backend the same way `should_stop` is, drives the progress bar); the result is rendered with `Plots.jl` (running headless, no display needed) straight to an in-memory PNG. See `GUI/app.jl` for the routes and `CONTRIBUTING.md` for general development setup.

`GUI/test/runtests.jl` exercises every route in-process via `Oxygen.internalrequest` — no real server or browser involved — against a throwaway results directory (`RESULTS_DIR` is overridable via the `MBM_GUI_RESULTS_DIR` environment variable for exactly this reason), covering the three backends end to end (including the `.xlsx`/`.jld2` exports and the citation stamp), input validation, job cancellation, the "already running" conflict, and both delete endpoints. It runs in CI on every change under `GUI/` or `src/` (see `.github/workflows/GUI.yml`); run it locally with `julia -t auto --project=GUI GUI/test/runtests.jl` (see `CONTRIBUTING.md`).
