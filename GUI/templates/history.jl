const MODE_DISPLAY_NAME = Dict(
    "single_crystal" => "Single crystal",
    "diffusion_couple" => "Diffusion couple",
    "thermo_growth" => "Thermodynamic growth",
)

"""
    history_page(runs) -> String

`runs` is a vector of `(; folder, mode, timestamp)` named tuples, already sorted
newest-first. Renders a table linking to each run's saved plot and results folder.
"""
function history_page(runs::Vector)
    rows = if isempty(runs)
        """<tr><td colspan="5">No runs yet - results will appear here after you run a simulation.</td></tr>"""
    else
        join([
            """
            <tr id="run-$(r.folder)">
              <td>$(get(MODE_DISPLAY_NAME, r.mode, r.mode))</td>
              <td>$(r.timestamp)</td>
              <td>$(isempty(r.run_name) ? "" : r.run_name)</td>
              <td><a href="/results/$(r.folder)/plot.png"><img class="thumb" src="/results/$(r.folder)/plot.png" alt="plot"></a></td>
              <td>
                <a href="/results/$(r.folder)/plot.png" download>PNG</a> ·
                <a href="/results/$(r.folder)/plot.pdf" download>PDF</a> ·
                <a href="/results/$(r.folder)/">Folder</a> ·
                <button type="button" class="delete-run" data-folder="$(r.folder)">Delete</button>
              </td>
            </tr>
            """ for r in runs
        ], "\n")
    end

    delete_all_html = isempty(runs) ? "" : """
    <button type="button" id="delete-all-runs" class="delete-run delete-all">Delete all runs</button>
    """

    body = """
    <p>Every simulation you run is saved to its own timestamped folder under <code>GUI/results/</code>,
    together with the exact inputs used (<code>inputs.toml</code>) and the raw profile data - nothing is
    ever overwritten. This page lists them all, newest first.</p>
    $(delete_all_html)
    <table class="history">
      <thead><tr><th>Mode</th><th>Run</th><th>Name</th><th>Plot</th><th>Downloads</th></tr></thead>
      <tbody>
      $(rows)
      </tbody>
    </table>
    """
    return page("Past runs", "history", body)
end
