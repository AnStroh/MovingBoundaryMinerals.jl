"""
    page(title, mode_slug, body_html) -> String

Wraps `body_html` in the shared page shell: nav between the 3 modes, shared CSS, and the
run/poll/cancel JavaScript used by every mode's form. `mode_slug` is the URL path segment
(e.g. "single-crystal") this page's form should POST to.
"""
function page(title::String, mode_slug::String, body_html::String)
    return """
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="UTF-8">
      <title>$(title) - MovingBoundaryMinerals.jl</title>
      <link rel="icon" href="/static/favicon.ico">
      <link rel="stylesheet" href="/static/style.css">
      <script>
      (function() {
        const saved = localStorage.getItem("theme");
        const theme = saved || (window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
        document.documentElement.setAttribute("data-theme", theme);
      })();
      </script>
    </head>
    <body>
      <nav class="nav">
        <a href="/single-crystal">Single crystal</a>
        <a href="/diffusion-couple">Diffusion couple</a>
        <a href="/thermo-growth">Thermodynamic growth</a>
        <a href="/history">Past runs</a>
        <button type="button" id="theme-toggle" class="theme-toggle" aria-label="Toggle dark/light theme"></button>
      </nav>
      <main>
        <h1>$(title)</h1>
        $(body_html)
        <div id="status" class="status"></div>
        <div id="result"></div>
      </main>
      <script>
      const modeSlug = "$(mode_slug)";
      const form = document.getElementById("run-form");
      const runButton = document.getElementById("run-button");
      const statusDiv = document.getElementById("status");
      const resultDiv = document.getElementById("result");
      let pollTimer = null;
      let elapsedTimer = null;
      let runStart = null;
      let currentJobId = null;
      let currentProgress = 0;

      const themeToggle = document.getElementById("theme-toggle");
      function updateThemeIcon() {
        const isDark = document.documentElement.getAttribute("data-theme") === "dark";
        themeToggle.textContent = isDark ? "☀️ Light" : "🌙 Dark";
      }
      updateThemeIcon();
      themeToggle.addEventListener("click", () => {
        const isDark = document.documentElement.getAttribute("data-theme") === "dark";
        const next = isDark ? "light" : "dark";
        document.documentElement.setAttribute("data-theme", next);
        localStorage.setItem("theme", next);
        updateThemeIcon();
      });

      // History page only - deleting a saved run. Wired unconditionally since
      // querySelectorAll simply finds nothing on pages without a history table.
      document.querySelectorAll(".delete-run[data-folder]").forEach((btn) => {
        btn.addEventListener("click", async () => {
          if (!confirm("Delete this run permanently? This cannot be undone.")) return;
          const folder = btn.dataset.folder;
          btn.disabled = true;
          const resp = await fetch(`/results/\${folder}/delete`, { method: "POST" });
          if (resp.ok) {
            document.getElementById(`run-\${folder}`)?.remove();
          } else {
            alert("Could not delete this run.");
            btn.disabled = false;
          }
        });
      });

      // History page only - deleting every saved run at once. Only present when there's
      // at least one run (see history.jl), so the optional chaining below is just a guard
      // for every other page, which has no history table at all.
      document.getElementById("delete-all-runs")?.addEventListener("click", async (e) => {
        if (!confirm("Delete ALL saved runs permanently? This cannot be undone.")) return;
        const btn = e.currentTarget;
        btn.disabled = true;
        const resp = await fetch("/results/delete-all", { method: "POST" });
        if (resp.ok) {
          location.reload();
        } else {
          alert("Could not delete all runs.");
          btn.disabled = false;
        }
      });

      function formatElapsed(seconds) {
        const m = Math.floor(seconds / 60);
        const s = Math.floor(seconds % 60);
        return m > 0 ? `\${m}m \${s}s` : `\${s}s`;
      }

      function renderRunning() {
        const elapsed = (Date.now() - runStart) / 1000;
        const pct = Math.round(currentProgress * 100);
        statusDiv.innerHTML = `Running... \${formatElapsed(elapsed)} elapsed
          <progress value="\${currentProgress}" max="1"></progress>
          <span class="progress-pct">\${pct}%</span>
          <button type="button" id="cancel-button" class="cancel">Cancel</button>`;
        document.getElementById("cancel-button").addEventListener("click", cancelRun);
      }

      async function cancelRun() {
        if (!currentJobId) return;
        statusDiv.textContent = "Cancelling...";
        await fetch(`/jobs/\${currentJobId}/cancel`, { method: "POST" });
      }

      // The history page has no form - only attach form behaviour when one exists.
      if (form) {
        // Remember the last-used values per mode, so repeat runs don't start from scratch.
        // `run_name` is deliberately excluded - each run typically wants a fresh name.
        const formStorageKey = (fieldName) => `mbm_form_\${modeSlug}_\${fieldName}`;
        for (const el of form.elements) {
          if (!el.name || el.name === "run_name") continue;
          const saved = localStorage.getItem(formStorageKey(el.name));
          if (saved !== null) el.value = saved;
        }
        form.addEventListener("change", () => {
          for (const el of form.elements) {
            if (!el.name || el.name === "run_name") continue;
            localStorage.setItem(formStorageKey(el.name), el.value);
          }
        });

        form.addEventListener("submit", async (e) => {
          e.preventDefault();
          runButton.disabled = true;
          resultDiv.innerHTML = "";
          statusDiv.textContent = "Starting...";
          const body = new URLSearchParams(new FormData(form));
          const resp = await fetch(`/\${modeSlug}/run`, { method: "POST", body });
          if (resp.status === 409) {
            statusDiv.textContent = "A simulation is already running - please wait for it to finish.";
            runButton.disabled = false;
            return;
          }
          if (!resp.ok) {
            const errData = await resp.json().catch(() => ({}));
            statusDiv.textContent = errData.error || "Could not start the simulation - please check your inputs.";
            runButton.disabled = false;
            return;
          }
          const { id } = await resp.json();
          currentJobId = id;
          runStart = Date.now();
          currentProgress = 0;
          renderRunning();
          elapsedTimer = setInterval(renderRunning, 1000);
          pollTimer = setInterval(() => pollStatus(id), 1000);
        });
      }

      async function pollStatus(id) {
        const resp = await fetch(`/jobs/\${id}/status`);
        const data = await resp.json();
        if (data.status === "running") {
          if (typeof data.progress === "number") currentProgress = data.progress;
          return;
        }
        clearInterval(pollTimer);
        clearInterval(elapsedTimer);
        runButton.disabled = false;
        if (data.status === "done") {
          statusDiv.textContent = data.summary || "Done.";
          resultDiv.innerHTML = `
            <img src="\${data.png_url}" alt="Result plot">
            <div class="downloads">
              <a href="\${data.png_url}" download>Download PNG</a>
              <a href="\${data.pdf_url}" download>Download PDF</a>
              <a href="\${data.data_url}" download>Download data (final profile)</a>
              <a href="\${data.xlsx_url}" download>Download .xlsx</a>
              <a href="\${data.jld2_url}" download>Download .jld2</a>
              <a href="\${data.folder_url}">Open full results folder</a>
            </div>
            <p class="hint">Saved under <code>GUI/results/\${data.folder}/</code> - includes the plot (PNG + PDF, 300 dpi),
            the initial and final profile data (tab-delimited, plus bundled together with every input parameter
            in <code>.xlsx</code> and <code>.jld2</code>), and <code>inputs.toml</code> recording
            every input used, so this run can be found again and reproduced exactly.</p>
          `;
        } else if (data.status === "cancelled") {
          statusDiv.textContent = "Cancelled.";
        } else {
          statusDiv.textContent = "Simulation failed: " + (data.error || "unknown error");
        }
      }
      </script>
    </body>
    </html>
    """
end

"""
    number_field(name, label, value; step="any", hint="", min=nothing) -> String

A labelled number input, pre-filled with `value`. `hint`, if given, is shown as a short
explanatory line under the label - use it for parameters whose name alone assumes background
the user may not have (activation energy, CFL, K_D, ...). `min`, if given, is enforced by the
browser - use it for quantities that can't physically be zero or negative (lengths, densities,
node counts, ...), so an obviously-wrong value is caught before the form is even submitted.
"""
function number_field(name::String, label::String, value; step = "any", hint::String = "", min = nothing)
    hint_html = isempty(hint) ? "" : """<small class="field-hint">$(hint)</small>"""
    min_html = min === nothing ? "" : " min=\"$(min)\""
    return """
    <label for="$(name)">$(label)</label>
    $(hint_html)
    <input type="number" id="$(name)" name="$(name)" value="$(value)" step="$(step)"$(min_html) required>
    """
end

"""
    text_field(name, label; placeholder="") -> String

A labelled, optional free-text input.
"""
function text_field(name::String, label::String; placeholder::String = "")
    return """
    <label for="$(name)">$(label)</label>
    <input type="text" id="$(name)" name="$(name)" placeholder="$(placeholder)" maxlength="60">
    """
end

"""
    textarea_field(name, label, value; hint="", rows=6) -> String

A labelled, multi-line text input, pre-filled with `value`. Used for the thermo-growth
temperature-time path, where a single number field isn't enough to describe an arbitrary
number of (time, temperature) nodes.
"""
function textarea_field(name::String, label::String, value::String; hint::String = "", rows::Int = 6)
    hint_html = isempty(hint) ? "" : """<small class="field-hint">$(hint)</small>"""
    return """
    <label for="$(name)">$(label)</label>
    $(hint_html)
    <textarea id="$(name)" name="$(name)" rows="$(rows)" required>$(value)</textarea>
    """
end

"""
    select_field(name, label, options, selected) -> String

A labelled dropdown; `options` is a vector of (value, display_label) pairs.
"""
function select_field(name::String, label::String, options::Vector{<:Pair}, selected)
    opts = join(
        ["""<option value="$(v)"$(v == selected ? " selected" : "")>$(l)</option>""" for (v, l) in options],
        "\n"
    )
    return """
    <label for="$(name)">$(label)</label>
    <select id="$(name)" name="$(name)">
    $(opts)
    </select>
    """
end
