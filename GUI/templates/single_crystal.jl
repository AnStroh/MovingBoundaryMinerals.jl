function single_crystal_page()
    fields = join([
        number_field("D0", "Pre-exponential factor D0 [m²/s]", 2.75e-6; min = 0,
                     hint = "The diffusion coefficient's high-temperature limit, before the temperature correction below is applied."),
        number_field("Ea1", "Activation energy [J/mol]", 292880.0; min = 0,
                     hint = "How strongly diffusion speeds up with temperature - higher values mean diffusion is more sensitive to temperature."),
        number_field("L", "Domain length [m]", 0.001; min = 0),
        number_field("Cstart", "Initial composition [-]", 4.0; min = 0),
        number_field("Cinf", "Composition at infinity [-]", 0.0; min = 0),
        select_field("n", "Geometry", [1 => "Planar", 2 => "Cylindrical", 3 => "Spherical"], 1),
        number_field("Tstart_C", "Starting temperature [°C]", 1000.0),
        number_field("Tstop_C", "End temperature [°C]", 650.0),
        number_field("t_tot_years", "Total time [years]", 100.0; min = 0),
        number_field("nx", "Number of nodes", 100; step = "1", min = 2),
        number_field("CFL", "CFL number", 0.99; min = 0,
                     hint = "Controls the numerical time-step size relative to the grid spacing. Lower = smaller, safer steps but slower to compute."),
        text_field("run_name", "Name this run (optional)"; placeholder = "e.g. test-1"),
    ], "\n")

    body = """
    <p>Diffusion of a single component into a homogeneous, semi-infinite crystal - no moving boundary, no second phase.</p>
    <form id="run-form">
      $(fields)
      <button type="submit" id="run-button">Run simulation</button>
    </form>
    """
    return page("Single-crystal diffusion", "single-crystal", body)
end
