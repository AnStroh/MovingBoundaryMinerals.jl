function diffusion_couple_page()
    fields = join([
        number_field("D0_left", "Pre-exponential factor, left [m²/s]", 2.75e-6; min = 0,
                     hint = "The diffusion coefficient's high-temperature limit, before the temperature correction (activation energy) is applied."),
        number_field("D0_right", "Pre-exponential factor, right [m²/s]", 3.9e-7; min = 0),
        number_field("Ea_left", "Activation energy, left [J/mol]", 292879.6767; min = 0,
                     hint = "How strongly diffusion speeds up with temperature - higher values mean diffusion is more sensitive to temperature."),
        number_field("Ea_right", "Activation energy, right [J/mol]", 360660.4018; min = 0),
        number_field("Ri_interface", "Initial interface position [m]", 0.0002; min = 0),
        number_field("Ri_total", "Total domain length [m]", 0.0005; min = 0),
        number_field("Cl_i", "Initial composition, left [-]", 0.6; min = 0),
        number_field("Cr_i", "Initial composition, right [-]", 0.3; min = 0),
        number_field("KD_start", "Distribution coefficient K_D at t=0", 1.0; min = 0,
                     hint = "The equilibrium partitioning ratio of the diffusing component between the two phases at the interface."),
        number_field("KD_end", "Distribution coefficient K_D at end", 0.7; min = 0),
        number_field("V_ip", "Interface velocity [m/s]", 0.0,
                     hint = "0 keeps the interface stationary (Diff_couple_Flux.jl). Positive grows the left phase at the right phase's expense; negative resorbs it (Diff_couple_Flux_growth.jl)."),
        select_field("n", "Geometry", [1 => "Planar", 2 => "Cylindrical", 3 => "Spherical"], 3),
        number_field("Tstart_C", "Starting temperature [°C]", 1000.0),
        number_field("Tstop_C", "End temperature [°C]", 700.0),
        number_field("t_tot_years", "Total time [years]", 1000.0; min = 0),
        number_field("nx_left", "Number of nodes, left", 100; step = "1", min = 2),
        number_field("nx_right", "Number of nodes, right", 150; step = "1", min = 2),
        number_field("CFL", "CFL number", 0.5; min = 0,
                     hint = "Controls the numerical time-step size relative to the grid spacing. Lower = smaller, safer steps but slower to compute."),
        text_field("run_name", "Name this run (optional)"; placeholder = "e.g. test-1"),
    ], "\n")

    body = """
    <p>A diffusion couple with a moving interface, using a flux-balance condition to link the two phases. Leave interface velocity at 0 for a stationary interface with pure ion exchange, or set it non-zero for simultaneous growth/resorption.</p>
    <form id="run-form">
      $(fields)
      <button type="submit" id="run-button">Run simulation</button>
    </form>
    """
    return page("Diffusion couple", "diffusion-couple", body)
end
