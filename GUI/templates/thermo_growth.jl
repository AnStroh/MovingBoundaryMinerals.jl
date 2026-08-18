function thermo_growth_page()
    fields = join([
        number_field("Ri", "Initial interface radius [m]", 0.0001; min = 0),
        select_field("n", "Geometry", [1 => "Planar", 2 => "Cylindrical", 3 => "Spherical"], 1),
        number_field("Tstart_C", "Starting temperature [°C]", 1400.0),
        number_field("Tstop_C", "End temperature [°C]", 1350.0),
        number_field("t_tot_days", "Total time [days]", 30.0; min = 0),
        number_field("P", "Pressure [Pa]", 1.0e6; min = 0),
        number_field("D0", "Pre-exponential factor D0 [m²/s]", 5.38e-9; min = 0,
                     hint = "The diffusion coefficient's high-temperature limit, before the temperature correction (activation energy) is applied."),
        number_field("Ea", "Activation energy [J/mol]", 226000.0; min = 0,
                     hint = "How strongly diffusion speeds up with temperature - higher values mean diffusion is more sensitive to temperature."),
        number_field("CompInt", "Bulk composition of interest [-]", 0.50; min = 0,
                     hint = "The overall (bulk) composition of the crystal + melt/fluid system, which sets the starting composition via the phase diagram."),
        text_field("run_name", "Name this run (optional)"; placeholder = "e.g. test-1"),
    ], "\n")

    body = """
    <p>Thermodynamically constrained crystal growth/resorption (olivine, following Dohmen &amp; Chakraborty 2007):
    interface compositions come from a digitized phase diagram rather than a fixed K<sub>D</sub>.
    Temperature varies linearly between the start and end value over the given time - a falling
    temperature (start &gt; end) drives growth, a rising one (start &lt; end) drives resorption.</p>
    <p><strong>This mode is much slower than the other two</strong> - with the default settings shown here,
    a run can take on the order of tens of minutes. Reducing the total time and/or increasing the CFL number
    (fixed at their example defaults in this form) would speed it up, but requires editing the underlying
    script directly - see the documentation.</p>
    <form id="run-form">
      $(fields)
      <button type="submit" id="run-button">Run simulation</button>
    </form>
    """
    return page("Thermodynamically constrained growth", "thermo-growth", body)
end
