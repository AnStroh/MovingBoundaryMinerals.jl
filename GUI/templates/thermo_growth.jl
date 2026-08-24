function thermo_growth_page()
    fields = join([
        number_field("Ri", "Initial interface radius [m]", 0.0001; min = 0),
        select_field("n", "Geometry", [1 => "Planar", 2 => "Cylindrical", 3 => "Spherical"], 1),
        textarea_field("path", "Temperature-time path", "0, 1400\n30, 1350";
                       hint = "One 'time [days], temperature [°C]' point per line, at least 2. Time must start " *
                              "at 0 and strictly increase; temperature can go up and down between lines - add " *
                              "more rows for a non-monotonic path (e.g. '0, 1400' / '15, 1420' / '30, 1350')."),
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
    Temperature is linearly interpolated between the nodes of the temperature-time path below -
    the simplest path is a straight line from a start to an end temperature (falling drives
    growth, rising drives resorption), but adding more nodes gives an arbitrary, non-monotonic
    path where temperature can go up and down over the course of the run.</p>
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
