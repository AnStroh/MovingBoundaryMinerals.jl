function thermo_growth_page()
    fields = join([
        number_field("Ri", "Initial interface radius [m]", 0.0005; min = 0),
        select_field("n", "Geometry", [1 => "Planar", 2 => "Cylindrical", 3 => "Spherical"], 1),
        textarea_field("path", "Temperature-time path", "0, 1400\n30, 1350";
                       hint = "One 'time [days], temperature [°C]' point per line, at least 2. Time must start " *
                              "at 0 and strictly increase; temperature can go up and down between lines - add " *
                              "more rows for a non-monotonic path (e.g. '0, 1400' / '15, 1420' / '30, 1350'). " *
                              "Only the relative spacing between these times matters - the field below sets how " *
                              "long the run actually lasts, stretching or compressing this path to fit."),
        number_field("t_tot_days", "Total time [days]", 30.0; min = 0,
                     hint = "The path above is rescaled to span exactly this many days, so you can speed up or slow down a run without rewriting every row."),
        number_field("P", "Pressure [Pa]", 1.0e6; min = 0),
        number_field("D0", "Pre-exponential factor D0, crystal side [m²/s]", 5.38e-9; min = 0,
                     hint = "The crystal (left) side's diffusion coefficient high-temperature limit, before the temperature/pressure/composition corrections are applied."),
        number_field("Ea", "Activation energy, crystal side [J/mol]", 226000.0; min = 0,
                     hint = "How strongly diffusion speeds up with temperature on the crystal side - higher values mean diffusion is more sensitive to temperature."),
        number_field("deltaV", "Activation volume [m³/mol]", 7e-6; min = 0,
                     hint = "How strongly diffusion responds to pressure, via the pressure-correction term in the Arrhenius equation."),
        number_field("D0_right", "Pre-exponential factor D0, melt/fluid side [m²/s]", 0.0003634023264950478; min = 0,
                     hint = "The melt/fluid (right) side's diffusion coefficient high-temperature limit."),
        number_field("Ea_right", "Activation energy, melt/fluid side [J/mol]", 218022.084784; min = 0,
                     hint = "How strongly diffusion speeds up with temperature on the melt/fluid side."),
        number_field("alpha", "Angle to [001] axis [°]", 0.0,
                     hint = "Angle between the crystal's [001] axis and the profile direction. Together with beta/gamma below, this sets how anisotropic diffusion projects onto the modelled direction."),
        number_field("beta", "Angle to [010] axis [°]", 90.0,
                     hint = "Angle between the crystal's [010] axis and the profile direction."),
        number_field("gamma", "Angle to [100] axis [°]", 90.0,
                     hint = "Angle between the crystal's [100] axis and the profile direction."),
        number_field("CompInt", "Bulk composition of interest [-]", 0.50; min = 0,
                     hint = "The overall (bulk) composition of the crystal + melt/fluid system, which sets the starting composition via the phase diagram."),
        number_field("nx_left", "Number of nodes, left", 50; step = "1", min = 2),
        number_field("nx_right", "Number of nodes, right", 50; step = "1", min = 2),
        number_field("CFL", "CFL number", 50.0; min = 0,
                     hint = "Controls the numerical time-step size. This mode's total time is set by the path above rather than a fixed CFL step count, so a higher CFL here (larger, faster steps) is the main way to speed up a slow run."),
        text_field("run_name", "Name this run (optional)"; placeholder = "e.g. test-1"),
    ], "\n")

    body = """
    <p>Thermodynamically constrained crystal growth/resorption (olivine, following Dohmen &amp; Chakraborty 2007):
    interface compositions come from a digitized phase diagram rather than a fixed K<sub>D</sub>.
    Diffusion is anisotropic: the three angles below set how the crystal's [001]/[010]/[100] axes
    are oriented relative to the modelled profile direction, which changes the effective diffusion
    coefficient (the defaults - 0°/90°/90° - align the profile with [001]).
    Temperature is linearly interpolated between the nodes of the temperature-time path below -
    the simplest path is a straight line from a start to an end temperature (falling drives
    growth, rising drives resorption), but adding more nodes gives an arbitrary, non-monotonic
    path where temperature can go up and down over the course of the run. Geometry defaults to
    Planar (n=1), matching the examples this mode is based on.</p>
    <p><strong>This mode is much slower than the other two</strong> - with the default settings shown here,
    a run can take on the order of tens of minutes. Shortening the path's total time, raising the
    CFL number, or lowering the node counts below (fewer, larger steps/cells) will speed it up, at
    some cost to accuracy.</p>
    <form id="run-form">
      $(fields)
      <button type="submit" id="run-button">Run simulation</button>
    </form>
    """
    return page("Thermodynamically constrained growth", "thermo-growth", body)
end
