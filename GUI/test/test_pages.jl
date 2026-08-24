@testset "Page routes" begin
    @test get_page("/").status == 302   # redirects to /single-crystal

    # Match against <h1>title</h1> specifically, not just any occurrence of the title text -
    # the shared nav bar repeats every page's title on every page (see templates/layout.jl),
    # so a bare `occursin` check here would pass even if the wrong page's body were returned.
    for (path, title) in (
        "/single-crystal" => "Single-crystal diffusion",
        "/diffusion-couple" => "Diffusion couple",
        "/thermo-growth" => "Thermodynamically constrained growth",
        "/history" => "Past runs",
    )
        resp = get_page(path)
        @test resp.status == 200
        @test occursin("<h1>$(title)</h1>", String(resp.body))
    end
end
