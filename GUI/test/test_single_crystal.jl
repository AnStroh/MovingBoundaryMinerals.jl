@testset "Single crystal: full run via HTTP" begin
    resp, data = post_form("/single-crystal/run", merge(VALID_SINGLE_CRYSTAL_PARAMS, Dict("run_name" => "gui-test")))
    @test resp.status == 200
    @test haskey(data, :id)

    result = poll_until_done(data.id)
    @test result.status == "done"
    @test occursin("Done.", result.summary)
    @test haskey(result, :xlsx_url)
    @test haskey(result, :jld2_url)

    folder = joinpath(TEST_RESULTS_DIR, result.folder)
    @test isdir(folder)
    check_result_files(folder, result)

    jl = JLD2.load(joinpath(folder, "profiles.jld2"))
    @test jl["source"] == CITATION_LINE
    @test length(jl["x_initial"]) == 100   # nx from VALID_SINGLE_CRYSTAL_PARAMS
    @test jl["parameters"]["L"] == 0.001

    XLSX.openxlsx(joinpath(folder, "profiles.xlsx")) do xf
        @test issetequal(XLSX.sheetnames(xf), ["Profile", "Parameters"])
        @test xf["Profile"]["A1"] == CITATION_LINE
        @test xf["Parameters"]["A1"] == CITATION_LINE
    end
end
