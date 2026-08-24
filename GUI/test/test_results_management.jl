@testset "Results management" begin
    # Two completed runs to work with, independent of whatever other test files left behind.
    # Distinct run_names matter here: save_run_outputs names folders by second-granularity
    # timestamp, so two runs started within the same wall-clock second (easily the case, since
    # these are fast) would otherwise collide on an unnamed folder and silently overwrite one
    # another - a real sharp edge in the app itself, not just a test artifact.
    _, data1 = post_form("/single-crystal/run", merge(VALID_SINGLE_CRYSTAL_PARAMS, Dict("run_name" => "delete-test-1")))
    result1 = poll_until_done(data1.id)
    folder1 = joinpath(TEST_RESULTS_DIR, result1.folder)
    @test isdir(folder1)

    _, data2 = post_form("/single-crystal/run", merge(VALID_SINGLE_CRYSTAL_PARAMS, Dict("run_name" => "delete-test-2")))
    result2 = poll_until_done(data2.id)
    folder2 = joinpath(TEST_RESULTS_DIR, result2.folder)
    @test isdir(folder2)
    @test folder1 != folder2

    @testset "/history lists both runs" begin
        resp = get_page("/history")
        body = String(resp.body)
        @test occursin(result1.folder, body)
        @test occursin(result2.folder, body)
        @test occursin("id=\"delete-all-runs\"", body)   # the button itself, only rendered when at least one run exists
    end

    @testset "delete a single run" begin
        resp, data = post_form("/results/$(result1.folder)/delete", Dict())
        @test resp.status == 200
        @test data.ok == true
        @test !isdir(folder1)
        @test isdir(folder2)   # untouched

        # Deleting the same folder again -> 404, not a silent success.
        resp2, data2b = post_form("/results/$(result1.folder)/delete", Dict())
        @test resp2.status == 404
        @test data2b.ok == false
    end

    @testset "delete all runs" begin
        resp, data = post_form("/results/delete-all", Dict())
        @test resp.status == 200
        @test data.ok == true
        @test !isdir(folder2)
        @test isempty(readdir(TEST_RESULTS_DIR))

        resp2 = get_page("/history")
        @test occursin("No runs yet", String(resp2.body))
        # The JS reference to #delete-all-runs is in the shared page layout unconditionally
        # (see templates/layout.jl) - only the button element itself is conditional on there
        # being runs to delete (see history.jl), so check for that specifically.
        @test !occursin("id=\"delete-all-runs\"", String(resp2.body))
    end

    @testset "deleting all when nothing is saved doesn't error" begin
        resp, data = post_form("/results/delete-all", Dict())
        @test resp.status == 200
        @test data.ok == true
    end
end
