using Documenter
using MOBILE

# Get MOBILE.jl root directory
DC_root_dir = dirname(@__DIR__)

license = read(joinpath(DC_root_dir, "LICENSE.md"), String)
write(joinpath(@__DIR__, "src", "man", "license.md"), license)

security = read(joinpath(DC_root_dir, "SECURITY.md"), String)
write(joinpath(@__DIR__, "src", "man", "security.md"), security)

# Copy list of authors to not need to synchronize it manually
authors_text = read(joinpath(DC_root_dir, "AUTHORS.md"), String)
# authors_text = replace(authors_text, "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)")
write(joinpath(@__DIR__, "src", "man", "authors.md"), authors_text)

# Copy some files from the repository root directory to the docs and modify them
# as necessary
# Based on: https://github.com/ranocha/SummationByPartsOperators.jl/blob/0206a74140d5c6eb9921ca5021cb7bf2da1a306d/docs/make.jl#L27-L41
open(joinpath(@__DIR__, "src", "man", "license.md"), "w") do io
  # Point to source license file
  println(io, """
  ```@meta
  EditURL = "https://github.com/AnStroh/MOBILE.jl/blob/main/LICENSE.md"
  ```
  """)
  # Write the modified contents
  println(io, "# [License](@id license)")
  println(io, "")
  for line in eachline(joinpath(dirname(@__DIR__), "LICENSE.md"))
    line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
    println(io, "> ", line)
  end
end

open(joinpath(@__DIR__, "src", "man", "code_of_conduct.md"), "w") do io
  # Point to source license file
  println(io, """
  ```@meta
  EditURL = "https://github.com/AnStroh/MOBILE.jl/blob/main/CODE_OF_CONDUCT.md"
  ```
  """)
  # Write the modified contents
  println(io, "# [Code of Conduct](@id code-of-conduct)")
  println(io, "")
  for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
    line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
    println(io, "> ", line)
  end
end

open(joinpath(@__DIR__, "src", "man", "contributing.md"), "w") do io
    # Point to source license file
    println(io, """
    ```@meta
    EditURL = "https://github.com/AnStroh/MOBILE.jl/blob/main/CONTRIBUTING.md"
    ```
    """)
    # Write the modified contents
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
      line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
      line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
      println(io,line)
    end
  end
@info "Making documentation..."
makedocs(;
    sitename="MOBILE.jl",
    authors="Annalena Stroh, Evangelos Moulas and contributors",
    modules=[MOBILE],
    format=Documenter.HTML(; assets = ["assets/favicon.ico"],
    prettyurls=get(ENV, "CI", nothing) == "true",
    size_threshold_ignore = ["man/listfunctions.md"]), # easier local build

    warnonly = Documenter.except(:footnote),
    pages=[
        "Home"      => "index.md",
        "List of functions" => "man/listfunctions.md",
        "Authors" => "man/authors.md",
        "Contributing" => "man/contributing.md",
        "Code of Conduct" => "man/code_of_conduct.md",
        "Security" => "man/security.md",
        "License" => "man/license.md"
    ],
)

withenv("GITHUB_REPOSITORY" => "AnStroh/MOBILE.jl") do
     deploydocs(
         repo = "github.com/AnStroh/MOBILE.jl",
         devbranch = "main",
         push_preview = true
     )
 end
