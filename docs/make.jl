using ModelingToolkitSampledData
using ModelingToolkit
using Documenter
using Plots
gr(fmt=:png)

ENV["JULIA_DEBUG"]=Documenter # Enable this for debugging
ENV["DOCS_BUILD"] = true # used to lower the default frame rate in animations for the docs

DocMeta.setdocmeta!(ModelingToolkitSampledData, :DocTestSetup, :(using ModelingToolkitSampledData); recursive = true)

makedocs(;
         modules = [ModelingToolkitSampledData],
        #  remotes = Dict(
        #     dirname(dirname(pathof(ModelingToolkit))) => (Remotes.GitHub("SciML", "ModelingToolkit.jl"), "9"),
        #     dirname(dirname(pathof(ModelingToolkitStandardLibrary))) => (Remotes.GitHub("SciML", "ModelingToolkitStandardLibrary.jl"), "2"),
        # ),
         authors = "Fredrik Bagge Carlson, JuliaHub Inc.",
         #  strict = [:example_block, :setup_block, :eval_block],
         sitename = "ModelingToolkitSampledData.jl",
         warnonly = [:missing_docs, :cross_references, :docs_block],
         pagesonly = true,
        #  draft = true,
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  edit_link = nothing),
         pages = [
             "Home" => "index.md",
             "Tutorials" => [
                "Getting started with Sampled-Data Systems" => "tutorials/SampledData.md",
                "Noise" => "tutorials/noise.md",
             ],
             "Examples" => [
                "Controlled DC motor" => "examples/dc_motor_pi.md",
             ],
             "Block library" => "blocks.md",
             "Functions" => "api.md",
         ],
)

deploydocs(;
           repo = "github.com/JuliaComputing/ModelingToolkitSampledData.jl",
           devbranch = "main")
