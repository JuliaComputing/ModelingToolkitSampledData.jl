using ModelingToolkitSampledData
using Documenter

ENV["JULIA_DEBUG"]=Documenter # Enable this for debugging
ENV["DOCS_BUILD"] = true # used to lower the default frame rate in animations for the docs

DocMeta.setdocmeta!(ModelingToolkitSampledData, :DocTestSetup, :(using ModelingToolkitSampledData); recursive = true)

makedocs(;
         modules = [ModelingToolkitSampledData],
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
                "Sampled-Data Systems" => "SampledData.md",
             ],
             "Examples" => [
             ],
             "Blocks" => "blocks.md",
         ])

deploydocs(;
           repo = "github.com/JuliaComputing/ModelingToolkitSampledData.jl",
           devbranch = "main")
