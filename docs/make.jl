using Documenter, SAIGtensor

makedocs(
    modules = [SAIGtensor],
    format  = :html,
    sitename = "SAIGtensor.jl",
    doctest  = true,
    strict   = true
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/WneleiGao/SAIGtensor.jl.git",
    julia  = "release",
    make   = nothing,
    target = "build"
)
