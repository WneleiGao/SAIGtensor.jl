using Documenter, SAIGtensor

makedocs(
    modules = [SAIGtensor],
    doctest  = true)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/WneleiGao/SAIGtensor.jl.git",
    julia  = "0.5",
    osname = "osx")
