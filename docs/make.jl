using Documenter, SAIGtensor

makedocs(
    modules = [SAIGtensor],
    clean = false,
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/WneleiGao/SAIGtensor.jl.git",
    julia  = "0.5",
    osname = "linux"
)
