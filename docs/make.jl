using Documenter
using QuantumCircuitOpt
using DocumenterTools: Themes

for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "src/assets/themes/style.scss"), String)
    theme = read(joinpath(@__DIR__, "src/assets/themes/$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "src/assets/themes/$(w).scss"), header*"\n"*theme)
end

Themes.compile(joinpath(@__DIR__, "src/assets/themes/light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "src/assets/themes/dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

makedocs(
    modules = [QuantumCircuitOpt],
    sitename = "QuantumCircuitOpt",
    format = Documenter.HTML(mathengine = Documenter.MathJax(), 
                             assets = [asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css)],
                             prettyurls = get(ENV, "CI", nothing) == "true"),
    strict = true,
    authors = "Harsha Nagarajan",
    pages = [
        "Introduction" => "index.md",
        # "Quick Start guide" => "quickguide.md",
        "Function References" => "function_references.md",
    ],
)

deploydocs(
    repo = "github.com/harshangrjn/QuantumCircuitOpt.jl.git",
    push_preview = true
)
