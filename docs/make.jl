import Documenter
using QuantumCircuitOpt
using DocumenterTools: Themes

for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "src/assets/themes/style.scss"), String)
    theme = read(joinpath(@__DIR__, "src/assets/themes/$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "src/assets/themes/$(w).scss"), header*"\n"*theme)
end

Themes.compile(
    joinpath(@__DIR__, "src/assets/themes/documenter-light.css"),
    joinpath(@__DIR__, "src/assets/themes/light.scss"),
)
Themes.compile(
    joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"),
    joinpath(@__DIR__, "src/assets/themes/dark.scss"),
)

Documenter.makedocs(
    modules = [QuantumCircuitOpt],
    sitename = "QuantumCircuitOpt.jl",
    format = Documenter.HTML(mathengine = Documenter.MathJax(), 
                             assets = [asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
                                             "assets/extra_styles.css"],
                             prettyurls = get(ENV, "CI", nothing) == "true",
                             sidebar_sitename=false
                             ),
    # strict = true,
    authors = "Harsha Nagarajan",
    # checkdocs = :none,
    pages = [
        "Introduction" => "index.md",
        "Quick Start guide" => "quickguide.md",
        "Quantum Gates Library" => [
            "1-qubit gates"     => "1_qubit_gates.md",
            "2-qubit gates"     => "2_qubit_gates.md",
            "3-qubit gates"     => "3_qubit_gates.md",
            "Multi-qubit gates" => "multi_qubit_gates.md",
        ],
        "Function References" => "function_references.md",
    ],
)

Documenter.deploydocs(
    repo = "github.com/harshangrjn/QuantumCircuitOpt.jl.git",
    push_preview = true
)
