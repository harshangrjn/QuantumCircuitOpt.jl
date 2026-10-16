using Documenter
using QuantumCircuitOpt

makedocs(
    modules = [QuantumCircuitOpt],
    sitename = "QuantumCircuitOpt.jl",
    checkdocs = :exports,
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        assets = ["assets/extra_styles.css", "assets/publications.js"],
        prettyurls = get(ENV, "CI", "false") == "true",
        sidebar_sitename = false,
    ),
    authors = "Harsha Nagarajan",
    pages = [
        "Introduction" => "index.md",
        "Quick Start guide" => "quickguide.md",
        "Quantum Gates Library" => [
            "1-qubit gates" => "1_qubit_gates.md",
            "2-qubit gates" => "2_qubit_gates.md",
            "3-qubit gates" => "3_qubit_gates.md",
            "Multi-qubit gates" => "multi_qubit_gates.md",
        ],
        "Function References" => "function_references.md",
        "Research & Publications" => "research.md",
    ],
)

deploydocs(
    repo = "github.com/harshangrjn/QuantumCircuitOpt.jl.git",
    devbranch = "master",
    push_preview = false,
)
