# Quick Start Guide

## Getting started

After the installation of LaplacianOpt and CPLEX (use GLPK if open-source is preferable) from the Julia package manager, and the cost matrix of the complete graph (for example, see [`examples/instances/5_nodes`](https://github.com/harshangrjn/LaplacianOpt.jl/tree/main/examples/instances/5_nodes)) is provided in the JSON format, an optimization model to maximize the algebraic connectivity of the weighted graph's Laplacian matrix can be executed with a few lines of code as follows:

```julia
using LaplacianOpt
using JuMP
using CPLEX

params = Dict{String, Any}(
"num_nodes" => 5,
"instance" => 1,
"optimizer" => "cplex"
)

lopt_optimizer = JuMP.optimizer_with_attributes(CPLEX.Optimizer) 
results = LaplacianOpt.run_LOpt_model(params, lopt_optimizer)
```

# Extracting results
The run commands (for example, `run_LOpt_model`) in LaplacianOpt return detailed results in the form of a dictionary. This dictionary can be saved for further processing as follows,

```julia
results = LaplacianOpt.run_LOpt_model(params, lopt_optimizer)
```

For example, for the given instance of a complete graph, the algorithm's runtime and the optimal objective value (maximum algebraic connectivity) can be accessed with,

```julia
results["solve_time"]
results["objective"]
```

The `"solution"` field contains detailed information about the solution produced by the optimization model.
For example, one can obtain the edges of the optimal graph toplogy from the symmetric adjacency matrix with,

```Julia
optimal_graph = LaplacianOpt.optimal_graph_edges(results["solution"]["z_var"])
```
and the Fiedler vector of the optimal graph topology from the edge weights and the adjacency matrix with,
```Julia
data = LaplacianOpt.get_data(params)
optimal_adjacency_matrix = results["solution"]["z_var"] .* data["edge_weights"] 
optimal_fiedler_vector = LaplacianOpt.fiedler_vector(optimal_adjacency_matrix)
```

# Visualizing results
LaplacianOpt also currently supports the visualization of optimal graphs layouts obtained from the `results` dictionary (from above. To do so, these are the two options: 
+ [TikzGraphs](https://github.com/JuliaTeX/TikzGraphs.jl) package for a simple and quick visualization of the graph layout without support to include edge weights, which can be executed with 

```julia
data = LaplacianOpt.get_data(params)
LaplacianOpt.visualize_solution(results, data, visualizing_tool = "tikz")
```

+ [Graphviz](https://graphviz.org) package for better visualization of weighted graphs. To this end, LaplacianOpt generates the raw `.dot` file, which can be further visualized using the Graphviz software either via the direct [installation](https://graphviz.org/download/) on the computer or using an online front-end visualization GUI (for example, see [Edotor](https://edotor.net)). Dot files can be generated in LaplacianOpt with 

```julia
data = LaplacianOpt.get_data(params)
LaplacianOpt.visualize_solution(results, data, visualizing_tool = "graphviz")
```
For example, on a weighted complete graph with 10 nodes in [instance #1](https://github.com/harshangrjn/LaplacianOpt.jl/blob/main/examples/instances/10_nodes/10_1.json), the optimal spanning tree with maximum algebraic connectivity, out of ``10^8`` feasible solutions, obtained by LaplacianOpt (using Graphviz visualization) is shown below 
