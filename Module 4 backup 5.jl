### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 79a50702-3f0f-11f1-a497-d5c2f97e9fe3
begin
	using Graphs
	using Karnak
	using Colors
	using StatsBase
	using Plots
	using Distributions
	using GLM
	using DataFrames
	using LinearAlgebra
	using PlutoUI
	using DifferentialEquations
end

# ╔═╡ e4367b45-f1d7-49af-92f8-b3a3c213eb4b
md"""
# Module 4
"""

# ╔═╡ 7273abcd-9fb3-410d-aa7d-dc3bafb2c6c7
begin

# 1. Setup Parameters
N = 100               
K = 0.05              # Reduced K since we removed 1/N scaling
target_degree = 4     
p_rewire = 0.1        

# 2. Network Topology
g = watts_strogatz(N, target_degree, p_rewire)

# 3. System Dynamics
const ω = randn(N)

function kuramoto_unscaled!(dθ, θ, p, t)
    K, g, N = p
    for i in 1:N
        coupling = 0.0
        # More efficient neighbor lookup using Graphs.jl adjacency list
        for j in neighbors(g, i)
            coupling += sin(θ[j] - θ[i])
        end
        dθ[i] = ω[i] + K * coupling
    end
end

# 4. Simulation
θ0 = 2π .* rand(N)
tspan = (0.0, 100.0)
p = (K, g, N)

prob = ODEProblem(kuramoto_unscaled!, θ0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-7, abstol=1e-7)

# 5. Order Parameter Calculation
function compute_order_parameter(sol, N)
    r = [abs(sum(exp.(im .* sol.u[i])) / N) for i in 1:length(sol.t)]
    return r
end

r_values = compute_order_parameter(sol, N)

# 6. Plotting
l = @layout [a; b]
p1 = plot(sol.t, mod2pi.(reduce(hcat, sol.u)'), 
          ylabel="θ (mod 2π)", legend=false, title="Phase Evolution")
p2 = plot(sol.t, r_values, 
          xlabel="Time", ylabel="Order Parameter (r)", 
          title="Synchronization Analysis (K = $K)", ylims=(0,1))

plot(p1, p2, layout=l)
end

# ╔═╡ Cell order:
# ╟─e4367b45-f1d7-49af-92f8-b3a3c213eb4b
# ╠═79a50702-3f0f-11f1-a497-d5c2f97e9fe3
# ╠═7273abcd-9fb3-410d-aa7d-dc3bafb2c6c7
