### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 28c123be-2e90-11f1-a502-1d7d2be52c84
begin
	using Plots
	using Random
	using PlutoUI
end

# ╔═╡ e871d352-6d18-4964-a92e-ca197a928988
function cml_step(x, r, ϵ, f)
    N = length(x)
    new_x = similar(x)
    
    fx = f.(x, r) 
    
    for i in 1:N
        left  = mod1(i - 1, N)
        right = mod1(i + 1, N)
        
        new_x[i] = (1 - ϵ) * fx[i] + (ϵ / 2) * (fx[left] + fx[right]) % 1
    end
    return new_x
end

# ╔═╡ 2aad3479-9337-4d52-a941-62b1659caa0f
begin
	function evolve(steps, r, ϵ, map, N)
		results = zeros(steps, N)
		results[1, :] = rand(Float64, N) 

		for t in 2:steps
		    results[t, :] = cml_step(results[t-1, :], r, ϵ, map)
		end

		return results
	end

	logistic_results = evolve(1500, 1.7522, 0.00115, f_logistic, 100)
	
	heatmap(logistic_results, xlabel="Lattice site", ylabel="Time", title="CML Evolution", c = :deep)
end

# ╔═╡ 14163eda-fda8-45dc-a50f-74bade6cc902
function(θ, Ω) = θ + Ω - sin(2*π*θ)/(2*π)

# ╔═╡ Cell order:
# ╠═28c123be-2e90-11f1-a502-1d7d2be52c84
# ╠═e871d352-6d18-4964-a92e-ca197a928988
# ╠═2aad3479-9337-4d52-a941-62b1659caa0f
# ╠═14163eda-fda8-45dc-a50f-74bade6cc902
