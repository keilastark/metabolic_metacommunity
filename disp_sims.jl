# Dispersal simulations
#cd("/data/landscape_data")

#script_dir = dirname(@__FILE__)

# Navigate to the data directory
#cd(joinpath(script_dir, "landscape_data"))

using Pkg

# Activate the current directory
Pkg.activate(".")

# Add necessary packages
packages = [
    "Distributions",
    "Distances",
    "DataFrames",
    "CSV",
    "StatsBase",
    "DataStructures",
    "ProgressMeter",
    "Distributed",
    "LinearAlgebra",
    "Random"
]

for pkg in packages
    Pkg.add(pkg)
end

# Import required packages
using Distributions
using Distances
using DataFrames
using CSV
using StatsBase
using DataStructures
using ProgressMeter
using Distributed
using LinearAlgebra
using Random

reps = 30
warming_V = [0, 2.5, 5, 7.5, 10]
r = 10
d_temp_V = ["unimodal", "monotonic", "independent"]
d_V = exp.(range(log(1e-5), stop = log(1), length = 3))
topt_V = [20,22.5,25,27.5,30]

@showprogress for b in eachindex(warming_V)
    warming_level = warming_V[b]
    warming_write = string(warming_level, "_C")

    @showprogress for rep = 1:reps
        env_df = CSV.read("./env_$rep.csv", DataFrame)
        env_df.env1 .+= warming_level
        disp_mat = CSV.read("./disp_mat_$rep.csv", DataFrame)
        disp_mat = Matrix(disp_mat)
        disp_mat = disp_mat[:, 2:end]
        landscape = CSV.read("./landscape_$rep.csv", DataFrame)
        time = CSV.read("./time_$rep.csv", DataFrame)

        S = 50 # species
        M = size(disp_mat)[1] # patches
        d = 0.5
        α = rand(Uniform(0, 1.0), S, S)
        α[diagind(α)] = repeat([1.0], outer = S)
        α = α * 0.01

        env_mat2 = zeros(maximum(env_df.time), M) # pulling temperature value for each time point and patch
        for patch in 1:50
            patch_indices = (env_df.patch .== patch)
            patch_values = env_df.env1[patch_indices]
            env_mat2[:, patch] = patch_values
        end

        burn_in = time.burn_in[1]
        generations = maximum(env_df.time) - burn_in
        Tmax = burn_in + generations

        for i in eachindex(d_temp_V)
            d_temp_function = d_temp_V[i]
            

            for j in eachindex(topt_V)
                topt = topt_V[j]


                for disp in eachindex(d_V)
                    d = d_V[disp]
                    global Model_df
                    N = rand(Poisson(2), M, S) * 1.0 # initializing at low density
                    sampV = convert(Array{Int64}, collect((burn_in+800):20:Tmax))
                    N_save = N
                    seedV = convert(
                        Array{Int64},
                        collect(burn_in/(10):burn_in/(10):burn_in/2), # subsample
                    )

                    @showprogress for gen = 1:Tmax
                        x = mean(env_mat2[gen, :])
                        # Random values for temperature-dependence parameters
                        exp_matrix = rand(Uniform(0.5,0.7), M, S)

                        # Generate disp_temp matrix
                        if d_temp_function == "unimodal"
                            disp_temp = d .* exp.(-exp_matrix ./ 8.62e-05 .* (1 / 293.15 .- 1 ./ (x .+ 273.15))) .* (1 ./ (1 .+ (exp_matrix ./ (2 .- exp_matrix)) .* exp.(2 ./ 8.62e-05 .* (1 ./ (topt .+ 273.15) .- 1 ./ (x .+ 273.15)))))
                        elseif d_temp_function == "monotonic"
                            disp_temp = d .* exp.(exp_matrix ./ 8.62e-05 .* (1 / 293.15 .- 1 ./ (x .+ 273.15)))
                        elseif d_temp_function == "independent"
                            disp_temp = fill(d, M, S)
                        end

                        # Ensure values are within bounds
                        disp_temp[disp_temp .< 0] .= 0
                        disp_temp[disp_temp .> 1] .= 1

                        r_temp = r * exp(0.65 / 8.62e-05 * (1 / 293.15 - 1 / (x + 273.15))) * (1 / (1 + (0.65 / (2 - 0.65)) * exp(2 / 8.62e-05 * (1 / (topt + 273.15) - 1 / (x + 273.15)))))

                        density = N .* α
                        λ_v = r_temp .* N .* (1.0 ./ (1.0 .+ density))
                        λ_v[λ_v .< 0.0] .= 0.0

                        N = [rand(Poisson(λ)) for λ in λ_v]
                        emigrants = [rand(Binomial(n, disp_temp[p])) for (p, n) in enumerate(N)]
                        immigrants_exp = disp_mat .* disp_temp .* emigrants # expected number of immigrants
                        immigrants_S = sum(emigrants, dims = 1) # number of immigrants per species
                        immigrants = zeros(Int, M, S)
                        for l in 1:S
                            probabilities = zeros(size(immigrants_exp[:, l]))
                            denominator = sum(immigrants_exp[:, l])
                            if denominator != 0
                                probabilities .= immigrants_exp[:, l] ./ denominator
                            end
                            sampled_patches = sample(1:M, Weights(probabilities), immigrants_S[l], replace=true)
                            count_map = countmap(sampled_patches)
                            sorted_patches = sort(collect(keys(count_map)))

                            for key in sorted_patches
                                immigrants[key, l] = count_map[key]
                            end
                        end

                        sum(emigrants, dims = 1)
                        sum(immigrants, dims = 1)
                        N = N .- emigrants .+ immigrants
                        N[N .< 0.0] .= 0.0

                        if any(y -> y == gen, sampV)
                            N_save = cat(dims = 3, N_save, N)
                        end
                        N = N .* 1.0
                    end

                    N_save = N_save[:, :, 2:end]

                    Model_df_1 = DataFrame(
                        N = N_save[:],
                        Species = repeat(1:S, inner = M, outer = length(sampV)),
                        Time = repeat(1:length(sampV), inner = S * M),
                        Patch = repeat(1:M, outer = length(sampV) * S),
                        topt = topt,
                        d_temp_func = string(d_temp_function),
                        warming = string(warming_write)
                    )

                    Model_df_1 = Model_df_1[Model_df_1[:, :N] .> 0, :]

                    if topt == topt_V[1] && d_temp_function == d_temp_V[1] 
                        Model_df = Model_df_1
                    else
                        Model_df = [Model_df; Model_df_1]
                    end
                end
            end
        end
        CSV.write("./outputs/disp_outputfile_$(warming_write)_$(rep).csv", Model_df)
    end
end
