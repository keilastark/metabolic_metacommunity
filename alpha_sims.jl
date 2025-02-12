#Alpha simulations

#cd("/Users/keilastark/Documents/metabolic_metacommunity/julia_workflow/landscape_data")


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
α_temp_V = ["exponential_epsilon","unimodal_epsilon","invariant"]
α_V = [0.5, 1.5] #[0, 0.5, 1, 1.5] #competition vector
topt_V = [20,22.5,25,27.5,30]
d = 0.0001

#N_save = Array{Float64}(undef, M, S, length(sampV))

## Topt loop

@showprogress for b in eachindex(warming_V)
    warming_level= warming_V[b]
    warming_write=string(warming_level,"_C")

    @showprogress for rep = 1:reps
        env_df = CSV.read("./env_$rep.csv", DataFrame)
        env_df.env1 .+= warming_level
        disp_mat = CSV.read("./disp_mat_$rep.csv",DataFrame)
        disp_mat = Matrix(disp_mat)
        disp_mat = disp_mat[:, 2:end]
        landscape = CSV.read("./landscape_$rep.csv", DataFrame)
        time = CSV.read("./time_$rep.csv", DataFrame)

        S = 30 #species
        #dominants = trunc(Int, round(S * 0.3)) #for comp-col scenario
        M = size(disp_mat)[1] #patches
        
        env_mat2 = zeros(maximum(env_df.time), M) #pulling temperature value for each time point and patch
        for patch in 1:50
            patch_indices = (env_df.patch .== patch)
            patch_values = env_df.env1[patch_indices]
            env_mat2[:, patch] = patch_values
        end

        burn_in = time.burn_in[1]
        generations = maximum(env_df.time) - burn_in
        Tmax = burn_in + generations
 
        for k = 1:(length(α_V)) 
            #α = rand(Normal(0.07,0.01),S,S)#rand(LogNormal(0.001, 0.1), S,S)
            α = rand(Uniform(0, α_V[k]), S, S)
            α_write = string(α_V[k])
            α[diagind(α)] = repeat([1.0], outer = S)
            α = α * 0.01

            for l in eachindex(α_temp_V)

                α_temp_function = α_temp_V[l]

                #for i in eachindex(d_temp_V)
                   # d_temp_function = d_temp_V[i]
                    
                    #for disp in eachindex(d_V)
                        #d = d_V[disp]

                        for j in eachindex(topt_V)
                            global Model_df

                            topt = topt_V[j] #species' topts relative to ambient temperatures (20-30C)

                            N = rand(Poisson(2), M, S) * 1.0 #initializing at low density   
                            sampV =convert(Array{Int64}, collect((burn_in+800):20:Tmax))
                            N_save = N
                            seedV = convert(
                            Array{Int64},
                            collect(burn_in/(10):burn_in/(10):burn_in/2), #subsample
                            )
                            
                            @showprogress for gen = 1:Tmax
                                #global N
                                #global N_save
                                #global λ_save
                                #global env_save
                                #global env_match_save
                                #global den_save
                                x = mean(env_mat2[gen, :])
                                #k = 8.62e-05
                                #tref = 283.15
                                #e = 0.65
                                #eh = 2
                                if (any(y -> y == gen, seedV))
                                    N = N + rand(Poisson(0.5), M, S) * 1.0
                                end

                                if α_temp_function == "exponential_epsilon"
                                    α_temp  =   α .* exp(-0.8/8.62e-05 * (1/293.15 - 1/(x+273.15)))
                                elseif α_temp_function == "unimodal_epsilon"
                                    α_temp =  α .* exp(0.8/8.62e-05 * (1/293.15 - 1/(x+273.15))) * (1/(1+(0.3/(2-0.3))*exp(2/8.62e-05*(1/((topt-2)+273.15)-1/(x+273.15)))))
                                elseif α_temp_function == "invariant" 
                                    α_temp = α
                                end
                                #println(α_temp)
                                r_temp =  r * exp(0.65/8.62e-05 * (1/293.15 - 1/(x+273.15))) * (1/(1+(0.65/(2-0.65))*exp(2/8.62e-05*(1/(topt+273.15)-1/(x+273.15)))))
                                #println(r_temp)
                                density = N * α_temp
                                #println(density)

                                λ_v = r_temp * N .* (1.0 ./ (1.0 .+ density)) 
                                λ_v[λ_v.<0.0] .= 0.0

                                N = [rand(Poisson(λ)) for λ in λ_v]
                                emigrants = [rand(Binomial(n, d)) for n in N]
                                immigrants_exp = disp_mat *d* emigrants #expected number of immigrants
                                immigrants_S = sum(emigrants, dims = 1) #number of immigrants per species
                                immigrants = zeros(Int, M, S)
                                for l in 1:S
                                    
                                    probabilities = zeros(size(immigrants_exp[:, l]))
                                    denominator = sum(immigrants_exp[:, l])
                                        if denominator != 0
                                            probabilities .= immigrants_exp[:, l] ./ denominator
                                        end 
                                    sampled_patches = sample(1:M, Weights(probabilities), immigrants_S[l], replace=true)
 
                                    # Count occurrences of each sampled patch
                                    count_map = countmap(sampled_patches)
                                    sorted_patches = sort(collect(keys(count_map)))
                                    
                                    for key in sorted_patches
                                        immigrants[key , l] = count_map[key]
                                    end 
                                end
                                sum(emigrants, dims = 1)
                                sum(immigrants, dims = 1)
                                N = N .- emigrants .+ immigrants
                                N[N.<0.0] .= 0.0

                                if (any(y -> y == gen, sampV))
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
                            dispersal = d,
                            topt = topt,
                            alpha = string(α_write),
                            alpha_temp_func=string(α_temp_function),
                            warming = string(warming_write),
                            )

                            Model_df_1 = Model_df_1[Model_df_1[!,:N].>0, :]

                           if topt == topt_V[1] && k == 1 && α_temp_function == α_temp_V[1]
                               Model_df = Model_df_1
                                else
                                Model_df = [Model_df; Model_df_1]    
                            end
                        end
                
                    #end 
               # end
            end
        end
        CSV.write("./outputs/alpha_outputfile2_$(warming_write)_$(rep).csv", Model_df)
    end
end