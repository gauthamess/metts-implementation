import tn_julia: svd, contract, updateLeft, identity, permutedims, spinlocalspace, reshape, rightcanonical, leftcanonical, mpo_expectation, applyHtoC, computeLeftEnvironment, computeRightEnvironment, sitecanonical, canonForm
using Plots
using JLD2
using Statistics
include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")

tol = 1e-10




L = 40
beta = 1
steps = 10
ensemble_size = 15
Nkeep = 30
Nsteps = beta*50

        


function createEnsemble()    
    for e in 9:ensemble_size
        println("starting ensemble $e")
        savemat_xz = Vector{Any}(undef, steps)
        savemat_z = Vector{Any}(undef, steps)
        success = false
        trial = 0
        max_trials = 10

        while !success && trial < max_trials
            try
                cps = cpsrandom(1,L) #new seed for each loop
                cps_xz = deepcopy(cps)
                cps_z = deepcopy(cps)
                for k in 1:steps
                    metts_xz = tdmrg(cps_xz,beta,Nsteps,Nkeep) #time evolve to create metts 
                    metts_z = tdmrg(cps_z,beta,Nsteps,Nkeep) #time evolve to create metts 
                    savemat_xz[k] = deepcopy(metts_xz) #save metts to row e, column k
                    savemat_z[k] = deepcopy(metts_z) #save metts to row e, column k
                    println("saved ensemble $e, metts $k")
                    rolls = rand()
                    if rolls > 0.5 #random collapse is new cps
                        cps_xz = cpscollapse(metts_xz,1)
                    else
                        cps_xz = cpscollapse(metts_xz,3)
                    end
                    cps_z = cpscollapse(metts_z,3) #collapse to z)
                end
                success = true
            catch err
                trial += 1
                println("Error: ", err)
                println("Error in ensemble_$e, trial $trial, retrying")
            end
        end
        @save "first_plot_x_$(beta)_LIs$(L)_$e.jld2" savemat_z
        @save "first_plot_xz_$(beta)_LIs$(L)_$e.jld2" savemat_xz

    end
end
createEnsemble()


# TO LOAD BETA = 1 ONLY!!:
# metts_matrix = Matrix{Any}(undef, ensemble_size, steps)
# filename = "output_betaIs1_LIs50.jld2" 
# @load filename savemat
# metts_matrix = savemat

metts_matrix_z = Matrix{Any}(undef, ensemble_size, steps)
metts_matrix_xz = Matrix{Any}(undef, ensemble_size, steps)
for e in 1:ensemble_size
    filename = "first_plot_x_$(beta)_LIs$(L)_$e.jld2" 
    @load filename savemat_z
    for k in 1:steps
        metts_matrix_z[e, k] = savemat_z[k]
    end
end
for e in 1:ensemble_size
    filename = "first_plot_xz_$(beta)_LIs$(L)_$e.jld2" 
    @load filename savemat_xz
    for k in 1:steps
        metts_matrix_xz[e, k] = savemat_xz[k]
    end
end


sz_arr = zeros(ensemble_size)
sz2_arr = zeros(ensemble_size)
mpo_sz = sz_tot_mpo(L)
mpo_sz2 = sz_squared_mpo(L)

for e in 1:ensemble_size
    sz_arr[e] = mpo_expectation(mpo_sz,metts_matrix[e,steps])
    sz2_arr[e] = mpo_expectation(mpo_sz2,metts_matrix[e,steps])
    println(sz2_arr[e], " is sz2 expectation for ensemble ", e)
end

sz_avg = sum(sz_arr)/ensemble_size
sz2_avg = sum(sz2_arr)/ensemble_size

susceptibility = beta*(sz2_avg - sz_avg^2)/L

