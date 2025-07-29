import tn_julia: svd, contract, updateLeft, identity, permutedims, spinlocalspace, reshape, rightcanonical, leftcanonical, mpo_expectation, applyHtoC, computeLeftEnvironment, computeRightEnvironment, sitecanonical, canonForm
using Plots
include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")

tol = 1e-10


using JLD2

L = 100
beta = 2
steps = 10
ensemble_size = 5
Nkeep = 25
Nsteps = beta*50

        


function createEnsemble()    
    for e in 1:ensemble_size
        println("starting ensemble $e")
        savemat = Vector{Any}(undef, steps)
        success = false
        trial = 0
        max_trials = 10

        while !success && trial < max_trials
            try
                cps = cpsrandom(1,L) #new seed for each loop
                cps_xz = deepcopy(cps)

                for k in 1:steps
                    metts = tdmrg(cps_xz,beta,Nsteps,Nkeep) #time evolve to create metts 
                    savemat[k] = deepcopy(metts) #save metts to row e, column k
                    println("saved ensemble $e, metts $k")
                    rolls = rand()
                    if rolls > 0.5 #random collapse is new cps
                        cps_xz = cpscollapse(metts,1)
                    else
                        cps_xz = cpscollapse(metts,3)
                    end
                end
                success = true
            catch err
                trial += 1
                println("Error: ", err)
                println("Error in ensemble_$e, trial $trial, retrying")
            end
        end
        @save "output_betaIs$(beta)_LIs_$(L)_$e.jld2" savemat
    end
end
createEnsemble()


#TO LOAD BETA = 1 ONLY!!:
# metts_matrix = Matrix{Any}(undef, ensemble_size, steps)
# filename = "output_betaIs1_LIs50.jld2" 
# @load filename savemat
# metts_matrix = savemat

metts_matrix = Matrix{Any}(undef, ensemble_size, steps)
for e in 1:ensemble_size
    filename = "output_betaIs$(beta)_LIs_$(L)_$e.jld2"
    @load filename savemat
    for k in 1:steps
        metts_matrix[e, k] = savemat[k]
    end
end


sz_arr = zeros(ensemble_size)
sz2_arr = zeros(ensemble_size)
mpo_sz = sz_tot_mpo(L)
mpo_sz2 = sz_squared_mpo(L)

for e in 1:ensemble_size
    sz_arr[e] = mpo_expectation(mpo_sz,metts_matrix[e,steps])
    sz2_arr[e] = mpo_expectation(mpo_sz2,metts_matrix[e,steps])
end

sz_avg = sum(sz_arr)/ensemble_size
sz2_avg = sum(sz2_arr)/ensemble_size

susceptibility = beta*(sz2_avg - sz_avg^2)/L

l=3
testcps = cps_z(1,l)

mpo_expectation(sz_tot_mpo(l), testcps)
mpo_expectation(sz_squared_mpo(l), testcps)




Sz = Sxyz(1, 3)
eigvals(Sz)
eigvecs(Sz)
