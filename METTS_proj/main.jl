import tn_julia: svd, contract, updateLeft, identity, permutedims, spinlocalspace, reshape, rightcanonical, leftcanonical, mpo_expectation, applyHtoC, computeLeftEnvironment, computeRightEnvironment, sitecanonical, canonForm
using Plots
include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")



using JLD2

L = 10
beta = 1
steps = 10
ensemble_size = 15
Nkeep = 30
Nsteps = beta*100


savemat = Matrix{Any}(undef, ensemble_size, steps)

function createEnsemble()    
    for e in 1:ensemble_size
        println("starting ensemble $e")
        success = false
        trial = 0
        max_trials = 10

        while !success && trial < max_trials
            try
                cps = cpsrandom(1,L) #new seed for each loop
                cps_xz = copy(cps)

                for k in 1:steps
                    metts = tdmrg(cps_xz,beta,Nsteps,Nkeep) #time evolve to create metts 
                    savemat[e,k] = metts #save metts to row e, column k
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
                println("Error in ensemble $e, trial $trial, retrying")
            end
        end
        @save "output_betaIs1_Lis10.jld2" savemat
    end
end

ensemble = mat = load("output_betais1_Lis10.jld2", "savemat")




energy_arr = zeros(10)  # one sum per column

for j in 1:10
    for i in 1:15
        energy_arr[j] += mpo_expectation(heisenbergmpo(10,1.0),ensemble[i,j])
    end
end

plot_arr = energy_arr ./ (L*ensemble_size)
plot(1:steps, [plot_arr], title="energy per site e vs Steps", label=["Z and X"], linewidth=3)
xlabel!("Step Number")
ylabel!("Energy per Site")


createEnsemble()
