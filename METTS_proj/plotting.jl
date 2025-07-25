import tn_julia: svd, contract, updateLeft, identity, permutedims, spinlocalspace, reshape, rightcanonical, leftcanonical, mpo_expectation, applyHtoC, computeLeftEnvironment, computeRightEnvironment, sitecanonical, canonForm
using Plots
include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")

L = 5
beta = 1
steps = 10
ensemble_size = 10
Nkeep = 20
Nsteps = 200


# initialise energy array for Sz Sx 
E_zx = zeros(steps)

# initialise energy array for only Sz
E_z = zeros(steps)

for j in 1:ensemble_size
    cps = cpsrandom(1, L)
    c = copy(cps)
    cz = copy(cps)

    #iteration for Sz Sx random
    for i in 1:steps
        metts = tdmrg(c, beta, Nsteps, Nkeep)

        # measure energy E for metts generated here
        E = mpo_expectation(heisenbergmpo(L, 1.0), copy(metts))
        E_zx[i] = E_zx[i] + E
        println("E: ", E)
        rolls = rand()
        if rolls > 0.5
            c = cpscollapse(metts, 1)
        else 
            c = cpscollapse(metts, 3)
        end
        println("Sz Sx ensemble state number, step number: ", j, ", ", i)
    end


    #iteration for Sz
    for i in 1:steps
        metts = tdmrg(cz, beta, Nkeep,Nsteps)

        # measure energy E for metts generated here
        Ez = mpo_expectation(heisenbergmpo(L, 1.0), copy(metts))
        E_z[i] = E_z[i] + Ez
        println("E: ", Ez)

        cz = cpscollapse(metts, 3)
        
        println("Sz ensemble state number, step number: ", j, ", ", i)
    end
end


# to get energy per state per site
zx_energy = E_zx / (ensemble_size * L) 
z_energy = E_z / (ensemble_size * L)

# plotting the energies (change title later)
plot(1:steps, [z_energy zx_energy], title="energy per site e vs Steps", label=["Z only" "Z and X"], linewidth=3)
xlabel!("Step Number")
ylabel!("Energy per Site")
print(zx_energy)


