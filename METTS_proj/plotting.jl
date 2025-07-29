import tn_julia: svd, contract, updateLeft, identity, permutedims, spinlocalspace, reshape, rightcanonical, leftcanonical, mpo_expectation, applyHtoC, computeLeftEnvironment, computeRightEnvironment, sitecanonical, canonForm
using Plots
include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")

tol = 1e-13

L = 5
beta = 5
steps = 10
ensemble_size = 5
Nkeep = 30
Nsteps = beta*50

# initialise energy array for Sz Sx 
E_zx = zeros(steps)

# initialise energy array for only Sz
E_z = zeros(steps)

for j in 1:ensemble_size
    cps = cpsrandom(1, L)
    c = copy(cps)
    cz = copy(cps)
    print("STARTING IS: ", mpo_expectation(heisenbergmpo(L, 1.0), copy(c)))

    #iteration for Sz Sx random
    for i in 1:steps
        metts = tdmrg(c, beta, Nsteps, Nkeep)

        # measure energy E for metts generated here
        E = mpo_expectation(heisenbergmpo(L, 1.0), copy(metts))
        println(E)
        E_zx[i] = E_zx[i] + E
        #println("E: ", E)
        rolls = rand()
        if rolls > 0.5
            c = cpscollapse(metts, 1)
        else 
            c = cpscollapse(metts, 3)
        end
        println("Sz Sx ensemble state number, step number: ", j, ", ", i)
    end

    #     for i in 1:steps
    #     metts = tdmrg(c, beta, Nsteps, Nkeep)

    #     # measure energy E for metts generated here
    #     E = mpo_expectation(heisenbergmpo(L, 1.0), copy(metts))
    #     E_z[i] = E_z[i] + E
    #     println("E: ", E)

    #     c = cpscollapse(metts, 3)
    #     println("Sz ensemble state number, step number: ", j, ", ", i)
    # end

end


# to get energy per state per site
zx_energy = E_zx / (ensemble_size * L) 
z_energy = E_z / (ensemble_size * L) 


# plotting the energies (change title later)
plot(1:steps, [zx_energy], title="energy per site e vs Steps", label=["Z and X" "Z"], linewidth=3)
xlabel!("Step Number")
ylabel!("Energy per Site")
print(zx_energy)


