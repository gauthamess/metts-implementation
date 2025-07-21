import tn_julia: svd, contract, updateLeft, identity, permutedims, spinlocalspace, reshape, rightcanonical, leftcanonical, mpo_expectation, applyHtoC, computeLeftEnvironment, computeRightEnvironment, sitecanonical, canonForm
using Plots
include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")

L = 5
S = 1
beta = 1
steps = 5
ensemble_size = 10
Nkeep = 8



# initialise energy array for Sz Sx 
E_zx = zeros(steps)

# initialise energy array for only Sz
E_z = zeros(steps)

for j in 1:ensemble_size
    cps = cps_z(1, L)
    c = copy(cps)
    cz = copy(cps)

    # # iteration for Sz Sx
    # for i in 1:steps
    #     metts = ctm(c, beta, Nkeep)

    #     # measure energy E for metts generated here
    #     E = mpo_expectation(heisenbergmpo(L, 1.0), copy(metts))
    #     E_zx[i] = E_zx[i] + E
    #     println("E: ", E)
    #     if isodd(i)
    #         c = cpscollapse(metts, 1)
    #     else 
    #         c = cpscollapse(metts, 3)
    #     end
    #     println("Sz Sx ensemble state number, step number: ", j, ", ", i)
    # end


    # iteration for Sz
    for i in 1:steps
        metts = ctm(cz, beta, Nkeep)

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

# cpscollapase(testerarrz[4],3)
# size(testerarrzx)
# testerarrzx[4][28][1,:,1]
# size(testerarrzx[6][10])
# input = testerarrzx[1]
# check = ctm(testerarrzx[1],1,15)
# check[9][1,:,1]

# plotting the energies (change title later)
plot(1:steps, [z_energy zx_energy], title="energy per site e vs Steps", label=["Z only" "Z and X"], linewidth=3)
xlabel!("Step Number")
ylabel!("Energy per Site")
print(zx_energy)