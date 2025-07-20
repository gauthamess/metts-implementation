include("basics.jl")
include("states.jl")
include("gates.jl")
include("main.jl")

cps = normalise(cps_z(1, 10))
out = normalise(ctm(copy(cps), 5, 10))
out2 = normalise(cpscollapse(copy(out),3))
println("Original CPS energy: ", mpo_expectation(heisenbergmpo(10,1.0), cps))
println("METTS energy: ", mpo_expectation(heisenbergmpo(10,1.0), out))
println("collapsed energy: ", mpo_expectation(heisenbergmpo(10,1.0), out2))




L = 10
steps = 40
beta = 1

cps = cps_z(1, L)
    # iteration for Sz
cz = copy(cps)

E = zeros(steps)

for i in 1:steps
    metts = ctm(cz, beta, 10) #get metts

    E[i] = mpo_expectation(heisenbergmpo(L, 1.0), metts) #energy of metts
    println(E[i])

    cz = cpscollapse(metts, 3) #new cps = collapsed metts
end


beta = 1

cps = cps_z(1,10)
metts1 = ctm(cps, 1, 15)
E1 = mpo_expectation(heisenbergmpo(10,1.0), metts1)



cps2 = cpscollapse(metts1, 3) #3 means Sz  
metts2 = ctm(cps2, beta, 15)
E2 = mpo_expectation(heisenbergmpo(10,1.0),metts2)

println("Energy before collapse: $E1")
println("Energy after collapse: $E2")
println("ΔE = $(E2-E1)")

#AVERAGE FIRST COLLAPSE Energy
ensemblesize = 5
E = zeros(5)
for i in 1:ensemblesize
    println(i)
    cps = cps_z(1,20)
    metts = ctm(cps,1,20)
    energy = mpo_expectation(heisenbergmpo(20,1.0),metts)
    E[i] = energy
    println(E)
    energy=0
end
sum(E)/100