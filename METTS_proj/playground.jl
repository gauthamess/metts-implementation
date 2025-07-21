include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")
include("main.jl")


for i in 1:10
    L=10
    cps = cps_z(1, 10)
    out = ctm(copy(cps), 1, 10)
    out2 = normalise(cpscollapsejoan(copy(out),3))
    println("Original CPS energy: ", mpo_expectation(heisenbergmpo(10,1.0), cps)/L)
    println("METTS energy: ", mpo_expectation(heisenbergmpo(10,1.0), out)/L)
    println("collapsed energy: ", mpo_expectation(heisenbergmpo(10,1.0), out2)/L)
end

#TESTING COLLAPSE VALS
for i in 1:3
    cps = cps_z(1, L)
    energy_before = mpo_expectation(heisenbergmpo(L,1.0), cps)/L
    metts = ctm(copy(cps), 1, L)
    energy_metts = mpo_expectation(heisenbergmpo(L,1.0), metts)/L
    collapsed = cpscollapse(metts, 3)  # Z collapse
    collapsedx = cpscollapse(metts,1)
    energy_afterz = mpo_expectation(heisenbergmpo(L,1.0), collapsed)/L
    energy_afterx = mpo_expectation(heisenbergmpo(L,1.0), collapsed)/L
    println("Before: $energy_before, METTS: $energy_metts,  Afterx: $energy_afterx, Afterz: $energy_afterz")
end






L = 10
steps = 40
beta = 1


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




