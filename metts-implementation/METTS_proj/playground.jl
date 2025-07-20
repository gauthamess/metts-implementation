include("basics.jl")
include("states.jl")
include("gates.jl")
include("main.jl")

function normalise(mps)
    L = length(mps)
    Id = proj_mpo(reshape(LinearAlgebra.I(3), 1, 3, 1, 3), L, 2)
    norm = mpo_expectation(Id, mps)
    mps[1] = mps[1] / sqrt(norm)
    return mps
end



cps = cps_z(1,10)
out = ctm(copy(cps), 2, 10)  # Use beta=1.0, not 0
println("Original CPS energy: ", mpo_expectation(heisenbergmpo(10,1.0), cps))
println("METTS energy: ", mpo_expectation(heisenbergmpo(10,1.0), out))

cps = normalise(cps_z(1, 10))
out = normalise(ctm(copy(cps), 5, 10))
println("Original CPS energy: ", mpo_expectation(heisenbergmpo(10,1.0), cps))
println("METTS energy: ", mpo_expectation(heisenbergmpo(10,1.0), out))

for β in [80]
    metts = ctm(copy(cps), β, 10)
    metts = normalise(metts)  # Make sure it's normalized
    E = mpo_expectation(heisenbergmpo(10,1.0), metts)
    println("β = $β: E = $E")
end

metts = ctm(copy(cps), 10, 30)

E = mpo_expectation(heisenbergmpo(10,1.0), cps)

Id = proj_mpo(reshape(LinearAlgebra.I(3), 1, 3, 1, 3), 10, 2)

mpo_expectation(Id,metts)
metts=normalise(metts)
mpo_expectation(Id,metts)