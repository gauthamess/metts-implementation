include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")
import LinearAlgebra: norm

for i in 1:4
    L=5
    cps = cps_z(1, 10)
    metts = ctm(copy(cps), 1, 10)
    collap = cpscollapse(copy(metts),3)
    println("Original CPS energy: ", mpo_expectation(heisenbergmpo(10,1.0), cps)/L, "original cps norm ", norm(cps))
    println("METTS energy: ", mpo_expectation(heisenbergmpo(10,1.0), metts)/L, "metts norm ", norm(metts))
    println("collapsed energy: ", mpo_expectation(heisenbergmpo(10,1.0), collap)/L, "collap norm ", norm(collap))
end

vec = [1,2,3]
norm(vec)

m = 0

cps = cps_z(1,10, base = 1)
cps[3][1,:,1]

contract(cps[3][1,:,1],1,cps[3][1,:,1],1)



newcps = copy(cps)

newcps[2][1,:,1]
projector = proj(-1,3)

mpo = proj_mpo(projector,10,2)

applyHtoC(mpo,newcps,2)