include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")

import LinearAlgebra: norm
import tn_julia: leftcanonical!, rightcanonical!
using LinearAlgebra

# for i in 1:4
#     L=5
#     cps = cps_z(1, 10)
#     metts = ctm(copy(cps), 1, 10)
#     collap = cpscollapse(copy(metts),3)
#     println("Original CPS energy: ", mpo_expectation(heisenbergmpo(10,1.0), cps)/L, "original cps norm ", norm(cps))
#     println("METTS energy: ", mpo_expectation(heisenbergmpo(10,1.0), metts)/L, "metts norm ", norm(metts))
#     println("collapsed energy: ", mpo_expectation(heisenbergmpo(10,1.0), collap)/L, "collap norm ", norm(collap))
# end


L= 7

randmps = cps_z(1,L)

# #Now, evolve it to metts
smallmetts1 = ctm(copy(randmps),0.1,1)
smallmetts20 = ctm(copy(randmps),20,1)

mpo_expectation(heisenbergmpo(L,1.0), smallmetts1)
mpo_expectation(heisenbergmpo(L,1.0), smallmetts20)




# firstcps = cps_z(1,L)
# steps = 10
# beta = 20
# cps = copy(firstcps)
# for i in 1:steps
#     metts = ctm(cps,beta,10)
#     println("ENERGY IS: ",mpo_expectation(heisenbergmpo(L,1.0), metts))
#     cps = cpscollapse(cps,3)
# end
#L=10 so i should get ~ -11
