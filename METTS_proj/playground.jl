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


L= 4

randmps = cps_z(1,L)

# #Now, evolve it to metts
smallmetts1 = ctm(copy(randmps),1,30)
smallmetts20 = ctm(copy(randmps),10,30)

mpo_expectation(heisenbergmpo(L,1.0),randmps)
mpo_expectation(heisenbergmpo(L,1.0), smallmetts1)
mpo_expectation(heisenbergmpo(L,1.0), smallmetts20)




for i in 10:-1:1
    print(i)
end