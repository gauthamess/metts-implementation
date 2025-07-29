import tn_julia: svd, contract, updateLeft, identity, permutedims, spinlocalspace, reshape, rightcanonical, leftcanonical, mpo_expectation, applyHtoC, computeLeftEnvironment, computeRightEnvironment, sitecanonical, canonForm
using Plots
using JLD2
include("basics.jl")
include("states.jl")
include("gates.jl")
include("collapse.jl")

tol = 1e-13

L = 40
beta = 1
steps = 10
ensemble_size = 2
Nkeep = 30
Nsteps = beta*50


metts_matrix_z = Matrix{Any}(undef, ensemble_size, steps)
metts_matrix_xz = Matrix{Any}(undef, ensemble_size, steps)
for e in 1:ensemble_size
    filename = "first_plot_x_$(beta)_LIs$(L)_$e.jld2" 
    @load filename savemat_z
    for k in 1:steps
        metts_matrix_z[e, k] = savemat_z[k]
    end
end
for e in 1:ensemble_size
    filename = "first_plot_xz_$(beta)_LIs$(L)_$e.jld2" 
    @load filename savemat_xz
    for k in 1:steps
        metts_matrix_xz[e, k] = savemat_xz[k]
    end
end

Ez_arr = zeros(steps)
Exz_arr = zeros(steps)

for i in 1:steps
    for e in 1:ensemble_size
        Ez_arr[i] += mpo_expectation(heisenbergmpo(L,1.0), metts_matrix_z[e, i])
        Exz_arr[i] += mpo_expectation(heisenbergmpo(L,1.0), metts_matrix_xz[e, i])
    end
end

zenergy = Ez_arr/(L*ensemble_size)
xzenergy = Exz_arr/(L*ensemble_size)

plot(1:steps, [zenergy, xzenergy], title="energy per site e vs Steps", label=["Z and X" "Z"], linewidth=3)
xlabel!("Step Number")
ylabel!("Energy per Site")
print(xzenergy-zenergy)#check if they are different

#NEXT -- SUSCEPTIBILITY
sz_arr = zeros(ensemble_size)
sz2_arr = zeros(ensemble_size)
mpo_sz = sz_tot_mpo(L)
mpo_sz2 = sz_squared_mpo(L)

for e in 1:ensemble_size
    sz_arr[e] = mpo_expectation(mpo_sz,metts_matrix_xz[e,steps])
    sz2_arr[e] = mpo_expectation(mpo_sz2,metts_matrix_xz[e,steps])
    println(sz2_arr[e], " is sz2 expectation for ensemble ", e)
end

sz_avg = sum(sz_arr)/ensemble_size
sz2_avg = sum(sz2_arr)/ensemble_size

susceptibility = beta*(sz2_avg - sz_avg^2)/L