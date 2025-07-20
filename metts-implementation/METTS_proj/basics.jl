import tn_julia: svd, contract, updateLeft, identity, permutedims, spinlocalspace, reshape, rightcanonical, leftcanonical, mpo_expectation, applyHtoC, computeLeftEnvironment, computeRightEnvironment, sitecanonical
using LinearAlgebra

function Sxyz(spin,index)

    S = spinlocalspace(rationalize(spin))
    Sxyz = []

    Sx = (S[1] + S[2])/2 # site tensor for Sx
    Sy = (S[1] - S[2]) / (2im) # site tensor for Sy

    push!(Sxyz, Sx)
    push!(Sxyz, Sy)
    push!(Sxyz, S[3])

    return Sxyz[index]
end

function proj(m, base) # m are the eigvenvalues, takes values -1, 0, 1, returns |m><m|, base is 1:x, 2:y, 3:z

    # getting set of states ket(m)
    mvecs = eigvecs(Sxyz(1,base))

    m_ind = m + 2 # so that we can use m_ind as index of mvecs
    mvec = mvecs[:, m_ind] # the state ket(m)

    # creating projector P_i = ket(m) bra(m)
    proj = LinearAlgebra.kron(conj(mvec)', mvec)

    # reshaping to produce local operator to be inserted in MPO
    proj = reshape(proj, 1, 3, 1, 3)

    return proj
end

function proj_mpo(P, L, i)
    d = size(P, 2)
    mpo = Array{Any}(undef, L)
    for n in 1:i-1
        mpo[n] = reshape(LinearAlgebra.I(d), 1, 3, 1, 3)
    end
    mpo[i] = P
    for n in i+1:L
        mpo[n] = reshape(LinearAlgebra.I(d), 1, 3, 1, 3)
    end

    return Vector{Array{ComplexF64,4}}(mpo) #force compatibility with expectation function
end

function mps_inner_product(A, B)
    # Start with 1x1 identity matrix
    env = ones(ComplexF64, 1, 1)
    
    for i in 1:length(A)
        # Contract environment with A site
        temp = contract(env, 2, A[i], 1)  # env × A[i]
        
        # Contract with conjugate of B site
        env = contract(temp, [1, 2], conj(B[i]), [1, 2])
    end
    
    return env[1, 1]  # Scalar value
end

function mps_norm(mps)
    return sqrt(real(mps_inner_product(mps, mps)))
end


function heisenbergmpo(L::Int, J::Float64 = 1.0)
    # Spin-1 local space operators
    Splus, Sminus, Sz, Id = spinlocalspace(rationalize(1))
    Sx = (Splus + Sminus) / sqrt(2)
    Sy = (Splus - Sminus) / (1im * sqrt(2))

    # MPO bond dimension: 5, physical dim: 3
    W = zeros(ComplexF64, 5, 3, 5, 3)  # [left, phys_in, right, phys_out]

    # Bulk MPO rules
    W[1, :, 1, :] = Id

    W[2, :, 1, :] = Sx
    W[3, :, 1, :] = Sy
    W[4, :, 1, :] = Sz

    W[5, :, 2, :] = J * Sx
    W[5, :, 3, :] = J * Sy
    W[5, :, 4, :] = J * Sz
    W[5, :, 5, :] = Id

    # Edges (correct shapes)
    W_first = W[[5], :, :, :]   # shape (1, 3, 5, 3)
    W_last = W[:, :, [1], :]    # shape (5, 3, 1, 3)

    # Build MPO
    mpo = [W_first]
    for _ in 2:(L - 1)
        push!(mpo, copy(W))
    end
    push!(mpo, W_last)

    return mpo
end