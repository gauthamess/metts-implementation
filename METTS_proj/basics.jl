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

    j = m + 2 # so that we can use m_ind as index of mvecs
    mvec = mvecs[:, j] # the state ket(m)

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

function normalise(mps)
    mps[1] = mps[1] / sqrt(normis(mps))
    return mps
end

function normis(mps)
    L = length(mps)
    Id = proj_mpo(reshape(LinearAlgebra.I(3), 1, 3, 1, 3), L, 1)
    normis = mpo_expectation(Id, mps)
    return normis  
end

function mpo_on_mps(mpo, mps)
    L = size(mps, 1)
    mps = sitecanonical(mps, 1)
    for i in 1:L
        mps[i] = applyHtoC(mpo, mps, i)
        if i < L
            mps = sitecanonical(mps, i+1)
        end
    end   
    return mps
end


function computeLeftEnvironment(MPO, MPS, ell_max)
    Lenv = ones((eltype(MPO[1])), (1, 1, 1))  # identity 3-leg tensor

    for ell in 1:ell_max
        Lenv = updateLeft(Lenv, MPS[ell], MPO[ell], MPS[ell])
    end

    return Lenv
end

function computeRightEnvironment(MPO, MPS, ell_min)
    L = length(MPO)
    Renv = ones((eltype(MPO[1])), (1, 1, 1))  # identity 3-leg tensor

    for ell in reverse(ell_min:L)
        Renv = updateLeft(Renv, permutedims(MPS[ell], (3,2,1)), permutedims(MPO[ell], (3,2,1,4)), permutedims(MPS[ell], (3,2,1)))
    end

    return Renv
end

function applyHtoC(W, MPS, ell)
    Lenv = computeLeftEnvironment(W, MPS, ell-1)
    Renv = computeRightEnvironment(W, MPS, ell+1)
    Cell = MPS[ell]
    Well = W[ell]

    HC = contract(Lenv,[3], Cell, [1])
    HC = contract(HC, [2,3], Well,[1,4])
    HC = contract(HC,[2,4], Renv,[3,2])



    return HC/LinearAlgebra.norm(HC)
end

function truncateRight(mps,Nkeep)
    lmps = size(mps,1)
    for i in lmps-1:-1:1
        mtemp = contract(mps[i],3,mps[i+1],1)
        U,S,Vd,_ = svd(mtemp,[1,2],Nkeep=Nkeep,tolerance = tol)
        mps[i] = contract(U,3,diagm(S),1)
        mps[i+1] = Vd
    end
    return mps
end

function mpo_expectation(W::Vector{<:AbstractArray{<:Number,4}}, 
                         MPS::Vector{<:AbstractArray{<:Number,3}})
    # Initialize environment as scalar identity in a 3-leg tensor form
    C = ones(eltype(W[1]), (1, 1, 1))

    for i in eachindex(W)
        C = updateLeft(C, MPS[i], W[i], MPS[i])
    end

    # At the end, C is (1,1,1), so extract scalar value
    return real(C[1,1,1])
end
