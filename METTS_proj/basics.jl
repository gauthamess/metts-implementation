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
    L = length(mps)
    Id = proj_mpo(reshape(LinearAlgebra.I(3), 1, 3, 1, 3), L, 1)
    norm = mpo_expectation(Id, mps)
    mps[1] = mps[1] / sqrt(norm)
    return mps
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

    return HC
end

function canonForm(M,id,Nkeep=Inf)

dw = zeros(length(M)-1,1); # discarded weights

# # Bring the left part of MPS into the left-canonical form
for it = (1:id)
    
    # reshape M[it] & SVD
    T = M[it]
    T = reshape(T,(size(T,1)*size(T,2),size(T,3)))
    svdT = LinearAlgebra.svd(T)
    U = svdT.U
    S = svdT.S
    V = svdT.Vt'     
    Svec = S; # vector of singular values
    
    # truncate singular values/vectors; keep up to Nkeep. Truncation at the
    # bond between M[id] & M[id+1] is performed later.
    if ~isinf(Nkeep) && (it < id)
        nk = min(length(Svec),Nkeep); # actual number of singular values/vectors to keep
        dw[it] = dw[it] + sum(Svec[nk+1:end].^2); # discarded weights
        U = U[:,(1:nk)]
        V = V[:,(1:nk)]
        Svec = Svec[1:nk]
    end
    
    S = diagm(Svec); # return to square matrix
    
    # reshape U into rank-3 tensor, & replace M[it] with it
    M[it] = reshape(U,(size(U,1)÷size(M[it],2),size(M[it],2),size(U,2)))
    
    if it < id
        # contract S & V' with M[it+1]
        M[it+1] = contract(S*V',2,2,M[it+1],3,1)
    else
        # R1: tensor which is the leftover after transforming the left
        #   part. It will be contracted with the counterpart R2 which is
        #   the leftover after transforming the right part. Then R1*R2 will
        #   be SVD-ed & its left/right singular vectors will be
        #   contracted with the neighbouring M-tensors.
        R1 = S*V'
    end
    
end

# # In case of fully right-canonical form; the above for-loop is not executed
if id == 0
    R1 = 1
end
    
# # Bring the right part into the right-canonical form
for it = (length(M):-1:id+1)
    
    # reshape M[it] & SVD
    T = M[it]
    T = reshape(T,(size(T,1),size(T,2)*size(T,3)))
    svdT = LinearAlgebra.svd(T)
    U = svdT.U
    S = svdT.S
    V = svdT.Vt' 
    Svec = S; # vector of singular values
    
    # truncate singular values/vectors; keep up to Nkeep. Truncation at the
    # bond between M[id] & M[id+1] is performed later.
    if ~isinf(Nkeep) && (it > (id+1))
        nk = min(length(Svec),Nkeep); # actual number of singular values/vectors to keep
        dw[it-1] = dw[it-1] + sum(Svec[nk+1:end].^2); # discarded weights
        U = U[:,(1:nk)]
        V = V[:,(1:nk)]
        Svec = Svec[1:nk]
    end
    
    S = diagm(Svec); # return to square matrix
    
    # reshape V' into rank-3 tensor, replace M[it] with it
    M[it] = reshape(V',(size(V,2),size(M[it],2),size(V,1)÷size(M[it],2)))
    
    if it > (id+1)
        # contract U & S with M[it-1]
        M[it-1] = contract(M[it-1],3,3,U*S,2,1)
    else
        # R2: tensor which is the leftover after transforming the right
        #   part. See the description of R1 above.
        R2 = U*S
    end
    
end

# # In case of fully left-canonical form; the above for-loop is not executed
if id == length(M)
    R2 = 1
end

# # SVD of R1*R2; & contract the left/right singular vectors to the tensors

    svdT = LinearAlgebra.svd(R1*R2)
    U = svdT.U
    S = svdT.S
    V = svdT.Vt' 

# truncate singular values/vectors; keep up to Nkeep. At the leftmost &
# rightmost legs [dummy legs], there should be no truncation, since they
# are already of size 1.
if ~isinf(Nkeep) && (id > 0) && (id < length(M))
    
    nk = min(length(S),Nkeep); # actual number of singular values/vectors
    dw[id] = dw[id] + sum(S[nk+1:end].^2); # discarded weights
    U = U[:,(1:nk)]
    V = V[:,(1:nk)]
    S = S[1:nk]
    
end

if id == 0 # fully right-canonical form
    # U is a single number which serves as the overall phase factor to the
    # total many-site state. So we can pass over U to V'.
    M[1] = contract(U*V',2,2,M[1],3,1)
elseif id == length(M) # fully left-canonical form
    # V' is a single number which serves as the overall phase factor to the
    # total many-site state. So we can pass over V' to U.
    M[end] = contract(M[end],3,3,U*V',2,1)
else
    M[id] = contract(M[id],3,3,U,2,1)
    M[id+1] = contract(V',2,2,M[id+1],3,1)
end

return M
end