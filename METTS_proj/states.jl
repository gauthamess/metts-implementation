
function cps_z(spin, N) #Sz product states # Eigenvectors in columns
    
    d = size(Sxyz(1,3), 1)  # Physical dimension (should be 2 for spin-1/2)
    
    mps = Vector{Array{ComplexF64,3}}(undef, N)
    
    for i in 1:N
        j = rand(1:d)  # Random eigenvector 
        v = eigvecs(Sxyz(1,3))[:,j]  # d-dimensional vector
        
        
        mps[i] = reshape(v, 1, d, 1) # Convert to MPS tensor: (1, d, l1)
    end
    
    return mps
end

function ctm(psi_i, beta, Nkeep)
    L = size(psi_i, 1)

    # time evolving the CPS by time beta/2
    mnew = tdmrg(psi_i, beta/2, 2000, Nkeep)

    # calculating P(i)
    Id = proj_mpo(reshape(LinearAlgebra.I(3), 1, 3, 1, 3), L, 1)

    P = mpo_expectation(Id, mnew)

    # calculating metts = P_i^(-1/2) * e^(- beta H / 2) psi_i
    mnew[1] = mnew[1]/sqrt(P)
    return mnew
end

Id = proj_mpo(reshape(LinearAlgebra.I(3), 1, 3, 1, 3), 10, 1)

P = mpo_expectation(Id, cps_z(1,10))