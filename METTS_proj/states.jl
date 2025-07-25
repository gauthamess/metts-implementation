
function cps_z(spin, N; base = 3) #Sz product states # Eigenvectors in columns
    
    d = size(Sxyz(1,base), 1)  # Physical dimension (should be 2 for spin-1/2)
    
    mps = Vector{Array{ComplexF64,3}}(undef, N)
    
    for i in 1:N
        j = rand(1:d)  # Random eigenvector 
        v = eigvecs(Sxyz(1,base))[:,j]  # d-dimensional vector
        
        
        mps[i] = reshape(v, 1, d, 1) # Convert to MPS tensor: (1, d, l1)
    end
    
    return mps
end

function cpsrandom(spin, N) #Sz product states # Eigenvectors in columns
    
    d = size(Sxyz(1,1), 1)  # Physical dimension (should be 2 for spin-1/2)
    
    mps = Vector{Array{ComplexF64,3}}(undef, N)

    roll = rand()
    for i in 1:N
        j = rand(1:d)
        if roll > 0.5
            v = eigvecs(Sxyz(1,3))[:,j] 
            mps[i] = reshape(v,1,d,1)
        else
            v = eigvecs(Sxyz(1,1))[:,j] 
            mps[i] = reshape(v,1,d,1)
        end
    end
    return mps
end

function ctm(psi_i, beta, Nkeep)
    L = size(psi_i, 1)

    # time evolving the CPS by time beta/2
    mnew = tdmrg(psi_i, beta, 2000, Nkeep)

    # calculating metts = P_i^(-1/2) * e^(- beta H / 2) psi_i
    #mnew = normalise(mnew)
    return mnew
end