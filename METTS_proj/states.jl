
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

function ctm(mps, beta, Nkeep) #2000
    # time evolving the CPS by time beta/2
    metts = tdmrg(mps, beta/2, 2000, Nkeep)
    return metts
end
