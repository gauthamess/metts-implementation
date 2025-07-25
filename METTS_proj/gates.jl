function trotter(spin, beta, Nsteps)
    tauT = beta / Nsteps # Trotter time step

    # Preparing the time evolution hamiltonian
    S = rationalize(spin)
    Sx = Sxyz(S,1)
    Sy = Sxyz(S,2)
    Sz = Sxyz(S,3)
    d = size(Sx)[1]

    # Construct two-site Hamiltonians via tensor contraction: <S S>
    Hloc = LinearAlgebra.kron(Sx, Sx) + LinearAlgebra.kron(Sy, Sy) + LinearAlgebra.kron(Sz, Sz)  # odd bond

    # Matrix exponentiation
    expH_mat = exp(-tauT * Hloc)

    # Reshape back to rank-4 tensors: (left1, left2, right1, right2)
    expH = reshape(expH_mat, d,d,d,d)

    return expH
end

function applygate(mps, l, Hgate)
    # contract the mps sites to produce a rank four tensor
    m = contract(mps[l], 3, mps[l+1], 1)

    # contracting the physical legs with the trotter gate
    mtau = contract(m, [2, 3], Hgate, [3,4])
    mtau = permutedims(mtau,(1,3,4,2))

    # svd to separate back into two sites
    U,S,Vd,_ = svd(mtau,[1, 2])

    # since the site index should move, contract S with Vd
    mps[l] = U
    mps[l+1] = contract(diagm(S), 2, Vd, 1)

    return mps
end

function applygate_backward(mps, l, Hgate)
    # contract the mps sites to produce a rank four tensor
    m = contract(mps[l], 3, mps[l+1], 1)

    # contracting the physical legs with the trotter gate
    mtau = contract(m, [2, 3], Hgate, [3,4])
    mtau = permutedims(mtau,(1,3,4,2))

    # svd to separate back into two sites
    U,S,Vd,_ = svd(mtau,[1, 2])

    # since the site index should move, contract S with Vd
    mps[l] = contract(U, 3, diagm(S), 1)
    mps[l+1] = Vd

    return mps
end

function tdmrg(mps, beta, Nsteps, Nkeep)
    lmps = size(mps, 1)
    # using time beta/2 for sites 1 to L-2 to get tau/2
    expH = trotter(1, beta/2, Nsteps)

    for n in 1:Nsteps
        #1. FORWARD SWEEP
        for i in 1:2:floor(Int, (2*lmps-1)/2) #ODD SWEEP
            mps = applygate(mps, i, expH)
        end
        for i in 2:2:floor(Int, (2*lmps-1)/2) #EVEN SWEEP
            mps = applygate(mps, i, expH)
        end
        mps = leftcanonical(mps,Nkeep)
        #2. BACKWARD SWEEP
        for i in (2*floor(Int, (lmps-2)/2) + 1):-2:1 #ODD SWEEP
            mps = applygate(mps, i, expH)
        end
        for i in (2*floor(Int, (lmps-1)/2)):-2:2 #EVEN SWEEP
            mps = applygate(mps, i, expH)
        end
    end
    mps = leftcanonical(mps,Nkeep=Nkeep)
    mps = normalise(mps)
    return mps
end