function trotter(spin, beta, Nsteps)
    tauT = beta / Nsteps # Trotter time step

    # Preparing the time evolution hamiltonian
    # Spin-1 local space tensor (rank-3: left-phys-right)
    S = spinlocalspace(rationalize(spin))
    Sx = (S[1] + S[2])/2 # site tensor for Sx
    Sy = (S[1] - S[2]) / (2im) # site tensor for Sy
    Sz = S[3] # site tensor for Sz
    size_s = size(Sx)[1]

    # Construct two-site Hamiltonians via tensor contraction: <S S>
    Ho = LinearAlgebra.kron(Sx, Sx) + LinearAlgebra.kron(Sy, Sy) + LinearAlgebra.kron(Sz, Sz)  # odd bond
    He = Ho # even bond

    # Matrix exponentiation
    expHo_mat = exp(-tauT/2 * Ho)
    expHe_mat = exp(-tauT/2 * He)

    # Reshape back to rank-4 tensors: (left1, left2, right1, right2)
    expHo = reshape(expHo_mat, size_s, size_s, size_s, size_s)

    expHe = reshape(expHe_mat, size_s, size_s, size_s, size_s)

    return expHo, expHe
end

function applygate(mps, l, Hgate, Nkeep)
    # contract the mps sites to produce a rank four tensor
    m = contract(mps[l], 3, mps[l+1], 1)

    # contracting the physical legs with the trotter gate
    mtau = contract(m, [2, 3], Hgate, [3,4])
    mtau = permutedims(mtau,(1,3,4,2))

    # svd to separate back into two sites
    U,S,Vd,_ = svd(mtau,[1, 2], Nkeep = Nkeep)

    # since the site index should move, contract S with Vd
    mps[l] = U
    mps[l+1] = contract(diagm(S), 2, Vd, 1)

    return mps
end

function applygate_backward(mps, l, Hgate, Nkeep)
    # contract the mps sites to produce a rank four tensor
    m = contract(mps[l], 3, mps[l+1], 1)

    # contracting the physical legs with the trotter gate
    mtau = contract(m, [2, 3], Hgate, [3,4])
    mtau = permutedims(mtau,(1,3,4,2))

    # svd to separate back into two sites
    U,S,Vd,_ = svd(mtau,[1, 2], Nkeep = Nkeep)

    # since the site index should move, contract S with Vd
    mps[l] = contract(U, 3, diagm(S), 1)
    mps[l+1] = Vd

    return mps
end

function tdmrg(mps, beta, Nsteps, Nkeep)
mps = rightcanonical(mps) # ensure mps is in right canonical form

    lmps = size(mps, 1)

    # using time beta/2 for sites 1 to L-2 to get tau/2
    expHo, expHe = trotter(1, beta, Nsteps)
    expHoL, expHeL = trotter(1, beta, Nsteps)

    for n in 1:Nsteps
        mps = normalise(mps)
        #println("TDMRG step $n / $Nsteps")
        # evolving sites using bond by bond compression
        for i in 1:(lmps-1)
            if isodd(i)
                mps = applygate(mps, i, expHo, Nkeep)
            else
                mps = applygate(mps, i, expHe, Nkeep)
            end
        
        end

        #sweeping backwards
        for i in (lmps-1):-1:1
            if isodd(i)
                mps = applygate_backward(mps, i, expHo, Nkeep)
            else
                mps = applygate_backward(mps, i, expHe, Nkeep)
            end
        end
    end

    return mps
end

