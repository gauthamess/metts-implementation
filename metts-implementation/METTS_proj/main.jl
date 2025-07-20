include("basics.jl")
include("states.jl")
include("gates.jl")

function cpscollapse(metts, base)
    L = size(metts, 1)

    for i in 1:L
        # finding probabilities pm
        pm = zeros(3)
        for m = -1:1:1
            P = proj(m, base)
            mpo = proj_mpo(P, L, i)
            pm[m+2] = mpo_expectation(mpo, metts) # probability for each collapse
        end
        pm = pm / sum(pm) # normalise probabilities just in case
        sumprob = cumsum(pm)

        proll = rand() # generates random in (0,1)

        for j in 1:length(pm)
          if proll < sumprob[j] #then outcome is j-2
            m = j-2
            metts = sitecanonical(metts,i) #sitecanonical at site i

            bond_right = metts[i][1,j,:]
            if i < L
                size_iplus1 = (1, size(metts[i+1])[2], size(metts[i+1])[3]) # new dimensions of the next site
                metts[i+1] = reshape(contract(bond_right, 1, metts[i+1], 1), size_iplus1) # contract bond with next site
            end

            c_vec = eigvecs(Sxyz(1,3))[:,j] #get |m>
            metts[i] = reshape(c_vec,1,3,1) #reshape to rank-3
            break
          end
        end
    end
    return metts
end