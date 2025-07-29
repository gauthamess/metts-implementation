include("basics.jl")
include("states.jl")
include("gates.jl")
tol = 1e-10
using Plots

function cpscollapse(metts, base)
    L = size(metts, 1)

    for i in 1:L
      flag = 0
        # finding probabilities pm
        pm = zeros(3)
        for m = 1:-1:-1
            P = proj(m, base)
            mpo = proj_mpo(P, L, i)
            pm[m+2] = mpo_expectation(mpo, metts) #probability for each collapse
        end
        pm = pm / sum(pm) #normalise just in case

        sumprob =  cumsum(pm)

        proll = rand() #generates random in (0,1)

        for j in 1:length(pm)
          if proll < sumprob[j]
            P = proj(j-2, base)
            mpo = proj_mpo(P, L, i)
            metts[i] = applyHtoC(mpo, metts, i)
            break
          end
        end
    end
    return metts
end


function cpscollapse2(metts, base) #BOND COLLAPSE, NO P(i)
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
            #metts = normalise(sitecanonical(metts,i)) #sitecanonical at site i

            p_i = pm[j]

            bond_right = metts[i][1,j,:]
            if i < L
                size_iplus1 = (1, size(metts[i+1])[2], size(metts[i+1])[3]) # new dimensions of the next site
                metts[i+1] = reshape(contract(bond_right, 1, metts[i+1], 1), size_iplus1)/sqrt(p_i) # contract bond with next site
            end

            c_vec = eigvecs(Sxyz(1,base))[:,j] #get |m> of Sz
            metts[i] = reshape(c_vec,1,3,1) #reshape to rank-3

            break
          end
        end
    end
    return metts
end


