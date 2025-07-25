L = 50
beta = 1
steps = 10
ensemble_size = 15
Nkeep = 30
Nsteps = beta*150

savemat = zeros(ensemble_size,steps)
function createEnsemble()
    for e in 1:ensemble_size
        cps = cpsrandom(1,L) #new seed for each loop
        cps_xz = copy(cps)
        cps_z = copy(cps)
        for k in 1:steps
            metts = tdmrg(cps_xz,beta,Nsteps,Nkeep)
            savemat[e,k] = metts
            writedlm("output_test.txt", savemat)
            if rolls > 0.5
                cps_xz = cpscollapse(metts,1)
            else
                cps_xz = cpscollapse(metts,3)
            end
        end
    end

end
