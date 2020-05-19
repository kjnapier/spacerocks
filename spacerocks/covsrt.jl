function covsrt(covar, ma, ia, mfit)

    for i in mfit+1:ma
        for j in 1:i
            covar[i][j] = covar[j][i] = 0.0
        end
    end

    k = mfit

    for j in ma:-1:1
        if ia[j] != 0
            for i in 1:ma
                covar[i][k], covar[i][j] = covar[i][j], covar[i][k]
            end
            for i in 1:ma
                covar[k][i], covar[j][i] = covar[j][i], covar[k][i]
            end
            k -= 1
        end
    end

    return covar

end
