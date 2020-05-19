function gaussj(a, n, b, m)

    ipiv = zeros(N)
    indxc = zeros(N)
    indxr = zeros(N)

    for i in 1:n
        big = 0.0
        for j in 1:n
            if ipiv[j] != 1
                for k in 1:n
                    if ipiv[k] == 0
                        if abs(a[j][k]) >= big
                            big = abs(a[j][k])
                            irow = j
                            icol = k
                        end
                    else if ipiv[k] > 1
                        error("gaussj: Singular Matrix-1")
                    end
                end
                ipiv[icol] += 1
                if irow != icol
                    for l in 1:n
                        a[irow][l], a[icol][l] = a[icol][l], a[irow][l]
                    end
                    for l in 1:m
                        b[irow][l], b[icol][l] = b[icol][l], b[irow][l]
                    end
                end
                indxr[i]=irow
		        indxc[i]=icol
                if a[icol][icol] == 0.0
                    error("gaussj: Singular Matrix-2")
                else
                    pivinv=1.0/a[icol][icol]
                end
                a[icol][icol] = 1.0
		        for l in 1:n
                    a[icol][l] *= pivinv
                end
                for l in 1:m
                    b[icol][l] *= pivinv
                end
                for ll in 1:n
                    if ll != icol
                        dum = a[ll][icol]
                        a[ll][icol]=0.0
                        for l in 1:n
                            a[ll][l] -= a[icol][l] * dum
                        end
                        for l in 1:n
                            b[ll][l] -= b[icol][l] * dum
                        end
                    end
                end
            end
            for l in n:-1:1
                if indxr[l] != indxc[l]
			       for k in 1:n
                       a[k][indxr[l]],a[k][indxc[l]] = a[k][indxc[l]], a[k][indxr[l]]
                   end
               end
           end

       end

   end

    return a, b

end
