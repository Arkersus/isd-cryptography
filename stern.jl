function decodeSTERN(H::Matrix{Bool}, s::Vector{Bool}, t::Int, p::Int, l::Int, niter::Int)
    r, n = size(H)

    It = 0
    while(It < niter || niter == -1)
        perm = randperm(n)
        H1 = H[:, perm[1:p]]
        H2 = H[:, perm[p+1:end]]

        L1 = [falses(p) for _ in 1:(1 << l)]
        L2 = [falses(n - p) for _ in 1:(1 << l)]

        for i in 1:(1 << l)
            for j in 1:l
                if i & (1 << (j - 1)) != 0
                    L1[i][j] = true
                else
                    L2[i][j] = true
                end
            end
        end

        for x in L1
            s1 = multiply(H1, Vector{Bool}(x))
            for y in L2
                s2 = multiply(H2, Vector{Bool}(y))
                diff = falses(r)
                for i in 1:r
                    diff[i] = s[i] ⊻ s1[i] ⊻ s2[i]
                end

                if sum(diff) == 0
                    e_out = falses(n)
                    for i in 1:p
                        e_out[perm[i]] = x[i]
                    end
                    for i in 1:(n - p)
                        e_out[perm[p + i]] = y[i]
                    end

                    if sum(e_out) != t
                        continue
                    end

                    s_new = multiply(H, Vector{Bool}(e_out))

                    if s_new != s
                        continue
                    end

                    return e_out
                end
            end
        end
    end

    return nothing
end
