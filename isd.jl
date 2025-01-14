using LinearAlgebra
using Random

function decodeISD(G::Array{Int, 2}, y::Array{Int, 1}, t::Int, niter::Int=-1)

    k, n = size(G)
    if n != length(y)
        throw(ArgumentError("Cannot multiply matrices"))
    end

    It = 0
    while(It < niter || niter == -1)
        It += 1
        I_ = sort(randperm(n)[1:k])
        G_I = G[:, I_]
        if rank(G_I) < k
            continue
        end

        m = mod.(round.(G_I \ y[I_]), 2)
        e = mod.(y - mod.(G' * m,2), 2)

        if sum(Int.(e)) == t && mod.(mod.(G' * m,2) + e, 2) == y
            return (Int.(m), Int.(e))
        end
    end
    return niter, -1
end

m = [0,1,0,1]
G = [1 0 0 0 0 0 0; 0 1 0 0 1 1 0; 0 0 1 0 1 0 1; 0 0 0 1 1 1 1]
y = [0, 1, 0, 0, 0, 0, 1]
t = 1

println("decodeISD")
m_new,e_new= @time decodeISD(G,y,t,1000)
if e_new!=-1
        println("Original m:",m)
	println("ISD decoded!");
	println("decoded m",m_new) 
	println("decoded e",e_new) 
	noisy_codeword_new=Int.(mod.(mod.(G'*m_new,2)+e_new,2))
else
	println("ISD failed to decode the message!")
end

