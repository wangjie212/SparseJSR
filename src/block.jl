function permTriang(A)
    n=size(A[1], 2)
    m=length(A)
    sA=sum(abs.(A[k]) for k=1:m)
    sA-=diagm(diag(sA))
    blocks=strongly_connected_components(DiGraph(sA.>0))
    if length(blocks)>1
        bA=Vector{Vector{Array{Float32, 2}}}(undef, length(blocks))
        for j=1:length(blocks)
            bA[j]=Vector{Array{Float32, 2}}(undef, m)
            for i=1:m
                bA[j][i]=A[i][blocks[j],blocks[j]]
            end
        end
        return bA,1
    else
        return A,0
    end
end

function checkblock(bA, m)
    nb = length(bA)
    Blb = zeros(nb)
    Bub = zeros(nb)
    for i = 1:nb
        if size(bA[i][1], 2)==1
            Blb[i] = maximum([abs(bA[i][j][1,1]) for j=1:m])
            Bub[i] = Blb[i]
        else
            Blb[i],Bub[i]=gripenberg(discreteswitchedsystem(bA[i]), δ=0.2, verbose = 0)
        end
    end
    mBlb=maximum(Blb)
    ind=[Bub[i]>mBlb for i = 1:nb]
    nnb=sum(ind)
    println("Split into $nb blocks. Preserve $nnb blocks.")
    return bA[ind],Blb[ind],Bub[ind],nnb
end
