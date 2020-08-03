function SpareseJSR0!(A,d;lb=0,ub=2,tol=1e-5,QUIET=false)
    n=size(A[1],2)
    m=length(A)
    # zvar=get_zcol(A,n)
    @polyvar x[1:n]
    init_poly=sum(x.^(2d))
    tsupp=Array(Diagonal([2d for i=1:n]))
    for i=1:m
        cons_poly=init_poly(x=>A[i]*x)
        mon=monomials(cons_poly)
        lm=length(mon)
        cons_supp=zeros(UInt8,n,lm)
        for j=1:lm, k=1:n
            cons_supp[k,j]=MultivariatePolynomials.degree(mon[j],x[k])
        end
        tsupp=[tsupp cons_supp]
    end
    tsupp=sortslices(tsupp,dims=2,rev=true)
    tsupp=unique(tsupp,dims=2)
    ltsupp=size(tsupp,2)
    ind=zeros(UInt16, n)
    for i=1:n
        temp=zeros(UInt8, n)
        temp[i]=2d
        ind[i]=bfind(tsupp,ltsupp,temp)
    end
    basis=get_hbasis(n,d)
    basis=generate_basis!(n,tsupp,basis)
    blocks,cl,blocksize=get_blocks(n,tsupp,basis,QUIET=QUIET)
    psupp=[prod(x.^tsupp[:,i]) for i=1:ltsupp]
    while ub-lb>tol
        gamma=(lb+ub)/2
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), true)
        coe0=@variable(model, [1:ltsupp])
        p=coe0'*psupp
        @constraint(model, coe0[ind].>=0.1)
        for k=1:m+1
            cons=[AffExpr(0) for i=1:ltsupp]
            for i=1:cl
                if blocksize[i]==1
                   pos=@variable(model, lower_bound=0)
                   bi=2*basis[:,blocks[i]]
                   Locb=bfind(tsupp,ltsupp,bi)
                   cons[Locb]+=pos
                else
                   pos=@variable(model, [1:blocksize[i], 1:blocksize[i]], PSD)
                   for j=1:blocksize[i], r=j:blocksize[i]
                       bi=basis[:,blocks[i][j]]+basis[:,blocks[i][r]]
                       Locb=bfind(tsupp,ltsupp,bi)
                       if j==r
                          cons[Locb]+=pos[j,r]
                       else
                          cons[Locb]+=2*pos[j,r]
                       end
                   end
                end
            end
            if k>1
                coe=coefficients(gamma^(2d)*p-p(x=>A[k-1]*x))
            else
                coe=coe0
            end
            @constraint(model, cons.==coe)
        end
        # if zvar!=[]
        #     temp=zeros(UInt8,n)
        #     temp[zvar[1]]=2d
        #     Locb=bfind(tsupp,ltsupp,temp)
        #     @constraint(model, coe0[Locb]==0.1)
        # else
            # @constraint(model, coe0[1]==0.1)
        # end
        optimize!(model)
        if termination_status(model)==MOI.OPTIMAL
            ub=gamma
        else
            lb=gamma
        end
    end
    return ub
end

function SpareseJSR!(A,d;lb=0,ub=2,tol=1e-5,TS="block",QUIET=false)
    n=size(A[1],2)
    m=length(A)
    @polyvar x[1:n]
    init_poly=sum(x.^(2d))
    supp=Vector{Array{UInt8, 2}}(undef, m+1)
    tsupp=Vector{Array{UInt8, 2}}(undef, m+1)
    gbasis=Vector{Array{UInt8, 2}}(undef, m+1)
    ltsupp=Vector{UInt16}(undef, m+1)
    supp[1]=Array(Diagonal([2d for i=1:n]))
    for i=1:m
        cons_poly=init_poly(x=>A[i]*x)
        mon=monomials(cons_poly)
        lm=length(mon)
        cons_supp=zeros(UInt8,n,lm)
        for j=1:lm, k=1:n
            cons_supp[k,j]=MultivariatePolynomials.degree(mon[j],x[k])
        end
        supp[1]=[supp[1] cons_supp]
    end
    supp[1]=sortslices(supp[1],dims=2,rev=true)
    supp[1]=unique(supp[1],dims=2)
    ind=zeros(UInt16, n)
    for i=1:n
        temp=zeros(UInt8, n)
        temp[i]=2d
        ind[i]=bfind(supp[1],size(supp[1],2),temp)
    end
    basis=get_hbasis(n,d)
    blocks=Vector{Vector{Vector{Int64}}}(undef, m+1)
    cl=Vector{Int64}(undef, m+1)
    blocksize=Vector{Vector{Int64}}(undef, m+1)
    sLocb=Vector{Vector{UInt16}}(undef, m+1)
    pcoe=ones(Float64, size(supp[1],2))
    psupp=[prod(x.^supp[1][:,i]) for i=1:size(supp[1],2)]
    p=pcoe'*psupp
    for i=1:m+1
        if i>1
            mon=monomials(p(x=>A[i-1]*x))
            lm=length(mon)
            supp[i]=zeros(UInt8,n,lm)
            for j=1:lm, k=1:n
                supp[i][k,j]=MultivariatePolynomials.degree(mon[j],x[k])
            end
            supp[i]=[supp[1] supp[i]]
            supp[i]=sortslices(supp[i],dims=2,rev=true)
            supp[i]=unique(supp[i],dims=2)
        end
        gbasis[i]=generate_basis!(n,supp[i],basis)
        blocks[i],cl[i],blocksize[i]=get_blocks(n,supp[i],gbasis[i],TS=TS,QUIET=QUIET)
        tsupp[i]=zeros(UInt8,n,Int(sum(blocksize[i].^2+blocksize[i])/2))
        s=1
        for k=1:cl[i], j=1:blocksize[i][k], r=j:blocksize[i][k]
            bi=gbasis[i][:,blocks[i][k][j]]+gbasis[i][:,blocks[i][k][r]]
            tsupp[i][:,s]=bi
            s+=1
        end
        tsupp[i]=sortslices(tsupp[i],dims=2,rev=true)
        tsupp[i]=unique(tsupp[i],dims=2)
        ltsupp[i]=size(tsupp[i],2)
        sLocb[i]=zeros(UInt16, size(supp[i],2))
        for k=1:size(supp[i],2)
            sLocb[i][k]=bfind(tsupp[i],ltsupp[i],supp[i][:,k])
        end
    end
    while ub-lb>tol
        gamma=(lb+ub)/2
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), true)
        pcoe=@variable(model, [1:size(supp[1],2)])
        p=pcoe'*psupp
        @constraint(model, pcoe[ind].>=0.1)
        for k=1:m+1
            cons=[AffExpr(0) for i=1:ltsupp[k]]
            for i=1:cl[k]
                if blocksize[k][i]==1
                   pos=@variable(model, lower_bound=0)
                   bi=2*gbasis[k][:,blocks[k][i]]
                   Locb=bfind(tsupp[k],ltsupp[k],bi)
                   cons[Locb]+=pos
                else
                   pos=@variable(model, [1:blocksize[k][i], 1:blocksize[k][i]], PSD)
                   for j=1:blocksize[k][i], r=j:blocksize[k][i]
                       bi=gbasis[k][:,blocks[k][i][j]]+gbasis[k][:,blocks[k][i][r]]
                       Locb=bfind(tsupp[k],ltsupp[k],bi)
                       if j==r
                          cons[Locb]+=pos[j,r]
                       else
                          cons[Locb]+=2*pos[j,r]
                       end
                   end
                end
            end
            if k>1
                coe=coefficients(gamma^(2d)*p-p(x=>A[k-1]*x))
            else
                coe=pcoe
            end
            bc=[AffExpr(0) for i=1:ltsupp[k]]
            bc[sLocb[k]].=coe
            @constraint(model, cons.==bc)
        end
        optimize!(model)
        if termination_status(model)==MOI.OPTIMAL
            ub=gamma
        else
            lb=gamma
        end
    end
    return ub
end

function JSR!(A,d;lb=0,ub=2,tol=1e-5)
    n=size(A[1],2)
    m=length(A)
    @polyvar x[1:n]
    basis=get_hbasis(n,d)
    lbasis=size(basis,2)
    supp=get_hbasis(n,2d)
    lp=size(supp,2)
    ind=zeros(UInt16, n)
    for i=1:n
        temp=zeros(UInt8, n)
        temp[i]=d
        ind[i]=bfind(basis,lbasis,temp)
    end
    while ub-lb>tol
        gamma=(lb+ub)/2
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), true)
        pos=@variable(model, [1:lbasis, 1:lbasis], PSD)
        xbasis=[prod(x.^basis[:,i]) for i=1:lbasis]
        p=xbasis'*pos*xbasis
        @constraint(model, [i in ind], pos[i, i]>=1)
        # @constraint(model, pos[1, 1]==0.1)
        for k=1:m
            cons=[AffExpr(0) for i=1:lp]
            pos=@variable(model, [1:lbasis, 1:lbasis], PSD)
            for j=1:lbasis, r=j:lbasis
                bi=basis[:,j]+basis[:,r]
                Locb=bfind(supp,lp,bi)
                if j==r
                   cons[Locb]+=pos[j,r]
                else
                   cons[Locb]+=2*pos[j,r]
                end
            end
            coe=coefficients(gamma^(2d)*p-p(x=>A[k]*x))
            @constraint(model, cons.==coe)
        end
        optimize!(model)
        # MOI.write_to_file(backend(model).optimizer.model, "E:\\Programs\\SparseJSR\\model3.cbf")
        if termination_status(model)==MOI.OPTIMAL
            ub=gamma
        else
            lb=gamma
        end
    end
    return ub
end

function get_hbasis(n,d)
    lb=binomial(n-1+d,d)
    basis=zeros(UInt8,n,lb)
    basis[1,1]=d
    t=1
    while t<lb
        t+=1
        j=findfirst(x->basis[x,t-1]!=0, 1:n)
        basis[:,t]=basis[:,t-1]
        if j==1
           basis[1,t]-=1
           basis[2,t]+=1
        else
           basis[1,t]=basis[j,t]-1
           basis[j,t]=0
           basis[j+1,t]+=1
        end
    end
    return basis[end:-1:1,end:-1:1]
end

function comp(a,b)
    i=1
    while i<=length(a)
        if a[i]<b[i]
            return -1
        elseif a[i]>b[i]
            return 1
        else
            i+=1
        end
    end
    return 0
end

function bfind(A,l,a)
    if l==0
        return 0
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if length(a)>1
            order=comp(A[:,mid],a)
        else
            order=comp(A[mid],a)
        end
        if order==0
           return mid
        elseif order>0
           low=mid+1
        else
           high=mid-1
        end
    end
    return 0
end

function generate_basis!(n,supp,basis)
    supp=sortslices(supp,dims=2,rev=true)
    supp=unique(supp,dims=2)
    lsupp=size(supp,2)
    lb=size(basis,2)
    indexb=UInt32[]
    for i = 1:lb, j = i:lb
        bi=basis[:,i]+basis[:,j]
        if bfind(supp,lsupp,bi)!=0
             push!(indexb,i,j)
        end
    end
    sort!(indexb)
    unique!(indexb)
    return basis[:,indexb]
end

function get_blocks(n,supp,basis;TS="block",QUIET=true)
    if TS==false
        blocksize=[size(basis,2)]
        blocks=[[i for i=1:size(basis,2)]]
        cl=1
    else
        lb=size(basis,2)
        G=SimpleGraph(lb)
        tsupp=[supp 2*basis]
        tsupp=sortslices(tsupp,dims=2,rev=true)
        tsupp=unique(tsupp,dims=2)
        ltsupp=size(tsupp,2)
        for i = 1:lb, j = i+1:lb
            bi=basis[:,i]+basis[:,j]
            if bfind(tsupp,ltsupp,bi)!=0
               add_edge!(G,i,j)
            end
        end
        if TS=="block"
            blocks=connected_components(G)
            blocksize=length.(blocks)
            cl=length(blocksize)
        else
            blocks,cl,blocksize=chordal_cliques!(G, method=TS, minimize=false)
        end
    end
    nub=unique(blocksize)
    nsizes=[sum(blocksize.== i) for i in nub]
    if QUIET==false
        println("------------------------------------------------------")
        println("The sizes of blocks:\n$nub\n$nsizes")
        println("------------------------------------------------------")
    end
    return blocks,cl,blocksize
end

# function get_zcol(A,n)
#     m=length(A)
#     zerocol=[UInt8[] for i=1:m]
#     for i=1:m
#         for j=1:n
#             if all(x->x==0, A[i][:,j])
#                 push!(zerocol[i],j)
#             end
#         end
#     end
#     zcol=zerocol[1]
#     for i=2:m
#         intersect!(zcol,zerocol[i])
#     end
#     return zcol
# end

# function JSR!(A,d;lb=0,ub=2,tol=1e-5)
#     n=size(A[1],2)
#     m=length(A)
#     @polyvar x[1:n]
#     basis=get_hbasis(n,d,rev=true)
#     lbasis=size(basis,2)
#     supp=get_hbasis(n,2d,rev=true)
#     lp=size(supp,2)
#     psupp=[prod(x.^supp[:,i]) for i=1:lp]
#     @polyvar coe[1:lp]
#     p=coe'*psupp
#     asub=Vector{Vector{UInt16}}(undef, (m+1)*lp)
#     aval=Vector{Vector{Float32}}(undef, (m+1)*lp)
#     for j=1:lp
#         asub[j]=[j]
#         aval[j]=[-1]
#     end
#     for i=1:m
#         for j=1:lp
#             coeff=coefficient(subs(p, x=>A[i]*x)-gamma^(2d)*p, psupp[j], x)
#             temp=sparse(coefficients(coeff, coe))
#             asub[i*lp+j]=temp.nzind
#             aval[i*lp+j]=temp.nzval
#         end
#     end
#     consi=[UInt16[] for i=1:lp]
#     consj=[UInt16[] for i=1:lp]
#     consk=[Float64[] for i=1:lp]
#     for j=1:lbasis
#         for r=j:lbasis
#             @inbounds bi=basis[:,j]+basis[:,r]
#             Locb=bfind(supp,lp,bi)
#             @inbounds push!(consi[Locb],r)
#             @inbounds push!(consj[Locb],j)
#             @inbounds push!(consk[Locb],1.0)
#         end
#     end
#     maketask() do task
#         printstream(msg)=print(msg)
#         putstreamfunc(task,MSK_STREAM_LOG,printstream)
#         appendvars(task,lp)
#         appendcons(task,lp*(m+1)+1)
#         appendbarvars(task, [lbasis for i=1:m+1])
#         for i=1:lp
#             putvarbound(task, i, MSK_BK_FR, -Inf, +Inf)
#         end
#         for i=1:(m+1)*lp
#             putconbound(task, i, MSK_BK_FX, 0, 0)
#             putarow(task, i, asub[i], aval[i])
#         end
#         putconbound(task, (m+1)*lp+1, MSK_BK_FX, 1, 1)
#         putarow(task, (m+1)*lp+1, [1], [1])
#         for j=1:m+1
#             for k=1:lp
#                 cons=appendsparsesymmat(task,lbasis,consi[k],consj[k],consk[k])
#                 putbaraij(task,(j-1)*lp+k,j,[cons],[1.0])
#             end
#         end
#         optimize(task)
#         solutionsummary(task, MSK_STREAM_MSG)
#         solsta=getsolsta(task, MSK_SOL_ITR)
#         if solsta==MSK_SOL_STA_OPTIMAL
#            println("Optimal!")
#         end
#     end
# end
#
# maketask() do task
#     printstream(msg)=print(msg)
#     putstreamfunc(task, MSK_STREAM_LOG, printstream)
#     readdata(task, "E:\\Programs\\SparseJSR\\model2.cbf")
#     optimize(task)
#     solutionsummary(task, MSK_STREAM_MSG)
#     solsta=getsolsta(task, MSK_SOL_ITR)
#     if solsta==MSK_SOL_STA_OPTIMAL
#        println("Optimal!")
#     end
# end

# function comp(a,b)
#     i=length(a)
#     while i>=1
#           if a[i]<b[i]
#              return -1
#           elseif a[i]>b[i]
#              return 1
#           else
#              i-=1
#           end
#     end
#     return 0
# end

# function SpareseJSR2(A,d;lb=0,ub=2,tol=1e-5)
#     n=size(A[1],2)
#     m=length(A)
#     var1=get_zcol(A,n)
#     lvar1=length(var1)
#     var2=UInt8[i for i=1:n]
#     setdiff!(var2,var1)
#     lvar2=length(var2)
#     @polyvar x[1:n]
#     coe=Vector{Vector{GenericAffExpr{Float64, VariableRef}}}(undef, m+2)
#     basis1=get_hbasis(lvar1,d)
#     basis2=get_hbasis(lvar2,d)
#     supp1=get_hbasis(lvar1,2d)
#     lp1=size(supp1,2)
#     supp2=get_hbasis(lvar2,2d)
#     lp2=size(supp2,2)
#     while ub-lb>tol
#         gamma=(lb+ub)/2
#         model=Model(optimizer_with_attributes(Mosek.Optimizer))
#         set_optimizer_attribute(model, MOI.Silent(), true)
#         xpsupp1=[prod(x[var1].^supp1[:,i]) for i=1:lp1]
#         coe[1]=@variable(model, [1:lp1])
#         xpsupp2=[prod(x[var2].^supp2[:,i]) for i=1:lp2]
#         coe[2]=@variable(model, [1:lp2])
#         p2=coe[2]'*xpsupp2
#         p=coe[1]'*xpsupp1+p2
#         for i=1:m
#             np=gamma^(2d)*p2-p(x=>A[i]*x)
#             coe[i+2]=coefficients(np)
#         end
#         for k=1:m+2
#             if k==1
#                 basis=basis1
#                 lp=lp1
#                 tsupp=supp1
#             else
#                 basis=basis2
#                 lp=lp2
#                 tsupp=supp2
#             end
#             lbasis=size(basis,2)
#             cons=[AffExpr(0) for i=1:lp]
#             pos=@variable(model, [1:lbasis, 1:lbasis], PSD)
#             for j=1:lbasis
#                 for r=j:lbasis
#                     bi=basis[:,j]+basis[:,r]
#                     Locb=bfind(tsupp,lp,bi)
#                     if j==r
#                        cons[Locb]+=pos[j,r]
#                     else
#                        cons[Locb]+=2*pos[j,r]
#                     end
#                 end
#             end
#             @constraint(model, cons.==coe[k])
#         end
#         @constraint(model, coe[1][1]==1)
#         optimize!(model)
#         # if termination_status(model)==MathOptInterface.OPTIMAL
#         if primal_status(model)==MathOptInterface.NO_SOLUTION
#             lb=gamma
#         else
#             ub=gamma
#         end
#     end
#     return gamma
# end
