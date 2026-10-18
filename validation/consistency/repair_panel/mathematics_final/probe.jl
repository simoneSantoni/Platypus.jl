# Independent 512-bit oracle; Base only. No package estimator implementation is
# imported. Inputs are the exact binary values represented by the saved CSV.
source=joinpath(@__DIR__,"../../repair_iterations/first_full/datasets/zero_inflated_2026.csv")
rows=[split(line,',') for line in readlines(source)[2:end]]
X=[parse(Float64,rows[i][4+j]) for i in eachindex(rows),j in 1:2]
n=size(X,1)
function bandwidths(X)
    setprecision(BigFloat,512) do
        [Float64(sqrt(sum((BigFloat(x)-sum(BigFloat.(X[:,j]))/n)^2 for x in X[:,j])/n)) *
            n^(-1/6) for j in 1:2]
    end
end
function float_scores(X; factored=false)
    h=bandwidths(X)
    result=zeros(n)
    for i in 1:n
        terms=Float64[]
        for j in 1:n
            i==j && continue
            lu=-.5*((X[i,1]-X[j,1])/h[1])^2
            lv=-.5*((X[i,2]-X[j,2])/h[2])^2
            if !factored
                lu-=log(h[1])+.5log(2pi)
                lv-=log(h[2])+.5log(2pi)
            end
            push!(terms,log(1/n)+(lu+lv))
        end
        m=maximum(terms)
        result[i]=m+log(sum(sort!(exp.(terms.-m))))-log(sum(fill(1/n,n-1)))
    end
    result
end
function high_scores(X; mathematical_bandwidth=false)
    setprecision(BigFloat,512) do
        Y=BigFloat.(X)
        h= mathematical_bandwidth ?
            [sqrt(sum((Y[i,j]-sum(Y[:,j])/n)^2 for i in 1:n)/n)*BigFloat(n)^(-BigFloat(1)/6) for j in 1:2] :
            BigFloat.(bandwidths(X))
        [log(sum(exp(-sum(((Y[i,c]-Y[j,c])/h[c])^2 for c in 1:2)/2) for j in 1:n if j!=i))-
            log(BigFloat(n-1))-log(2BigFloat(pi))-sum(log.(h)) for i in 1:n]
    end
end
function dimensionless_sum(X)
    h=bandwidths(X)
    [sum(sort!([exp(-.5*(((X[i,1]-X[j,1])/h[1])^2+
                         ((X[i,2]-X[j,2])/h[2])^2)) for j in 1:n if j!=i]))/(n-1)
     for i in 1:n]
end
function midpercentiles(s)
    [(count(<(s[i]),s)+count(==(s[i]),s)/2)/n for i in 1:n]
end
open(joinpath(@__DIR__,"oracle.txt"),"w") do io
    baseline=float_scores(X)
    oracle=high_scores(X)
    ideal=high_scores(X;mathematical_bandwidth=true)
    permutation=sortperm(abs.(oracle.-oracle[41]))
    println(io,"nearest oracle scores to entity41; bandwidth=",bandwidths(X))
    for i in permutation[1:8]
        println(io,"entity=",i," xy=",X[i,:]," float_raw=",baseline[i]," high_raw=",oracle[i],
            " oracle_difference_from41=",oracle[i]-oracle[41]," ideal_difference=",ideal[i]-ideal[41])
    end
    for power in (-200,0,200)
        Y=ldexp.(X,power)
        raw=float_scores(Y)
        factored=float_scores(Y;factored=true)
        linear=dimensionless_sum(Y)
        high=high_scores(Y)
        println(io,"scale2^",power," raw41=",raw[41]," percentile41=",midpercentiles(raw)[41],
            " factored41=",factored[41]," factoredpercentile41=",midpercentiles(factored)[41],
            " highpercentile41=",midpercentiles(high)[41],
            " mismatched_float_oracle=",findall(midpercentiles(raw).!=midpercentiles(high)),
            " mismatched_factored_oracle=",findall(midpercentiles(factored).!=midpercentiles(high)))
        println(io,"linear41=",linear[41]," linear49=",linear[49]," linear_difference=",linear[41]-linear[49],
            " linearpercentile41=",midpercentiles(linear)[41],
            " mismatched_linear_oracle=",findall(midpercentiles(linear).!=midpercentiles(high)))
    end
end
print(read(joinpath(@__DIR__,"oracle.txt"),String))
