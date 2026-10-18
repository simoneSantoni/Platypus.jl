# Independent acceptance of represented scaled-mass ordering. Does not import
# package tests or mutate source. BigFloat oracle setup is serial and separate
# from the comparator, which must leave ambient precision/rounding unchanged.
using Test, SHA
module ComparatorUnderReview
include(joinpath(@__DIR__,"../../../../src/numerics.jl"))
include(joinpath(@__DIR__,"../../../../src/density.jl"))
end
const K=ComparatorUnderReview
const Q=Rational{BigInt}
q(x)=Q(x)

@testset "Scaled Gaussian kernel comparator independent acceptance" begin
    @testset "Directed rational enclosures including negative inputs" begin
        rationals=Q[0,1,-1,1//3,-1//3,(big(2)^4096+1)//3,
                    -(big(2)^4096+1)//3,1//(big(2)^4096+1),-1//(big(2)^4096+1)]
        for bits in (24,64,128,257),r in rationals
            before=(precision(BigFloat),rounding(BigFloat))
            lo,hi=K._density_rational_bounds(r,bits)
            @test Q(lo) <= r <= Q(hi)
            @test precision(lo)==bits && precision(hi)==bits
            @test (precision(BigFloat),rounding(BigFloat))==before
        end
    end

    # Two dyadic rationals bracketing e at 400 fractional bits require more
    # precision than the comparator's initial128-bit enclosure.
    below_e,above_e=setprecision(BigFloat,1024) do
        scale=big(2)^400
        lower=floor(BigInt,exp(BigFloat(1))*BigFloat(scale))//scale
        upper=lower+1//scale
        @assert BigFloat(lower)<exp(BigFloat(1))<BigFloat(upper)
        (lower,upper)
    end
    @testset "Adaptive separation and unchanged ambient context" begin
        fixtures=[(q(1),q(0),0),(q(1),q(2),-1),(q(1),q(-2),1),
                  (q(2),q(0),1),(q(1//2),q(0),-1),
                  (below_e,q(1),-1),(above_e,q(1),1),
                  (inv(below_e),q(-1),1),(inv(above_e),q(-1),-1)]
        for precision_bits in (64,256),mode in (RoundDown,RoundUp,RoundNearest)
            setprecision(BigFloat,precision_bits) do
                setrounding(BigFloat,mode) do
                    before=(precision(BigFloat),rounding(BigFloat))
                    for (ratio,target,expected) in fixtures
                        @test K._density_logratio_compare(ratio,target)==expected
                        @test (precision(BigFloat),rounding(BigFloat))==before
                    end
                end
            end
        end
    end

    @testset "One strict order across scaling representations" begin
        keys=[(q(0),0.),(q(0),-900.),(q(1),0.),(q(1),0.),
              (q(1),-1000.),(q(1//3),0.),(q(1//3),-1.),
              (q(3),-1000.),(q(1)//big(2)^2000,0.),
              (q(big(2)^2000),-1400.),(below_e,0.),(above_e,0.),
              (q(1),1.),(q(1)+q(1)//big(2)^1000,0.)]
        expected=setprecision(BigFloat,2048) do
            values=[iszero(a) ? BigFloat(-Inf) : log(BigFloat(a))+BigFloat(o) for (a,o) in keys]
            [cmp(a,b) for a in values,b in values]
        end
        actual=[K._density_kernel_compare(a,oa,b,ob) for (a,oa) in keys,(b,ob) in keys]
        @test actual==expected
        for i in eachindex(keys),j in eachindex(keys)
            @test actual[i,j]==-actual[j,i]
            @test (actual[i,j]==0)==(iszero(keys[i][1]) && iszero(keys[j][1]) || keys[i]==keys[j])
        end
        for i in eachindex(keys),j in eachindex(keys),k in eachindex(keys)
            if actual[i,j]<0 && actual[j,k]<0
                @test actual[i,k]<0
            end
        end
        original=K._density_kernel_order(first.(keys),last.(keys))
        permutation=reverse(collect(eachindex(keys)))
        @test K._density_kernel_order(first.(keys[permutation]),last.(keys[permutation]))[invperm(permutation)]==original
        results=Vector{Int}(undef,256)
        before=(precision(BigFloat),rounding(BigFloat))
        Threads.@threads for i in eachindex(results)
            results[i]=K._density_logratio_compare(isodd(i) ? below_e : above_e,q(1))
        end
        @test results==[isodd(i) ? -1 : 1 for i in eachindex(results)]
        @test (precision(BigFloat),rounding(BigFloat))==before
    end

    @testset "Weighted mixed-underflow oracle" begin
        x=[0.,100.,38.,139.]
        observation_lines=["case,entity,x,weight,normalization_order,oracle_gap_from_row4,package_percentile,oracle_percentile"]
        for (case,weights) in (("wide_weights",[1e-300,1.,0.,1e-200]),("moderate_weights",[1e-20,1.,0.,1e-10]))
            _,_,order=K._density_scores(x,zeros(4),weights;support=:real,
                bandwidth=(1.,1.),logscale=true,return_ordering=true)
            oracle=setprecision(BigFloat,2048) do
                w=BigFloat.(weights)
                [log(sum(w[j]*exp(-(BigFloat(x[i])-BigFloat(x[j]))^2/2)
                    for j in 1:4 if j!=i)/sum(w[j] for j in 1:4 if j!=i)) for i in 1:4]
            end
            # First verify the well-separated query placement, then the complete
            # oracle including the smaller reference-to-reference contrast.
            @test oracle[1]<oracle[3]<oracle[2]
            @test oracle[3]<oracle[4]
            @test order[1]<order[3]<order[2]
            @test order[3]<order[4]
            # The near reference-to-reference contrast is also resolvable here:
            # exact LOO denominator differences dominate exponentially tiny
            # distant terms. Require the complete four-row oracle, not merely
            # the originally reported query placement.
            @test [cmp(a,b) for a in order,b in order]==[cmp(a,b) for a in oracle,b in oracle]
            exact_w=Q.(weights)
            package_percentiles=[sum(exact_w[j]*(order[j]<order[i] ? 1 : order[j]==order[i] ? 1//2 : 0)
                for j in 1:4)/sum(exact_w) for i in 1:4]
            oracle_percentiles=[sum(exact_w[j]*(oracle[j]<oracle[i] ? 1 : oracle[j]==oracle[i] ? 1//2 : 0)
                for j in 1:4)/sum(exact_w) for i in 1:4]
            @test package_percentiles==oracle_percentiles
            for i in 1:4
                push!(observation_lines,join((case,i,x[i],weights[i],order[i],
                    Float64(oracle[i]-oracle[4]),K._rational_float(package_percentiles[i]),
                    K._rational_float(oracle_percentiles[i])),','))
            end
            for power in (-200,200)
                scale=2.0^power
                _,_,changed=K._density_scores(x.*scale,zeros(4),weights;support=:real,
                    bandwidth=(scale,scale),logscale=true,return_ordering=true)
                @test changed==order
            end
        end
        write(joinpath(@__DIR__,"scaled_comparator_mixed_oracle.csv"),join(observation_lines,'\n')*"\n")
    end

    @testset "Genuine41/49 distinction survives extreme common unit changes" begin
        path=joinpath(@__DIR__,"../../repair_iterations/first_full/datasets/zero_inflated_2026.csv")
        rows=[split(line,',') for line in readlines(path)[2:end]]
        X=[parse(Float64,rows[i][4+j]) for i in eachindex(rows),j in 1:2]
        original=nothing
        for power in (0,-200,200)
            Y=ldexp.(X,power)
            _,_,order=K._density_scores(Y[:,1],Y[:,2],ones(size(Y,1));support=:real,
                                       logscale=true,return_ordering=true)
            @test order[41]>order[49]
            @test (count(<(order[41]),order)+count(==(order[41]),order)/2)/96==8.5/96
            if original===nothing
                original=order
            else
                @test order==original
            end
        end
    end
end
open(joinpath(@__DIR__,"scaled_comparator_environment.txt"),"w") do io
    println(io,"Julia ",VERSION,"; threads=",Threads.nthreads())
    println(io,"density.jl SHA256=",bytes2hex(sha256(read(joinpath(@__DIR__,"../../../../src/density.jl")))))
    println(io,"No package source changes. Independent oracle precision1024/2048bits; adaptive near comparisons require>128bits.")
end
