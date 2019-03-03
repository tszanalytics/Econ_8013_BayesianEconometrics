# gsres function example

## data here:
using Distributions, Plots, StatPlots
n = 100
x = randn(n)
y = 1.0 .+ 1.0.*x .+ 0.2.*randn(n)
plot(x,y,st=:scatter)


include("gsreg.jl")

## gsreg(y::Array{Float64},X::Array{Float64}; M=10000::Int64, burnin = Int(floor(M/10.0))::Int64,tau=[1.0]::Array{Float64},iB0=[0.0001]::Array{Float64},b0=[0.0]::Array{Float64},d0=0.0001::Float64, a0=3.0::Float64)
# uninformative prior
X = [ones(n) x]
bdraws,s2draws = gsreg(y,X)

plot(bdraws[:,2],st=:density)
mean(bdraws[:,2])
std(bdraws[:,2])
quantile(bdraws[:,2],[0.025,0.975])


# Informative prior (still uninformative for variance parameter)
b = [0.0; 1.0]    # prior coeff. means
iB = inv([1000.0 0.0; 0.0 0.001])  # prior cov matrix
bdraws,s2draws = gsreg(y,X, b0=b, iB0=iB)

plot!(bdraws[:,2],st=:density,linecolor=:red,label="informative")
mean(bdraws[:,2])
std(bdraws[:,2])
quantile(bdraws[:,2],[0.025,0.975])
