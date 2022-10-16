# using Pkg

# try
#     Pkg.activate("example")
#     Pkg.instantiate()
# catch
#     error("Could not initialize environment.")
# end

using RDC
using StatsBase
using Plots

# line 
npoints = 500
x = rand(npoints)
y = 2*x .- 1
# plt = scatter(x, y)
# display(plt)

rdc1 = rdc(x, y)
@show rdc1

# circle
npoints = 500

theta = π*(2*rand(npoints) .- 1)
r = 1 .+ 0.05*randn(npoints)
xs, ys = 0., 0.
x = xs .+ r.*cos.(theta)
y = ys .+ r.*sin.(theta)

# plt = scatter(x, y)
# display(plt)

rdc2 = rdc(x, y, sin)
@show rdc2

# plt = scatter(plt, RDC.copula_transform(x)[:], RDC.copula_transform(y)[:])
# display(plt)
# std mv normal
npoints = 500

x = randn(npoints)
y = randn(npoints)

# plt = scatter(x, y)
# display(plt)

rdc3 = rdc(x, y)
@show rdc3
