
module RDC

using StatsBase: ecdf, mean
using LinearAlgebra

export rdc

"""
Empirical copula transformation
"""
function copula_transform(data::Matrix)
    mapslices(u->ecdf(u)(u), data, dims=2)
end

"""
Random non-linear projection

data - d x n Matrix
"""
function random_projection(data::Matrix, f::Function=sin, k::Int=20, s::Float64=1/6)
    d, n = size(data)
    data = vcat(data, ones(1, n))  
    w_t =  randn(k, d+1) * s / (d+1) # +1 for added dim
    return f.(w_t * data)
end

"""
Canonical correlation coefficient
"""

function cancor(x::Matrix, y::Matrix)
    @assert size(x, 2) == size(y, 2) "Unequal num of obs"

    xt = Matrix(x')
    yt = Matrix(y')

    # center data
    xt = xt .- mean(xt, dims=1)
    yt = yt .- mean(yt, dims=1)

    # QR decomposition of data
    Qx, Rx = qr(xt)
    Qy, Ry = qr(yt)

    dx = rank(Rx)
    dy = rank(Ry)

    @assert dx > 0 "x has rank 0"
    @assert dy > 0 "y has rank 0"

    F = svd( (Qx' * Qy)[1:dx, 1:dy]; full=:false )

    return F.S[1]
end


"""
The randomized dependence coefficient
"""
function rdc(x::Matrix, y::Matrix, f::Function=sin, k::Int=20, s::Float64=1/6)
    
    Px = copula_transform(x)
    Py = copula_transform(y)

    Φx = random_projection(Px, f, k, s)
    Φy = random_projection(Py, f, k, s)

    n = size(Φx, 2)
    Φx = vcat(Φx, ones(1, n))
    Φy = vcat(Φy, ones(1, n))

    cor = cancor(Φx, Φy)

    return cor
end

copula_transform(data::Vector) = copula_transform(reshape(data, 1, :))
random_projection(data::Vector, f::Function=sin, k::Int=10, s::Float64=1/6) = random_projection(reshape(data, 1, :), f, k, s)
rdc(x::Vector, y::Vector, f::Function=sin, k::Int=20, s::Float64=1/6) = rdc(reshape(x, 1, :), reshape(y, 1, :), f, k, s) 

end # module RDC
