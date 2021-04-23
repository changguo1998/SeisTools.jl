abstract type Mechanism end

struct DoubleCouple <: Mechanism
    strike::AbstractFloat
    dip::AbstractFloat
    rake::AbstractFloat
    magnitude::AbstractFloat
    function DoubleCouple(s::AbstractFloat=0.0, d::AbstractFloat=90.0, r::AbstractFloat=0.0, m::AbstractFloat=0.0)
        if (s < 0) || (s >= 360)
            error("strike should satisfy 0<=strike<360.")
        end
        if (d < 0) || (d > 90)
            error("dip should satisfy 0<=dip<=90.")
        end
        if (r <= -180) || (r > 180)
            error("rake should satisfy -180<rake<=180.")
        end
        return new(s, d, r, m)
    end
end

function DoubleCouple(s=0, d=90, r=0, m=0)
    return DoubleCouple(Float64(s), Float64(d), Float64(r), Float64(m))
end

function DoubleCouple(;strike=0, dip=0, rake=0, mag=0)
    return DoubleCouple(Float64(strike), Float64(dip), Float64(rake), Float64(mag))
end

function DoubleCouple(sdr::AbstractArray, m::AbstractFloat=0.0)
    return DoubleCouple(sdr[1], sdr[2], sdr[3], m)
end

function DoubleCouple(sdr::Tuple, m::AbstractArray=0.0)
    return DoubleCouple(sdr[1], sdr[2], sdr[3], m)
end

function DoubleCouple(sdrm::AbstractArray)
    return DoubleCouple(sdrm[1], sdrm[2], sdrm[3], sdrm[4])
end

struct MomentTensor <: Mechanism
    value::AbstractArray
    coor::Symbol
    MomentTensor(v::AbstractArray=zeros(3, 3), coor::Symbol=:NED) = new(v, coor)
end