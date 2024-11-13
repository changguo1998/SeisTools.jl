using DelimitedFiles

include(abspath(@__DIR__, "../../../src/basic.jl"))
# include( "../../../src/basic.jl")

x = readdlm("IASP91.csv", ',')

flag_diff = map(axes(x, 1)) do i
    if i == 1
        return true
    end
    return !(x[i,3:end] == x[i-1,3:end])
end

y = x[flag_diff, :]

dep = y[:, 1]
r = y[:, 2]
vp = y[:, 3]
vs = y[:, 4]
rho = 1.743 .* (vp .^ 0.25)
Qk = fill(57823.0, length(r))
Qm = fill(600.0, length(r))


io1 = open("flat.bin", "w")
_write_f64_vector!(io1, dep)
_write_f64_vector!(io1, vp)
_write_f64_vector!(io1, vs)
_write_f64_vector!(io1, rho)
_write_f64_vector!(io1, Qk)
_write_f64_vector!(io1, Qm)
close(io1)

io2 = open("sphere.bin", "w")
_write_f64_vector!(io2, reverse(r))
_write_f64_vector!(io2, reverse(vp))
_write_f64_vector!(io2, reverse(vs))
_write_f64_vector!(io2, reverse(rho))
_write_f64_vector!(io2, reverse(Qk))
_write_f64_vector!(io2, reverse(Qm))
close(io2)
