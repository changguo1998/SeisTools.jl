# Source

```@meta
CurrentModule = SeisTools.Source
```

Module including the representation and calculation functions on source.
Currently there are two kinds of representation: `SDR`(fault plane solution) and
`MomentTensor`

## Moment Tensor

Moment tensor can be expressed as an 3 by 3 matrix:

```math
M = \begin{pmatrix}
m11 & m12 & m13 \\
m12 & m22 & m23 \\
m13 & m23 & m33
\end{pmatrix}
```

There are several ways to initial a `MomentTensor` object:

- using a vector or tuple with 6 elements or input them seperately
- using a 3 by 3 matrix of the total moment tensor. the matrix must be symetric
- using `strike, dip, rake, m0` to initial a `MomentTensor` with only double-couple component or
  convert from a `SDR` object

Available functions are:

- [M0](@ref M0) calculate scalar moment of the moment tensor using formular
  ```math
  M_0=\frac{1}{\sqrt{2}}\sqrt{\sum_{i,j=1,2,3}m_{ij}^2}
  ```
- [decompose](@ref decompose) decompose the full moment tensor into ISO, DC and CLVD component
- [kagan](@ref kagan(mt1::MomentTensor, mt2::MomentTensor)) calculate the Kagan angle between two
  `MomentTensor`s' principle axis
- [beachball_sdrline](@ref beachball_sdrline(m::MomentTensor, dtheta::Real = 1.0; innerdecompose::Bool = true))
  return lines to draw the edge of beachball of the DC component of `m`
- [beachball_bitmap](@ref beachball_bitmap(m::MomentTensor; resolution::Tuple{<:Integer,<:Integer} = (201, 201)))
  return a `Bool` matrix of beachball.

```@docs
MomentTensor(m11::Real, m22::Real, m33::Real, m12::Real, m13::Real, m23::Real)
MomentTensor(x::Union{AbstractVector{<:Real},Tuple{<:Real,<:Real,<:Real,<:Real,<:Real,<:Real}})
MomentTensor(x::AbstractMatrix{<:Real})
MomentTensor(strike::Real, dip::Real, rake::Real; m0::Real = 1.0)
MomentTensor(sdr::SDR)
M0
decompose
kagan(mt1::MomentTensor, mt2::MomentTensor)
beachball_sdrline(m::MomentTensor, dtheta::Real = 1.0; innerdecompose::Bool = true)
beachball_bitmap(m::MomentTensor; resolution::Tuple{<:Integer,<:Integer} = (201, 201))
```
