"""
 We follow a similar structure as used in ITensorMPS.jl
 https://github.com/ITensor/ITensorMPS.jl/blob/main/src/lattices/lattices.jl
"""


"""
A LatticeBond is a struct which represents
a single bond in a geometrical lattice or
else on interaction graph defining a physical
model such as a quantum Hamiltonian.

LatticeBond has the following data fields:

  - s1::Int -- number of site 1
  - s2::Int -- number of site 2
  - x1::Float64 -- x coordinate of site 1
  - y1::Float64 -- y coordinate of site 1
  - x2::Float64 -- x coordinate of site 2
  - y2::Float64 -- y coordinate of site 2
  - type::String -- optional description of bond type
"""
struct LatticeBond
  s1::Int
  s2::Int
  x1::Float64
  y1::Float64
  x2::Float64
  y2::Float64
  type::String
end

"""
    LatticeBond(s1::Int,s2::Int)

    LatticeBond(s1::Int,s2::Int,
                x1::Real,y1::Real,
                x2::Real,y2::Real,
                type::String="")

Construct a LatticeBond struct by
specifying just the numbers of sites
1 and 2, or additional details including
the (x,y) coordinates of the two sites and
an optional type string.
"""
function LatticeBond(s1::Int, s2::Int)
  return LatticeBond(s1, s2, 0.0, 0.0, 0.0, 0.0, "")
end

#function LatticeBond(s1::Int, s2::Int, bondtype::String="")
#  return LatticeBond(s1, s2, 0.0, 0.0, 0.0, 0.0, bondtype)
#end

function LatticeBond(
  s1::Int, s2::Int, x1::Real, y1::Real, x2::Real, y2::Real, bondtype::String=""
)
  cf(x) = convert(Float64, x)
  return LatticeBond(s1, s2, cf(x1), cf(y1), cf(x2), cf(y2), bondtype)
end

"""
Lattice is an alias for Vector{LatticeBond}
"""
const Lattice = Vector{LatticeBond}

