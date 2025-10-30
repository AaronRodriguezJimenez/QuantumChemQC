"""
    square_lattice(Nx::Int,
                   Ny::Int;
                   kwargs...)::Lattice

Return a Lattice (array of LatticeBond
objects) corresponding to the two-dimensional
square lattice of dimensions (Nx,Ny).
By default the lattice has open boundaries,
but can be made periodic in the y direction
by specifying the keyword argument
`yperiodic=true`.
"""
function square_lattice(Nx::Int, Ny::Int; yperiodic=false)::Lattice
  yperiodic = yperiodic && (Ny > 2)
  N = Nx * Ny
  Nbond = 2N - Ny + (yperiodic ? 0 : -Nx)
  latt = Lattice(undef, Nbond)
  b = 0
  for n in 1:N
    x = div(n - 1, Ny) + 1
    y = mod(n - 1, Ny) + 1
    if x < Nx
      latt[b += 1] = LatticeBond(n, n + Ny, x, y, x + 1, y)
    end
    if Ny > 1
      if y < Ny
        latt[b += 1] = LatticeBond(n, n + 1, x, y, x, y + 1)
      end
      if yperiodic && y == 1
        latt[b += 1] = LatticeBond(n, n + Ny - 1, x, y, x, y + Ny - 1)
      end
    end
  end
  return latt
end

"""
    triangular_lattice(Nx::Int,
                       Ny::Int;
                       kwargs...)::Lattice

Return a Lattice (array of LatticeBond
objects) corresponding to the two-dimensional
triangular lattice of dimensions (Nx,Ny).
By default the lattice has open boundaries,
but can be made periodic in the y direction
by specifying the keyword argument
`yperiodic=true`.
"""
function triangular_lattice(Nx::Int, Ny::Int; yperiodic=false)::Lattice
  yperiodic = yperiodic && (Ny > 2)
  N = Nx * Ny
  Nbond = 3N - 2Ny + (yperiodic ? 0 : -2Nx + 1)
  latt = Lattice(undef, Nbond)
  b = 0
  for n in 1:N
    x = div(n - 1, Ny) + 1
    y = mod(n - 1, Ny) + 1

    # x-direction bonds
    if x < Nx
      latt[b += 1] = LatticeBond(n, n + Ny)
    end

    # 2d bonds
    if Ny > 1
      # vertical / y-periodic diagonal bond
      if (n + 1 <= N) && ((y < Ny) || yperiodic)
        latt[b += 1] = LatticeBond(n, n + 1)
      end
      # periodic vertical bond
      if yperiodic && y == 1
        latt[b += 1] = LatticeBond(n, n + Ny - 1)
      end
      # diagonal bonds
      if x < Nx && y < Ny
        latt[b += 1] = LatticeBond(n, n + Ny + 1)
      end
    end
  end
  return latt
end