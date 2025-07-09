mutable struct Param
    mol::String
    natoms::Int64
    nbf::Int64
    nbfaux::Int64
    nalpha::Int64
    nbeta::Int64
    ne::Int64
    mul::Int64
    ndoc::Int64
    nsoc::Int64
    ndns::Int64
    nvir::Int64
    closed::Bool
    alpha::Float64
    spherical::Bool
    h_cut::Float64
    RI::Bool
end

function Param(bset, mul, charge)

    mol = ""
    natoms = bset.natoms
    nbf = bset.nbas

    ne = 0
    for i = 1:natoms
        ne += bset.atoms[i].Z
    end
    ne -= charge

    nalpha = (ne + mul - 1) / 2
    nbeta = (ne - mul + 1) / 2

    nbfaux = 0
    no1 = 0

    ndoc = nbeta - no1
    nsoc = nalpha - nbeta
    ndns = ndoc + nsoc
    nvir = nbf - nalpha

    closed = (nbeta == (ne + mul - 1) / 2 && nalpha == (ne - mul + 1) / 2)

    alpha = 0.01
    spherical = false
    h_cut = 0.02 * sqrt(2.0)
    RI = false

    return Param(
        mol,
        natoms,
        nbf,
        nbfaux,
        nalpha,
        nbeta,
        ne,
        mul,
        ndoc,
        nsoc,
        ndns,
        nvir,
        closed,
        alpha,
        spherical,
        h_cut,
        RI
    )
end