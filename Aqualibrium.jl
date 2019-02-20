# Aqualibrium.jl: Equilibrium solver for aqueous solutions
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

module aq

include("../Pitzer.jl/Pitzer.jl")
include("Dissociation.jl")

using Optim

# Get data with Pitzer
filepath = "../aqualibrium/data/"
filename = "NaCl3.csv"
mols, tempK, ions, ioncharges, data = Pitzer.IO.getCSV(filepath,filename)

# Add H-OH to ions list
ions = append!(ions,[:H,:OH])
ioncharges = Pitzer.IO.Properties.ioncharges(ions)

# Define solver functions
function Gsolve_mols(pmX,mols,charges)

    Gmols = deepcopy(mols)

    zbalance = sum(Gmols .* charges)

    if zbalance < 0 # acidic solution

        mOH = 10.0^-pmX
        Gmols[ions .== :OH] .= mOH

        mH = -sum(Gmols .* charges)
        Gmols[ions .== :H] .= mH

    else # if zbalance >= 0; basic or neutral solution

        mH = 10.0^-pmX
        Gmols[ions .== :H] .= mH

        mOH = sum(Gmols .* charges)
        Gmols[ions .== :OH] .= mOH

    end

    return Gmols, mH, mOH

end # function Gsolve_mols


function Gsolve(pmX,mols,tempK,ioncharges,lnkH2O)

    Gmols, mH, mOH = Gsolve_mols(pmX,mols,ioncharges[1])

    # Calculate activity coefficients
    lnγs = Pitzer.Model.lnγ(Gmols,tempK,ioncharges)

    lnγH  = lnγs[ions .== :H ][1]
    lnγOH = lnγs[ions .== :OH][1]

    # Evaluate equilibrium state
    GH2O = lnγH + log(mH) + lnγOH + log(mOH) - lnkH2O

    return GH2O^2

end # function Gsolve

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Start loop
for i in eachindex(tempK)

    # Calculate H2O thermodynamic dissociation constant
    lnkH2O = Dissociation.H2O(tempK[i])

    # Add H-OH to molalities
    mols[i] = append!(mols[i],[0,0])

    # Guess pH or pOH
    pmX = 7.0

    # Solve!
    Goptim = optimize(pmX -> Gsolve(pmX[1],mols[i],tempK[i],ioncharges,lnkH2O),
        [pmX],BFGS())

    println(Optim.minimizer(Goptim)[1])

    mols[i] = Gsolve_mols(Optim.minimizer(Goptim)[1],mols[i],ioncharges[1])[1]

end # for i

end # module aq
