# Aqualibrium.jl: Equilibrium solver for aqueous solutions
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

module Dissociation

# M88
H2O(T) =
     1.04031130e+3 +
     4.86092851e-1 * T +
    -3.26224352e+4 / T +
    -1.90877133e+2 * log(T) +
    -5.35204850e-1 / (T - 263) +
    -2.32009393e-4 * T^2 +
     5.20549183e+1 / (680 - T)

end # module Dissociation
