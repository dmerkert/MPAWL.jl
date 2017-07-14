module MPAWL

using ArgCheck
using IntegerSmithNormalForm
using Base.fft

export modM!,
modM,
Lattice,
frequencyLatticeBasisDecomp,
getSamplingPoint,
getFrequencyPoint,
getSamplingIterator,
getFrequencyIterator,
getUnitCell

include("Lattice.jl")
include("LatticeFunctions.jl")
include("FFT.jl")


end # module
