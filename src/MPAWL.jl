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
getUnitCell,
getMaxIndex

include("Lattice.jl")
include("LatticeFunctions.jl")
include("FFT.jl")
include("frequencySampling.jl")
include("TISpaces.jl")
include("dlVP.jl")



end # module
