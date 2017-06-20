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
FFT!,
FFT,
IFFT!,
IFFT,
setZerothFourierCoefficient!,
getUnitCell

include("Lattice.jl")
include("LatticeFunctions.jl")


end # module
