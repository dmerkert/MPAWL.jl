"""
    _modM!(k,M; target="unit")

Compute k mod M, where the `target` specifies the set of congruence class
representants, i.e. either "unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function modM!{R <: AbstractFloat, I <: Integer}(
               k :: Array{R,1},
               M :: Array{I,2},
               target
              )

  @argcheck length(k) == getd(M)
  @argcheck target == "unit" || target == "symmetric"

  if target == "unit"
    k = M * mod.(M\k,1);
  else
    k = M * (mod.(M\k+0.5,1)-0.5);
  end
  k
end

modM!{I <: Integer}(
               k :: Array{I,1},
               M,
               target
              ) = round.(I,modM!(convert(Array{Float64,1},k),M,target))

"""
    h = _modM(k,M; target="unit")

Compute h = k mod M, where the `target` specifies the set of congruence class
representants, i.e. either "unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function modM(
     k,
     M,
     target
    )
  h = copy(k)
  modM!(h,M,target)
end




"""
    modM!(k,M)

Compute k mod M.
"""
modM!(k, M :: Lattice) = modM!(k,M.M,M.target) 
"""
    h = modM(k,M)

Compute h = k mod M.
"""
function modM(k, M)
  h = copy(k)
  modM!(h,M)
end


"""
    frequencyLatticeBasisDecomp(k,M; target="unit")
Decompose the vector k with respect to the frequency lattice basis of pattern matrix
M, i.e. k = modM(V.v,M), where V = frequencyLatticeBasisDecomp(M)

The `target` specifies the set of congruence class representants, i.e. either
"unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function frequencyLatticeBasisDecomp{I <: Integer}(
                                                k :: Array{I,1},
                                                M :: Lattice
                                               )
  @argcheck length(k) == M.d

  d = M.d
  dM = M.dimension
  (U,S,V) = SNF(M.M)
  # decomposing w.r.t. M^T we use the U and S of M, following the thesis
  # from above: aBV is V^-T and thats used as inverse -> transpose and
  # multiply
  v = modM(V.'*k,S,M.target)
  v = v[(d-dM+1):d]
  v
end

function getSamplingIterator(M :: Lattice)
  CartesianRange(M.size)
end

function getFrequencyIterator(M :: Lattice)
  CartesianRange(M.size)
end

function getSamplingPoint(
                          M :: Lattice,
                          coord :: CartesianIndex
                         )
  @argcheck length(coord.I) == M.dimension
  2.0*pi*modM(M.samplingLatticeBasis*(collect(coord.I)-1),
              eye(Int64,M.d),
              M.target
             )
end

function getFrequencyPoint(
                          M :: Lattice,
                          coord :: CartesianIndex
                         )
  @argcheck length(coord.I) == M.dimension
  modM(M.frequencyLatticeBasis*(collect(coord.I)-1),M.M',M.target)
end

function getUnitCell(M :: Lattice)
  [
   (0.5*M.M\[-1,-1])';
   (0.5*M.M\[-1, 1])';
   (0.5*M.M\[ 1, 1])';
   (0.5*M.M\[ 1,-1])'
  ]
end

function FFT!{R <: AbstractFloat, C <: Complex}(
              frequencyData :: Array{C},
              samplingData :: Array{R},
              M :: Lattice
             )
  @argcheck size(samplingData) == size(frequencyData)
  @argcheck size(samplingData)[1:M.dimension] == M.size

  frequencyData = fft(samplingData,1:M.dimension)
  frequencyData
end

function FFT(samplingData, M)
  frequencyData = similar(samplingData,Complex128)
  FFT!(frequencyData,samplingData,M)
end

function IFFT!{R <: AbstractFloat, C <: Complex}(
              samplingData :: Array{R},
              frequencyData :: Array{C},
              M :: Lattice
             )
  @argcheck size(samplingData) == size(frequencyData)
  @argcheck size(samplingData)[1:M.dimension] == M.size

  samplingData = real(ifft(frequencyData,1:M.dimension))
  samplingData
end

function IFFT(frequencyData, M)
  samplingData = similar(frequencyData,Float64)
  IFFT!(samplingData,frequencyData,M)
end

function setZerothFourierCoefficient!{C <: Complex}(
                                                   frequencyData :: Array{C},
                                                   M :: Lattice,
                                                   data :: Array{C}
                                                  )
  @argcheck size(frequencyData)[1:M.dimension] == M.size
  @argcheck size(frequencyData)[(M.dimension+1):end] == size(data)
  frequencyData[ones(Int,M.dimension)...,:] = data
  frequencyData
end

setZerothFourierCoefficient!{R <: AbstractFloat}(
                                                 frequencyData,
                                                 M,
                                                 data :: Array{R}
                                                ) =
setZerothFourierCoefficient!(frequencyData,M,convert(Array{Complex128},data))

