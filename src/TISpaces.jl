export bracketSums,
       coeffsSpace2Fourier!,
       changeBasis!

function bracketSums{C <: Complex}(L :: Lattice,
                                        data :: Array{C}
                                       )
  @argcheck length(size(data)) == L.d

  N = size(data)
  bSq = zeros(Float64,L.size)
  LUnit = Lattice(L.M,target="unit") :: Lattice{Int}
  bSq = _bracketSums!(LUnit,data,bSq,N,CartesianRange(N))
  bSq
end

function _bracketSums!{I,R,C,N1,N2}(LUnit :: Lattice{I},
                                    data :: Array{C,N1},
                                    bSq :: Array{R,N2},
                                    N,
                                    range :: CartesianRange
                                   )

  for coord in range
    frequency = getFrequency(N,coord)
    h = frequencyLatticeBasisDecomp(frequency,LUnit)
    bSq[(h+1)...] :: R += abs2(data[coord]) :: R
  end
  bSq
end

function coeffsSpace2Fourier!{R <: AbstractFloat,
                              C <: Complex}(data :: Array{C},
                                            L :: Lattice,
                                            hata :: Array{R},
                                            ckPhi :: Array{C}
                                           )
  @argcheck size(hata) == L.size
  @argcheck size(ckPhi) == size(data)

  N = size(data)
  LUnit = Lattice(L.M,target="unit")
  for coord in getFrequencyIterator(N)
    frequency = getFrequency(N,coord)
    h = frequencyLatticeBasisDecomp(frequency,LUnit)
    data[coord] = hata[(h+1)...] * ckPhi[coord]
  end
  data
end

function changeBasis!{
                      R <: AbstractFloat,
                      R1 <: Union{AbstractFloat,Complex},
                      R2 <: Union{AbstractFloat,Complex},
                      I <: Integer,
                      N,
                      M
                     }(
                       inputData :: Array{R1,N},
                       outputData :: Array{R2,N},
                       L :: Lattice,
                       bracketSums :: Array{R,M},
                       dims :: Array{I,1};
                       inputDomain :: String = "Space",
                       outputDomain :: String = "Space"
                      )
  @argcheck (inputDomain == "Space") || (inputDomain == "Fourier")
  @argcheck (outputDomain == "Space") || (outputDomain == "Fourier")
  @argcheck size(inputData) == L.size
  @argcheck size(outputData) == L.size
  @argcheck size(bracketSums) == L.size

  dataFourier = (inputDomain == "Space") ?
    FFT(inputData,L,dims) :
    inputData

  dataFourier .*= 1.0/L.m.*bracketSums

  outputData = (outputDomain == "Space") ?
    IFFT!(outputData,dataFourier,L,dims) :
    dataFourier
end
