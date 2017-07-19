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

function changeBasis!{R <: AbstractFloat}(data :: Array{R},
                                          L :: Lattice,
                                          bracketSums :: Array{R},
                                          t
                                         )
  @argcheck size(data) == L.size
  @argcheck size(bracketSums) == L.size

  hatData = FFT(data,L,t)
  hatData .*= 1.0/L.m.*bracketSums
  IFFT!(data,hatData,L,t)
end

function changeBasis!{R <: AbstractFloat, C <: Complex}(hatData :: Array{C},
                                                        L :: Lattice,
                                                        bracketSums :: Array{R},
                                                        t
                                                       )
  @argcheck size(hatData) == L.size
  @argcheck size(bracketSums) == L.size

  hatData .*= 1.0/L.m.*bracketSums
end

function changeBasis!{R <: AbstractFloat, C <: Complex}(inputData :: Array{R},
                                                        outputData :: Array{C},
                                                        L :: Lattice,
                                                        bracketSums :: Array{R},
                                                        t
                                                       )
  @argcheck size(inputData) == L.size
  @argcheck size(outputData) == L.size
  @argcheck size(bracketSums) == L.size

  outputData = FFT!(outputData,inputData,L,t)
  outputData .*= 1.0/L.m.*bracketSums
end

function changeBasis!{R <: AbstractFloat, C <: Complex}(inputData :: Array{C},
                                                        outputData :: Array{R},
                                                        L :: Lattice,
                                                        bracketSums :: Array{R},
                                                        t
                                                       )
  @argcheck size(inputData) == L.size
  @argcheck size(outputData) == L.size
  @argcheck size(bracketSums) == L.size

  inputData .*= 1.0/L.m.*bracketSums
  IFFT!(outputData,inputData,L,t)
end
