export FFT!,
FFT,
IFFT!,
IFFT,
setFourierCoefficient!

function FFT!{
              R <: Union{Complex,AbstractFloat},
              C <: Complex,
              I <: Integer,
              N
             }(
               frequencyData :: Array{C,N},
               samplingData :: Array{R,N},
               L :: Lattice,
               dims :: Array{I,1}
              )
  @argcheck size(samplingData) == size(frequencyData)
  @argcheck size(samplingData)[dims] == L.size

  frequencyData = fft(samplingData,dims)
  frequencyData
end

function FFT(samplingData, L, dims)
  frequencyData = similar(samplingData,Complex128)
  FFT!(frequencyData,samplingData,L,dims)
end

function IFFT!{
               R <: AbstractFloat,
               C <: Complex,
               I <: Integer,
               N
              }(
                samplingData :: Array{R,N},
                frequencyData :: Array{C,N},
                L :: Lattice,
                dims :: Array{I,1}
               )
  @argcheck size(samplingData) == size(frequencyData)
  @argcheck size(samplingData)[dims] == L.size

  samplingData = real(ifft(frequencyData,dims))
  samplingData
end

function IFFT!{
               C <: Complex,
               I <: Integer,
               N
              }(
                samplingData :: Array{C,N},
                frequencyData :: Array{C,N},
                L :: Lattice,
                dims :: Array{I,1}
               )
  @argcheck size(samplingData) == size(frequencyData)
  @argcheck size(samplingData)[dims] == L.size

  samplingData = ifft(frequencyData,dims)
  samplingData
end

function IFFT(frequencyData, L, dims; realData :: Bool = true)
  if realData
    samplingData = similar(frequencyData,Float64)
    return IFFT!(samplingData,frequencyData,L, dims)
  else
    samplingData = similar(frequencyData,Complex128)
    return IFFT!(samplingData,frequencyData,L, dims)
  end
end

function setFourierCoefficient!{
                                C <: Complex,
                                I <: Integer,
                                N,
                                M
                               }(
                                 frequencyData :: Array{C,N},
                                 L :: Lattice,
                                 data :: Array{C,M},
                                 dims :: Array{I,1},
                                 k :: Array{I,1} = ones(I,length(dims))
                                )

  @argcheck size(frequencyData)[dims] == L.size

  d = length(size(frequencyData))
  dataDimensions = [i for i=1:d if findfirst(dims,i) == 0]

  @argcheck size(frequencyData)[dataDimensions] == size(data)
  @argcheck length(k) == length(dims)

  beginIndex = ones(I,d)
  beginIndex[dims] = k
  endIndex = ones(I,d)
  endIndex[dims] = k
  endIndex[dataDimensions] = collect(size(frequencyData)[dataDimensions])

  for coord in CartesianRange(
                              CartesianIndex((beginIndex...)),
                              CartesianIndex((endIndex...))
                             )
    dataIndex = CartesianIndex(coord.I[dataDimensions])
    frequencyData[coord] = data[dataIndex]
  end
  frequencyData
end

setFourierCoefficient!{
                       C <: Complex,
                       R <: AbstractFloat,
                       I <: Integer,
                       N,
                       M
                      }(
                        frequencyData :: Array{C,N},
                        L :: Lattice,
                        data :: Array{R,M},
                        dims :: Array{I,1},
                        k :: Array{I,1} = ones(I,length(dims))
                       ) = setFourierCoefficient!(frequencyData,
                                                  L,
                                                  convert(Array{C,M},data),
                                                  dims,
                                                  k)

