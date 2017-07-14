export FFT!,
       FFT,
       IFFT!,
       IFFT,
       setZerothFourierCoefficient!,
       setZerothFourierCoefficient,
       FirstDimensionsFFT,
       LastDimensionsFFT


       abstract type FirstDimensionsFFT end
       abstract type LastDimensionsFFT end



function FFT!{R <: AbstractFloat, C <: Complex}(
              frequencyData :: Array{C},
              samplingData :: Array{R},
              M :: Lattice,
              t :: Type{FirstDimensionsFFT}
             )
  @argcheck size(samplingData) == size(frequencyData)
  dims = 1:M.dimension
  @argcheck size(samplingData)[dims] == M.size

  frequencyData = fft(samplingData,dims)
  frequencyData
end

function FFT!{R <: AbstractFloat, C <: Complex}(
              frequencyData :: Array{C},
              samplingData :: Array{R},
              M :: Lattice,
              t :: Type{LastDimensionsFFT}
             )
  @argcheck size(samplingData) == size(frequencyData)
  dims = (length(size(samplingData))-M.dimension+1):length(size(samplingData))
  @argcheck size(samplingData)[dims] == M.size

  frequencyData = fft(samplingData,dims)
  frequencyData
end

function FFT(samplingData, M, t)
  frequencyData = similar(samplingData,Complex128)
  FFT!(frequencyData,samplingData,M,t)
end

function IFFT!{R <: AbstractFloat, C <: Complex}(
              samplingData :: Array{R},
              frequencyData :: Array{C},
              M :: Lattice,
              t :: Type{FirstDimensionsFFT}
             )
  @argcheck size(samplingData) == size(frequencyData)
  dims = 1:M.dimension
  @argcheck size(samplingData)[dims] == M.size

  samplingData = real(ifft(frequencyData,dims))
  samplingData
end

function IFFT!{R <: AbstractFloat, C <: Complex}(
              samplingData :: Array{R},
              frequencyData :: Array{C},
              M :: Lattice,
              t :: Type{LastDimensionsFFT}
             )
  @argcheck size(samplingData) == size(frequencyData)
  dims = (length(size(samplingData))-M.dimension+1):length(size(samplingData))
  @argcheck size(samplingData)[dims] == M.size

  samplingData = real(ifft(frequencyData,dims))
  samplingData
end

function IFFT(frequencyData, M, t)
  samplingData = similar(frequencyData,Float64)
  IFFT!(samplingData,frequencyData,M, t)
end

function setZerothFourierCoefficient!{C <: Complex}(
                                                   frequencyData :: Array{C},
                                                   M :: Lattice,
                                                   data :: Array{C},
                                                   t :: Type{FirstDimensionsFFT}
                                                  )
  dims = 1:M.dimension
  dataDimensions = (dims[end]+1):length(size(frequencyData))
  @argcheck size(frequencyData)[dims] == M.size
  @argcheck size(frequencyData)[dataDimensions] == size(data)
  frequencyData[ones(Int,length(dims))...,:] = data
  frequencyData
end

function setZerothFourierCoefficient!{C <: Complex}(
                                                   frequencyData :: Array{C},
                                                   M :: Lattice,
                                                   data :: Array{C},
                                                   t :: Type{LastDimensionsFFT}
                                                  )
  dims = (length(size(frequencyData))-M.dimension+1):length(size(frequencyData))
  dataDimensions = 1:(dims[1]-1)
  @argcheck size(frequencyData)[dims] == M.size
  @argcheck size(frequencyData)[dataDimensions] == size(data)
  frequencyData[:,ones(Int,length(dims))...] = data
  frequencyData
end

setZerothFourierCoefficient!{R <: AbstractFloat}(
                                                 frequencyData,
                                                 M,
                                                 data :: Array{R},
                                                 t
                                                ) =
setZerothFourierCoefficient!(frequencyData,M,convert(Array{Complex128},data),t)

