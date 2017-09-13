export bracketSums,
       coeffsSpace2Fourier!,
       changeBasis!


function coeffsSpace2Fourier!{R <: Union{AbstractFloat,Complex},
                              C <: Complex
                             }(
                               data :: Array{C},
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

#= function changeBasis!{ =#
#=                       R <: AbstractFloat, =#
#=                       R1 <: Union{AbstractFloat,Complex}, =#
#=                       R2 <: Union{AbstractFloat,Complex}, =#
#=                       I <: Integer, =#
#=                       N, =#
#=                       M =#
#=                      }( =#
#=                        inputData :: Array{R1,N}, =#
#=                        outputData :: Array{R2,N}, =#
#=                        L :: Lattice, =#
#=                        bracketSums :: Array{R,M}, =#
#=                        dims :: Array{I,1}; =#
#=                        inputDomain :: String = "Space", =#
#=                        outputDomain :: String = "Space" =#
#=                       ) =#
#=   @argcheck (inputDomain == "Space") || (inputDomain == "Fourier") =#
#=   @argcheck (outputDomain == "Space") || (outputDomain == "Fourier") =#
#=   @argcheck size(inputData) == L.size =#
#=   @argcheck size(outputData) == L.size =#
#=   @argcheck size(bracketSums) == L.size =#

#=   dataFourier = (inputDomain == "Space") ? =#
#=     patternfft(inputData,L,dims) : =#
#=     inputData =#

#=   dataFourier .*= 1.0/L.m.*bracketSums =#

#=   outputData = (outputDomain == "Space") ? =#
#=     patternifft(dataFourier,L,dims) : =#
#=     dataFourier =#
#= end =#

################ NEW #############################
function changeBasis!(
                      inputData :: Array{R1,N},
                      outputData :: Array{R2,N},
                      L :: Lattice,
                      bracketSum :: Function,
                      dims :: Array{I,1};
                      inputDomain :: String = "Space",
                      outputDomain :: String = "Space"
                     ) where {
                              R1 <: Union{AbstractFloat,Complex},
                              R2 <: Union{AbstractFloat,Complex},
                              I <: Integer,
                              N
                             }

  @argcheck (inputDomain == "Space") || (inputDomain == "Fourier")
  @argcheck (outputDomain == "Space") || (outputDomain == "Fourier")
  @argcheck size(inputData) == L.size
  @argcheck size(outputData) == L.size

  dataFourier = (inputDomain == "Space") ?
    patternfft(inputData,L,dims) :
    inputData

    tmp = zeros(Float64,L.d)
    point = zeros(I,L.d)
    coordExtended = Array{Union{Int,Colon},1}(N)
    coordExtended .= Colon()
    for coord in getFrequencyIterator(L)
      h = getFrequencyPoint!(L,coord,point,tmp)
      coordExtended[dims] = collect(coord.I)
      dataFourier[coordExtended] .*= 1.0/L.m*bracketSum(h)
    end

  outputData = (outputDomain == "Space") ?
    patternifft(dataFourier,L,dims) :
    dataFourier
end

