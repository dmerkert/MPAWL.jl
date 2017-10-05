export patternfft!,
patternfft,
patternifft!,
patternifft,
setFourierCoefficient!
"""
    patternfft!(d, L, [dims])

computes the pattern discrete fast Fourier transformation on a pattern of the
lattice L. An array d is transformed in the dimensions specified by dims,
such that the sizes of these dimensions (in the order they are given in dims) correspond to the elementary divisors of the pattern matrix. If no dims are given, all of d are taken into account

# Example
for an array d of size 8x8 and a lattice of the matrix [8,0 0,8] calling

    pattern!(d,L)

performs the clacssical 2D discrete Fourier transform of the 8x8 array in place.
"""
function patternfft!{
              R <: Union{Complex,AbstractFloat},
              I <: Integer,
              N
             }(
               d :: Array{R,N},
               L :: Lattice,
               dims :: Array{I,1} = collect(1:ndims(d))
              )
  @argcheck size(d)[dims] == L.size
  fft!(d,dims)
end
"""
    patternfft(d, L, [dims])

computes the pattern discrete fast Fourier transformation on a pattern of the
lattice L. An array d is transformed in the dimensions specified by dims,
such that the sizes of these dimensions (in the order they are given in dims) correspond to the elementary divisors of the pattern matrix. If no dims are given, all of d are taken into account.

# Example
for an array d of size 8x8 and a lattice of the matrix [8,0 0,8] calling

    dhat = pattern(d,L)

performs the clacssical 2D discrete Fourier transform of the 8x8 array.
"""
function patternfft{
              R <: Union{Complex,AbstractFloat},
              I <: Integer,
              N
             }(
               d :: Array{R,N},
               L :: Lattice,
               dims :: Array{I,1} = collect(1:ndims(d))
              )
  @argcheck size(d)[dims] == L.size
  return fft(d,dims)
end
"""
patternifft!(dhat, L, [dims])

computes the pattern discrete fast Fourier transformation on a pattern of the
lattice L in place. An array d is transformed in the dimensions specified by
dims, such that the sizes of these dimensions (in the order they are given in
dims) correspond to the elementary divisors of the pattern matrix. If no dims
are given, all of d are taken into account.
"""
function patternifft!{
              R <: Union{Complex,AbstractFloat},
              I <: Integer,
              N
             }(
               dhat :: Array{R,N},
               L :: Lattice,
               dims :: Array{I,1} = collect(1:ndims(dhat))
              )
  @argcheck size(dhat)[dims] == L.size
  ifft!(dhat,dims)
end
"""
patternifft(dhat, L, [dims])

computes the pattern discrete fast Fourier transformation on a pattern of the
lattice L. An array d is transformed in the dimensions specified by
dims, such that the sizes of these dimensions (in the order they are given in
dims) correspond to the elementary divisors of the pattern matrix. If no dims
are given, all of d are taken into account.
"""
function patternifft{
              R <: Union{Complex,AbstractFloat},
              I <: Integer,
              N
             }(
               dhat :: Array{R,N},
               L :: Lattice,
               dims :: Array{I,1} = collect(1:ndims(dhat))
              )
  @argcheck size(dhat)[dims] == L.size
  return ifft(dhat,dims)
end
"""
  setFourierCoefficient!(dhat,L,value,[k,dims])

set the kth (0th) Fourier coefficient in dhat
to value, where the Fourier dimensions are specified by dims.
If no dimension is given, the value is complex-valued
"""
function setFourierCoefficient{
                                C <: Complex,
                                I <: Integer,
                                N,
                               }(
                                 dhat :: Array{C,N},
                                 L :: Lattice,
                                 value :: C,
                                 k :: Array{I,1} = ones(I,length(dims))
                                )
    print(k)
    setFourierCoefficient(dhat,L,value,[k],collect(1:ndims(dhat)))
end
function setFourierCoefficient!{
                                C <: Complex,
                                I <: Integer,
                                N,
                                M
                               }(
                                 dhat :: Array{C,N},
                                 L :: Lattice,
                                 value :: Array{C,M},
                                 k :: Array{I,1} = ones(I,length(dims)),
                                 dims :: Array{I,1} = collect(1:ndims(dhat))
                                )
  @argcheck size(dhat)[dims] == L.size
  d = length(size(dhat))
  dataDimensions = [i for i=1:d if findfirst(dims,i) == 0]

  @argcheck size(dhat)[dataDimensions] == size(value)
  @argcheck length(k) == length(dims)

  beginIndex = ones(I,d)
  beginIndex[dims] = k
  endIndex = ones(I,d)
  endIndex[dims] = k
  endIndex[dataDimensions] = collect(size(dhat)[dataDimensions])

  for coord in CartesianRange(
                              CartesianIndex((beginIndex...)),
                              CartesianIndex((endIndex...))
                             )
    valueIndex = CartesianIndex(coord.I[dataDimensions])
    dhat[coord] = value[valueIndex]
  end
  dhat
end
function setFourierCoefficient!{
                       C <: Complex,
                       R <: AbstractFloat,
                       I <: Integer,
                       N,
                       M
                      }(
                        dhat :: Array{C,N},
                        L :: Lattice,
                        value :: Array{R,M},
                        k :: Array{I,1} = ones(I,length(dims)),
                        dims :: Array{I,1} = collect(1:ndims(dhat))
                       )
  setFourierCoefficient!(dhat,L,convert(Array{C,M},value),k,dims)
end
