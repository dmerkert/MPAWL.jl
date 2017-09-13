export pyramidFunction,
       delaValleePoussinMean,
       delaValleePoussinMeanBracketSum
"""
    pyramidFunction(alpha,x)

 evaluate the d-dimensional analog of the pyramid-function, which is a
 basis for the one-dimensional de la Vallée Poussin scaling function. For
 generality, this function is evaluated w.r.t. the unit cube, i.e. from
 -0.5-alpha to 0.5+alpha, where alpha is a vector or a number (meaning it
 works in any dimension the same) and x is a point or a set of points

 INPUT
   alpha : describing the length of the linear decay, if it is a number,
           each dimension has the same decay, and each dimension has to be
           nonnegative and less or equal to 0.5
       x : a point or a matrix consisting of a point per column

 OUTPUT
       y : result(s) of the pyramid function at the point(s) x
"""
function pyramidFunction(alpha::Array{Float64,1},
                         x::Array{Float64,1}
                        )

  #= C = similar(x) =#
  #= broadcast!(pyramidFunction,C,alpha,x) =#
  #= prod(C) =#
  #= prod(pyramidFunction.(alpha,x)) =#
  y = 1.0
  for i in 1:length(x)
    y *= pyramidFunction(alpha[i],x[i])
  end
  y
end

"""
    pyramidFunction(alpha::Float64, x::Float64)

    Pyramid function in 1D
"""
@inline function pyramidFunction(
                                  alpha :: Float64,
                                  x :: Float64
                                 )
  y = 0.0

  if alpha == 0.0
    if abs(x) < 1.0/2.0
      y = 1.0
    elseif abs(x) == 1.0/2.0
      y = 1.0/2.0
    end
  else
    if abs(x) < 1.0/2.0 - alpha
      y = 1.0
    elseif (abs(x) >= 1.0/2.0 - alpha) && (abs(x) <= 1.0/2.0 + alpha)
      y = (0.5 + alpha-abs(x))/(2.0alpha)
    end
  end
  y
end

pyramidFunction(alpha::Float64,
                x::Array{Float64,1}
               ) = pyramidFunction(alpha*ones(x),x)

function delaValleePoussinMean(L :: Lattice,
                               g :: Array{Float64,1};
                               orthonormalize :: Bool = true
                              )
  @argcheck length(g) == L.d

  ind = max.(1.0+2.0g,ones(length(g)))
  tmax = getMaxDualLatticeIndex(L,ind)+1
  ckphi = zeros(Complex128,((2tmax+1)...))
  N = size(ckphi)
  for coord in getFrequencyIterator(N)
    frequency = getFrequency(N,coord)
    ckphi[coord] = pyramidFunction(g,transpose(L.M)\frequency)
  end
  ckBSq = bracketSums(L,ckphi)

  if orthonormalize
    ckphi = coeffsSpace2Fourier!(ckphi,L,1.0./sqrt.(ckBSq),ckphi)
    ckBSq = bracketSums(L,ckphi)
  end

  (ckphi,ckBSq)
end

function delaValleePoussinMean(
                               frequency :: Array{I,1},
                               L :: Lattice{I},
                               g :: Array{R,1},
                               orthonormalize :: Bool = true
                              ) :: R where {
                                       I <: Integer,
                                       R <: AbstractFloat
                                      }
  ckφ = pyramidFunction(g,L.MTFactorize\frequency)

  if orthonormalize
    bs = delaValleePoussinMeanBracketSum(
                                         frequency,
                                         L,
                                         g,
                                         false
                                        )
    if bs > 0.0
      ckφ /= sqrt(bs)
    end
  end
  ckφ
end

delaValleePoussinMean(
                      coord :: CartesianIndex{N},
                      L :: Lattice{I},
                      g :: Array{R,1},
                      orthonormalize :: Bool = true
                     ) where {
                              N,
                              I <: Integer,
                              R <: AbstractFloat
                             } = 
delaValleePoussinMean(getFrequency(N,coord),L,g,orthonormalize)

function delaValleePoussinMeanBracketSum(
                               frequency :: Array{I,1},
                               L :: Lattice{I},
                               g :: Array{R,1},
                               orthonormalize :: Bool = true
                              ) :: R where {
                                       I <: Integer,
                                       R <: AbstractFloat
                                      }
  if orthonormalize
    return 1.0
  else
    bSq = 0.0
    for i in BracketSumIterator(frequency, -2, 2, L)
      bSq += abs2(pyramidFunction(g,L.MTFactorize\i))
    end
    bSq
  end
end
