using IntegerSmithNormalForm

export Lattice
"""
    Lattice

Lattice represents a regular integer matrix parametrizing a regular
anisotropic pattern on the torus. Has field `M` containing the matrix.


"""
immutable Lattice{I <: Integer, MF, MF2}
  M :: Array{I,2}
  MFactorize :: MF
  MTFactorize :: MF2
  target :: String
  d :: I
  m :: I
  size :: Tuple
  rank :: I
  samplingLatticeBasis :: Array{Float64,2}
  frequencyLatticeBasis :: Array{I,2}
  patternNormalForm :: Array{I,2}
  SNF :: Tuple{Array{I,2},Array{I,2},Array{I,2}}

  function Lattice{I}(M :: Array{I,2}; target="symmetric") where I <: Integer
    @argcheck size(M,1) == size(M,2)
    @argcheck det(M) != 0
    @argcheck target == "unit" || target == "symmetric"

    _MF = factorize(M)
    _MFT = factorize(M')

    d = getd(M)
    m = getm(M)
    _SNF = SNF(M)
    _patternSize = patternSize(M,_SNF)
    _patternRank = patternRank(M,_SNF)
    _patternBasis = patternBasis(M,target,_SNF)
    _generatingSetBasis = generatingSetBasis(M.',target,_SNF)
    _patternNormalForm = patternNormalForm(M)

    @assert d > 0
    @assert m > 0
    @assert _patternRank <= d


    #= S = diag(_SNF[2]) =#
    #=  for i in 1:size(_patternBasis,2) =#
    #=    for j in 1:size(_patternBasis,2) =#
    #=      y = _patternBasis[:,i] =#
    #=      h = _generatingSetBasis[:,j] =#

    #=      i == j && @assert abs( =#
    #=                            mod(dot(y,h),1.0) - 1.0/S[d-_patternRank+i] =#
    #=                           )/abs(1.0/S[d-_patternRank+i]) < 1e-10 =#
    #=      i != j && @assert ( =#
    #=                         ( =#
    #=                          abs( =#
    #=                              mod(dot(y,h),1.0) =#
    #=                             ) < 1e-12 =#
    #=                         ) || =#
    #=                         ( =#
    #=                          abs( =#
    #=                              mod(dot(y,h),1.0) - 1.0 =#
    #=                             ) < 1e-12 =#
    #=                         ) =#
    #=                        ) =#
    #=    end =#
    #=  end =#

    new{I,typeof(_MF),typeof(_MFT)}(
        M,
        _MF,
        _MFT,
        target,
        d,
        m,
        _patternSize,
        _patternRank,
        _patternBasis,
        _generatingSetBasis,
        _patternNormalForm,
        _SNF
       )
  end
end
"""
    Lattice(M): constructs a lattice the integer matrix `M`
"""
Lattice{I}(M :: Array{I,2}; target="symmetric") =
Lattice{I}(M,target=target)
"""
    Lattice(v): constructs a diagonal pattern matrix from the integer
    vector `v`
"""
Lattice{I}(v :: Array{I,1}; target="symmetric") =
Lattice{I}(diagm(v),target=target)

"""
    m = getm(M): calculates the number of pattern points given by m = abs(det(M))
"""
function getm{I <: Integer}(M :: Array{I,2})
  round(I,abs(det(M)))
end

"""
    d = getd(M)

returns the dimension d = size(M,1) of the pattern
"""
function getd{I <: Integer}(M :: Array{I,2})
  size(M,1)
end

"""
    s = patternSize(M)
 Returns the elementary divisors of M, that are greater than 1.
 This corresponds to the size of the lattice.
 """
 function patternSize{I <: Integer}(M :: Array{I,2},
                                    SNF :: Tuple{Array{I,2},Array{I,2},Array{I,2}}
                                   )
   d = diag(SNF[2])
   s = tuple(d[d .> 1]...)
end

"""
    dM = patternRank(M)
 Returns the number of elementary divisors of M, that are greater than 1.
 This corresponds to the number of basis vectors of the pattern and hence
 the dimension of the corresponding lattice.
"""
function patternRank{I <: Integer}(M :: Array{I,2},
                                   SNF :: Tuple{Array{I,2},Array{I,2},Array{I,2}}
                                  )
  sum(diag(SNF[2]) .> 1)
end

"""
    V = patternBasis(M; target="unit")
 creates a set of vectors (columns of V), whose integer multiples (up to
 the corresponding elementary divisor - 1) create the pattern of M. The vectors
 are ordered w. r. t. the ordering of the elementary divisors, i.e. w.r.t.
 their nondecreasing ordering from the diagonal of the Smith normal form of M.
 V is a d x patternDimension(M) matrix.

 The `target` specifies the set of congruence class representants, i.e. either
 "unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
 function patternBasis{I <: Integer}(M :: Array{I,2},
                                     target,
                                     SNF :: Tuple{Array{I,2},Array{I,2},Array{I,2}}
   )
   @argcheck target == "unit" || target == "symmetric"

   d = getd(M)
   rank = patternRank(M,SNF)
   _V = SNF[3]*(diagm(1./diag(SNF[2])))
   _V = _V[:,d-rank+1:d]
   for i in 1:rank
     _V[:,i] = modM(_V[:,i],eye(I,d),target)
   end
   _V
 end

"""
    generatingSetBasis(M; target="unit")
 creates a set of vectors (columns of V), whose integer multiples (up to the
 corresponding elementary divisor - 1) create the generating set of M. The
 vectors are ordered w. r. t. the ordering of the elementary divisors, i.e.
 w.r.t. their nondecreasing ordering from the diagonal of the Smith normal form
 of M. V is a d x patternDimension(M) matrix.

 The `target` specifies the set of congruence class representants, i.e. either
 "unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function generatingSetBasis{I <: Integer}(M :: Array{I,2},
                                          target,
                                          _SNF :: Tuple{Array{I,2},Array{I,2},Array{I,2}}
  )
  @argcheck target == "unit" || target == "symmetric"

  d = getd(M)
  rank = patternRank(M,_SNF)
  (U,S,V) = SNF(transpose(M))

  _V = transpose(inv(V))
  _V = _V[:,d-rank+1:d]
  for i = 1:rank
    _V[:,i] = modM(_V[:,i],M,target)
  end
  round.(I,_V)
end


"""
    patternNormalform(M)

computes the pattern normal form of the patter matrix M.

The pattern normal form of the matrix, i.e. a matrix obtained by
computing the gaussian elimination (in integers) to obtain an upper
triangular matrix, where all elements above the diagonal are smaller
than the diagonal element
"""
function patternNormalForm{I <: Integer}(M :: Array{I,2})
  pMf = copy(M)
  d = getd(M)
  # form upper triangular matrix
  for col in 1:(d-1)
    for row in (col+1):d
      #perform a GCD on complete rows of the matrix A
      if pMf[row, col] != 0
        # Modify by row addition, such that col,col is nonnegative
        if pMf[col,col] < 0
          pMf[col,:] = pMf[col,:] - trunc(pMf[col,col]/pMf[row,col])*pMf[row,:]
        end
        # get it positive
        if pMf[col,col] == 0
          pMf[col,:] = pMf[col,:] + sign(pMf[row,col])*pMf[row,:]
        end
        # get the second entry in that column also positive
        if pMf[row,col] < 0
          f = ceil(pMf[row,col]/pMf[col,col])
          if f==0
            f=f+1
          end
          pMf[row,:] = pMf[row,:] - f*pMf[col,:]
        end
        # Euclidean algorowthm on rows in order to get (row,col) to zero
        while pMf[row,col] != 0
          if abs(pMf[col,col]) > abs(pMf[row,col])
            f = floor(pMf[col,col]/pMf[row,col])
            if mod(pMf[col,col],pMf[row,col])==0
              f = f-sign(pMf[row,col])*sign(pMf[col,col])
            end
            pMf[col,:] = pMf[col,:]-f*pMf[row,:]
          else
            f = floor(pMf[row,col]/pMf[col,col])
            pMf[row,:] = pMf[row,:] - f*pMf[col,:]
          end
        end
      end
    end
    # get diagonal positive
    for row in 1:d
      if pMf[row,row]<0
        pMf[row,:] = -pMf[row,:]
      end
    end
    # get upper nonzero values of a colum to lie between 0 and
    for col in 2:d
      for row=(col-1):-1:1
        f=0
        if pMf[col,col] != 0
          if (pMf[row,col]<0) || (pMf[row,col]>=pMf[col,col])
            f = -floor(pMf[row,col]/pMf[col,col])
          end
        end
        pMf[row,:] = pMf[row,:] + f*pMf[col,:]
      end
    end
  end
  pMf
end
