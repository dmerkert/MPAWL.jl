export maxInd
"""
    Lattice

Lattice represents a regular integer matrix parametrizing a regular
anisotropic pattern on the torus. Has field `M` containing the matrix.


"""
immutable Lattice{I <: Integer}
  M :: Array{I,2}
  target :: String
  d :: I
  m :: I
  size :: Tuple
  dimension :: I
  samplingLatticeBasis :: Array{Float64,2}
  frequencyLatticeBasis :: Array{I,2}
  patternNormalForm :: Array{I,2}
  snfU :: Array{I,2}
  snfS :: Array{I,2}
  snfV :: Array{I,2}

  function Lattice{I}(M :: Array{I,2}; target="symmetric") where I <: Integer
    @argcheck size(M,1) == size(M,2)
    @argcheck det(M) != 0
    @argcheck target == "unit" || target == "symmetric"

    d = getd(M)
    m = getm(M)
    _patternSize = patternSize(M)
    _patternDimension = patternDimension(M)
    _patternBasis = patternBasis(M,target)
    _generatingSetBasis = generatingSetBasis(M.',target)
    _patternNormalForm = patternNormalForm(M)

    @assert d > 0
    @assert m > 0
    @assert _patternDimension <= d

    (snfU,snfS,snfV) = SNF(M)

    S = diag(snfS)
    #@show M
    for i in 1:size(_patternBasis,2)
      for j in 1:size(_patternBasis,2)
        y = _patternBasis[:,i]
        h = _generatingSetBasis[:,j]
        #@show dot(y,h)
        #@show 1.0/S[d-_patternDimension+i]

        i == j && @assert mod(dot(y,h),1.0) ≈ 1.0/S[d-_patternDimension+i]
        i != j && @assert ((mod(dot(y,h),1.0) ≈ 0.0) || (mod(dot(y,h),1.0) ≈ 1.0))
      end
    end

    new(
        M,
        target,
        d,
        m,
        _patternSize,
        _patternDimension,
        _patternBasis,
        _generatingSetBasis,
        _patternNormalForm,
        snfU,
        snfS,
        snfV
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
  convert(I,abs(det(M)))
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
function patternSize{I <: Integer}(M :: Array{I,2})
  d = diag(SNFWithoutTransform(M))
  s = tuple(d[d .> 1]...)
end

"""
    dM = patternDimension(M)
 Returns the number of elementary divisors of M, that are greater than 1.
 This corresponds to the number of basis vectors of the pattern and hence
 the dimension of the corresponding lattice.
"""
function patternDimension{I <: Integer}(M :: Array{I,2})
  sum(diag(SNFWithoutTransform(M)) .> 1)
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
 function patternBasis{I <: Integer}(M :: Array{I,2} ,target)
   @argcheck target == "unit" || target == "symmetric"

   d = getd(M)
   dM = patternDimension(M)
   (U,S,V) = SNF(M)
   V = V*(diagm(1./diag(S)))
   V = V[:,d-dM+1:d]
   for i in 1:dM
     V[:,i] = modM(V[:,i],eye(Int,d),target)
   end
   V
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
function generatingSetBasis{I <: Integer}(M :: Array{I,2},target)
  @argcheck target == "unit" || target == "symmetric"

  d = getd(M)
  dM = patternDimension(M)
  (U,S,V) = SNF(transpose(M))
  #(_U,S,_V) = SNF(M)
  #V = _U.'
  #U = _V.'

  V = transpose(inv(V))
  V = V[:,d-dM+1:d]
  for i = 1:dM
    V[:,i] = modM(V[:,i],M,target)
  end
  round.(Int,V)
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
          pMf[col,:] = pMf[col,:] - fix(pMf[col,col]/pMf[row,col])*pMf[row,:]
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
        if (pMf[row,col]<0) || (pMf[row,col]>=pMf[col,col])
          f = -floor(pMf[row,col]/pMf[col,col])
        end
        pMf[row,:] = pMf[row,:] + f*pMf[col,:]
      end
    end
  end
  pMf
end
