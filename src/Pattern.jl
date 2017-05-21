"""
    PatternMatrix

PatternMatrix represents a regular integer matrix parametrizing a regular
anisotropic pattern on the torus. Has field `M` containing the matrix.


"""
type PatternMatrix{I <: Integer}
  M :: Array{I,2}

  function PatternMatrix{I}(M :: Array{I,2})
    @argcheck size(M,1) == size(M,2)
    @argcheck det(M) != 0

    new(M)
  end
end
"""
    PatternMatrix(M): constructs a pattern matrix from the integer matrix `M`
"""
PatternMatrix{I}(M :: Array{I,2}) = PatternMatrix{I}(M)
"""
    PatternMatrix(v): constructs a diagonal pattern matrix from the integer
    vector `v`
"""
PatternMatrix{I}(v :: Array{I,1}) = PatternMatrix{I}(diagm(v))

"""
    m = getm(M): calculates the number of pattern points given by m = abs(det(M))
"""
function getm(M :: PatternMatrix)
  abs(det(M.M))
end

"""
    d = getd(M)

returns the dimension d = size(M.M,1) of the pattern
"""
function getd(M :: PatternMatrix)
  size(M.M,1)
end

"""
    modM!(k,M; target="unit")

Compute k mod M, where the `target` specifies the set of congruence class
representants, i.e. either "unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function modM!{R <: AbstractFloat}(
               k :: Array{R,1},
               M :: PatternMatrix;
               target="unit"
              )

  @argcheck length(k) == size(M.M,1)
  @argcheck target == "unit" || target == "symmetric"

  if target == "unit"
    k = M.M * mod(M.M\k,1);
  else
    k = M.M * (mod(M.M\k+0.5,1)-0.5);
  end
  k
end

modM!{I <: Integer}(
               k :: Array{I,1},
               M;
               target="unit"
              ) = round(I,modM!(convert(Float64(k),M,target=target)))

"""
    h = modM(k,M; target="unit")

Compute h = k mod M, where the `target` specifies the set of congruence class
representants, i.e. either "unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function modM(
     k,
     M; target="unit"
    )
  h = copy(k)
  modM!(h,M,target=target)
end

"""
    s = patternSize(M)
 Returns the elementary divisors of M, that are greater than 1.
 This corresponds to the size of the lattice.
"""
function patternSize(M :: PatternMatrix)
  d = diag(SNFWithoutTransform(M.M))
  s = d[d .> 1]
end

"""
    dM = patternDimension(M)
 Returns the number of elementary divisors of M, that are greater than 1.
 This corresponds to the number of basis vectors of the pattern and hence
 the dimension of the corresponding lattice.
"""
function patternDimension(M :: PatternMatrix)
  sum(diag(SNFWithoutTransform(M.M)) .> 1)
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
 function patternBasis(M :: PatternMatrix; target="unit")
   @argcheck target == "unit" || target == "symmetric"

   d = getd(M)
   dM = patternDimension(M)
   (U,S,V) = SNF(M.M)
   V = V*(diagm(1./diag(S)))
   V = V[:,d-dM+1:d]
   for i in 1:dM
     V[:,i] = modM(V[:,i],PatternMatrix(eye(Int,d)),target=target)
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
function generatingSetBasis(M :: PatternMatrix; target="unit")
  @argcheck target == "unit" || target == "symmetric"

d = getd(M)
dM = patternDimension(M)
(U,S,V) = SNF(transpose(M.M))
V = transpose(inv(V))
V = V[:,d-dM+1:d]
for i = 1:dM
  V[:,i] = modM(V[:,i],M,target=target)
end
round(Int,V)
end


"""
    generatingSet(M; target="unit") 

returns the set of p, where the `target` specifies the set of congruence class
representants, i.e. either "unit" for [0,1)^d or "symmetric" for
[-0.5,0.5)^d.oints belonging to a generating set of the PatternMatrix M.
"""
#= function generatingSet( =#

#=                       ) =#
