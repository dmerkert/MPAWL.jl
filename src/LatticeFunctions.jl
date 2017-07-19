export modM!,
modM,
frequencyLatticeBasisDecomp,
getSamplingIterator,
getFrequencyIterator,
getSamplingPoint,
getFrequencyPoint,
getUnitCell,
getMaxDualLatticeIndex

"""
    _modM!(k,M; target="unit")

Compute k mod M, where the `target` specifies the set of congruence class
representants, i.e. either "unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function modM!{R <: AbstractFloat, I <: Integer}(
               k :: Array{R,1},
               M :: Array{I,2},
               target
              )

  @argcheck length(k) == getd(M)
  @argcheck target == "unit" || target == "symmetric"

  if target == "unit"
    k = M * mod.(M\k,1.0);
  else
    k = M * (mod.(M\k+0.5,1.0)-0.5);
  end
  k
end

modM!{I <: Integer}(
               k :: Array{I,1},
               M,
               target
              ) = round.(I,modM!(convert(Array{Float64,1},k),M,target))

"""
    h = _modM(k,M; target="unit")

Compute h = k mod M, where the `target` specifies the set of congruence class
representants, i.e. either "unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function modM(
     k,
     M,
     target
    )
  h = copy(k)
  modM!(h,M,target)
end




"""
    modM!(k,M)

Compute k mod M.
"""
modM!(k, M :: Lattice) = modM!(k,M.M,M.target) 
"""
    h = modM(k,M)

Compute h = k mod M.
"""
function modM(k, M)
  h = copy(k)
  modM!(h,M)
end


"""
    frequencyLatticeBasisDecomp(k,M; target="unit")
Decompose the vector k with respect to the frequency lattice basis of pattern matrix
M, i.e. k = modM(V.v,M), where V = frequencyLatticeBasis(M)

The `target` specifies the set of congruence class representants, i.e. either
"unit" for [0,1)^d or "symmetric" for [-0.5,0.5)^d.
"""
function frequencyLatticeBasisDecomp{I <: Integer}(
                                                k :: Array{I,1},
                                                L :: Lattice
                                               )
  @argcheck length(k) == L.d

  d = L.d
  dM = L.rank
  # decomposing w.r.t. M^T we use the U and S of M, following the thesis
  # from above: aBV is V^-T and thats used as inverse -> transpose and
  # multiply
  v = modM(L.SNF[3].'*k,L.SNF[2],L.target)
  v = v[(d-dM+1):d]
  v
end

function getSamplingIterator(L :: Lattice)
  CartesianRange(L.size)
end

function getFrequencyIterator(L :: Lattice)
  CartesianRange(L.size)
end

function getSamplingPoint(
                          L :: Lattice,
                          coord :: CartesianIndex
                         )
  @argcheck length(coord.I) == L.rank
  2.0*pi*modM(L.samplingLatticeBasis*(collect(coord.I)-1),
              eye(Int64,L.d),
              L.target
             )
end

function getFrequencyPoint(
                          L :: Lattice,
                          coord :: CartesianIndex
                         )
  @argcheck length(coord.I) == L.rank
  modM(L.frequencyLatticeBasis*(collect(coord.I)-1),L.M',L.target)
end

function getUnitCell(L :: Lattice)
  [
   (0.5*L.M\[-1,-1])';
   (0.5*L.M\[-1, 1])';
   (0.5*L.M\[ 1, 1])';
   (0.5*L.M\[ 1,-1])'
  ]
end

"""
    getMaxDualLatticeIndex(M)

returns the maximal index (vector) an entry of the (symmetric)
generating Set uses when looking at the generating set of `M`.

"""
function getMaxDualLatticeIndex(L :: Lattice,
                                cubeSize :: Array{Float64,1}
                               )

  @argcheck length(cubeSize) == L.d

  shift = 0.5cubeSize
  if L.target == "unit"
    shift = cubeSize
  end

  maxInd = zeros(Int,L.d)
  for coord in CartesianRange(
                              CartesianIndex((zeros(Int,L.d)...)),
                              CartesianIndex((ones(Int,L.d)...))
                             )

    maxInd = max.(maxInd,
                  abs.(ceil.(Int,L.M'*(collect(coord.I) - shift)))+1
                 )
  end
  maxInd
end

getMaxDualLatticeIndex(L :: Lattice;
                       cubeSize :: Float64 = 1.0
                      ) = getMaxDualLatticeIndex(L,cubeSize*ones(L.d))

