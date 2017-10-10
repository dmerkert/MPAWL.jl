export modM!,
modM,
frequencyLatticeBasisDecomp,
getSamplingIterator,
getFrequencyIterator,
getSamplingPoint,
getFrequencyPoint,
getFrequencyPoint!,
getUnitCell,
getMaxDualLatticeIndex,
isSublattice,
isLatticeDecomposition,
getPatternBasisDecomp,
getSuperpatternPoint,
getSuperpatternIndex

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

  #= @argcheck length(k) == getd(M) =#
  #= @argcheck target == "unit" || target == "symmetric" =#

  if target == "unit"
    k[:] = M * mod.(M\k,1.0);
  else
    k[:] = M * (mod.(M\k+0.5,1.0)-0.5);
  end
  k
end

function modM!{I <: Integer, MF}(
                             k :: Array{I,1},
                             M :: MF,
                             target
                            )

  #= @argcheck length(k) == getd(M) =1# =#
  #= @argcheck target == "unit" || target == "symmetric" =#

  if target == "unit"
    k[:] = round.(I,M * mod.(M\k,1.0))
  else
    k[:] = round.(I,M * (mod.(M\k+0.5,1.0)-0.5))
  end
  k
end

function modM!{I <: Integer, R, MF}(
                             k :: Array{I,1},
                             tmp :: Array{R,1},
                             M :: MF,
                             target
                            )
  modM!(k,M,target)

  #= @argcheck length(k) == getd(M) =1# =#
  #= @argcheck target == "unit" || target == "symmetric" =#

  #= tmp .= k =#
  #= A_ldiv_B!(M,tmp) =#
  #= if target == "unit" =#
  #=   tmp .= mod.(tmp,1.0) =#
  #=   A_mul_B!(tmp,M,tmp) =#
  #=   k .= round.(I,tmp) =#

  #=   #1= k[:] = round.(I,M * mod.(M\k,1.0)) =1# =#
  #= else =#
  #=   tmp .+= 0.5 =#
  #=   tmp .= mod.(tmp,1.0) =#
  #=   tmp .-= 0.5 =#
  #=   A_mul_B!(tmp,M,tmp) =#
  #=   k .= round.(I,tmp) =#
  #=   #1= k[:] = round.(I,M * (mod.(M\k+0.5,1.0)-0.5)) =1# =#
  #= end =#
  #= k =#
end

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
function frequencyLatticeBasisDecomp{I <: Integer,LI,MF,MF2}(
                                                k :: Array{I,1},
                                                L :: Lattice{LI,MF,MF2}
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
                          L :: Lattice{I,MF1,MF2},
                          coord :: CartesianIndex{N}
                         ) where {
                                  I,
                                  MF1,
                                  MF2,
                                  N
                                 }
  @argcheck length(coord.I) == L.rank
  (
   2.0*pi*modM(L.samplingLatticeBasis*(collect(coord.I)-1),
              eye(Int64,L.d),
              L.target
             )
  ) :: Array{Float64,1}
end

function getFrequencyPoint(
                           L :: Lattice{I,MF,MF2},
                          coord :: CartesianIndex{N}
                         ) where {N,I,MF,MF2}
  @argcheck length(coord.I) == L.rank
  modM(L.frequencyLatticeBasis*(collect(coord.I)-1),L.M',L.target)
end

function getFrequencyPoint!(
                          L :: Lattice,
                          coord :: CartesianIndex{N},
                          point :: Array{I,1},
                          tmp :: Array{R,1}
                         ) where {N,I,R}
  @argcheck length(coord.I) == L.rank

  A_mul_B!(point,L.frequencyLatticeBasis,[coord.I...]-1)
  #= modM!(point,tmp,L.MTFactorize,L.target) :: Array{I,1} =#
  modM!(point,L.M',L.target) :: Array{I,1}
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


"""
    isSublattice(LSuper,LSub)

returns true if LSub is a sublattice of LSuper, i.e. if there is a regular
integer matrix J such that LSuper.M = J*LSub.M
"""
function isSublattice(
                      LSuper :: Lattice{I,MF11,MF12},
                      LSub :: Lattice{I,MF21,MF22}
                     ) where {
                              I <: Integer,
                              MF11,MF12,
                              MF21,MF22
                             }

  J = round.(I,LSuper.M * inv(LSub.M))

  LSuper.M == J*LSub.M
end

function isLatticeDecomposition(
                      LSuper :: Lattice{I,MF11,MF12},
                      J :: Lattice{I,MF21,MF22},
                      LSub :: Lattice{I,MF31,MF32}
                     ) where {
                              I <: Integer,
                              MF11,MF12,
                              MF21,MF22,
                              MF31,MF32
                             }

  LSuper.M == J.M * LSub.M
end

function getPatternBasisDecomp(
                               L :: Lattice{I,MF11,MF12},
                               point :: Array{R,1}
                              ) where {I,MF11,MF12, R <: Real}

  @show coordTmp = round.(
                      I,
                      L.SNF[3]\
                      (L.SNF[2]*point/(2.0pi))
                     )
  @show coord = modM(
               coordTmp,
               L.M,
               "unit"
              )
  k = coordTmp
  @show mod.(L.M\k,1.0)
  @show L.M * mod.(L.M\k,1.0)
  @show round.(I,L.M * mod.(L.M\k,1.0))
  @show L.M

  @assert all(coord .>= 0)
  @assert all(coord .< L.size)

  coord
end

function getSuperpatternPoint(
                              coordDecomposition :: CartesianIndex{N},
                              pointSub :: Array{Float64,1},
                              LSub :: Lattice{I,MF11,MF12},
                              LDecomposition :: Lattice{I,MF31,MF32}
                             ) where {
                                      N,
                                      I,
                                      MF11,
                                      MF12,
                                      MF31,
                                      MF32
                                     }

  pointSub + LSub.M\getSamplingPoint(
                                     LDecomposition,
                                     coordDecomposition
                                    )
end

function getSuperpatternIndex(
                              coordDecomposition :: CartesianIndex{N},
                              coordSub :: CartesianIndex{N2},
                              LSub :: Lattice{I,MF11,MF12},
                              LSuper :: Lattice{I,MF21,MF22},
                              LDecomposition :: Lattice{I,MF31,MF32}
                             ) where {
                                      I,
                                      MF11,
                                      MF12,
                                      MF21,
                                      MF22,
                                      MF31,
                                      MF32,
                                      N,N2
                                     }

  if isdiag(LDecomposition.M)
    return coordSub + coordDecomposition
  else

    pointSub = getSamplingPoint(
                                LSub,
                                coordSub
                               )
    superpatternPoint = getSuperpatternPoint(
                                             coordDecomposition,
                                             pointSub,
                                             LSub,
                                             LDecomposition
                                            )

    coordSuper = CartesianIndex(((getPatternBasisDecomp(LSuper,superpatternPoint)+1)...))

    @assert all(coordSuper.I .>= 1)
    @assert all(coordSuper.I .<= LSuper.size)

    return coordSuper
  end
end
