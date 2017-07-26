export BracketSumIterator,
BracketSumIteratorState,
start,
done,
next,
iteratorsize,
size,
length,
iteratoreltype

type BracketSumIterator{N, I, MF, MF2}
  h :: Array{I,1}  # h+M^t z
  support :: CartesianRange{CartesianIndex{N}} # für z, e.g. für dlvp: [-1,1] oder so
  L :: Lattice{I,MF,MF2}
end

BracketSumIterator(
                   h :: Array{I,1},
                   supportLower :: Array{I,1},
                   supportUpper :: Array{I,1},
                   L :: Lattice{I,MF,MF2}
                  ) where {I,MF,MF2} =
BracketSumIterator(
                   h,
                   CartesianRange(
                                  CartesianIndex(Tuple(supportLower)),
                                  CartesianIndex(Tuple(supportUpper))
                                 ),
                   L
                  )

BracketSumIterator(
                   h :: Array{I,1},
                   supportLower :: I,
                   supportUpper :: I,
                   L :: Lattice{I,MF,MF2}
                  ) where {I,MF,MF2} =
BracketSumIterator(h,
                   supportLower*ones(I,L.d),
                   supportUpper*ones(I,L.d),
                   L
                  )


type BracketSumIteratorState{N}
  cartesianState :: CartesianIndex{N}
end

function Base.start(iter :: BracketSumIterator{N,I,MF,MF2}) where {N,I,MF,MF2}
  BracketSumIteratorState(start(iter.support))
end

function Base.done(
              iter :: BracketSumIterator{N,I,MF,MF2},
              state :: BracketSumIteratorState{N}
             ) where {N,I,MF,MF2}
  done(iter.support,state.cartesianState)
end

function Base.next(
              iter :: BracketSumIterator{N,I,MF,MF2},
              state :: BracketSumIteratorState{N}
             ) where {N,I,MF,MF2}
  (z,cartesianState) = next(iter.support,state.cartesianState)
  state.cartesianState = cartesianState

  (iter.h + iter.L.MTFactorize*[z.I...], state)
end

Base.iteratorsize(iter :: BracketSumIterator{N,I,MF,MF2}) where {N,I,MF,MF2} =
iteratorsize(iter.support)

Base.size(iter :: BracketSumIterator{N,I,MF,MF2}) where {N,I,MF,MF2} = size(iter.support)
Base.length(iter :: BracketSumIterator{N,I,MF,MF2}) where {N,I,MF,MF2} = length(iter.support)
Base.iteratoreltype(iter :: BracketSumIterator{N,I,MF,MF2}) where {N,I,MF,MF2} =
iteratoreltype(iter.support)
