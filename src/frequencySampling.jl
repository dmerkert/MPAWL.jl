export getFrequencyIterator,
       getFrequency

function getFrequencyIterator(N :: Tuple)
 CartesianRange(N)
end

function getFrequency(N :: Tuple, coord :: CartesianIndex)
  modM(collect(coord.I)-1,diagm(collect(N)),"symmetric")
end
