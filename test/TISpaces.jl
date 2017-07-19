using MPAWL
using Base.Test

@testset "TI Spaces" begin
  L = Lattice([8 1;0 16])
  data = rand(Complex128,(64,128))
  ckPhi = rand(Complex128,(64,128))

  bSq = bracketSums(L,data)

  data2 = copy(data)
  coeffsSpace2Fourier!(data2,L,bSq,ckPhi)

  dataC = rand(Complex128,128)
  dataR = rand(Float64,128)

  changeBasis!(dataC,dataC,L,bSq,[1],inputDomain="Space",outputDomain="Space")
  changeBasis!(dataC,dataC,L,bSq,[1],inputDomain="Space",outputDomain="Fourier")
  changeBasis!(dataC,dataC,L,bSq,[1],inputDomain="Fourier",outputDomain="Space")
  changeBasis!(dataC,dataC,L,bSq,[1],inputDomain="Fourier",outputDomain="Fourier")

  changeBasis!(dataR,dataC,L,bSq,[1],inputDomain="Space",outputDomain="Space")
  changeBasis!(dataC,dataR,L,bSq,[1],inputDomain="Space",outputDomain="Space")

  changeBasis!(dataR,dataC,L,bSq,[1],inputDomain="Space",outputDomain="Fourier")
  changeBasis!(dataC,dataR,L,bSq,[1],inputDomain="Fourier",outputDomain="Space")
end
