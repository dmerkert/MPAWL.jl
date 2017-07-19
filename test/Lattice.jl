using MPAWL
using Base.Test

@testset "Lattice" begin
  MDu = Lattice([64,64],target="unit")
  MDs = Lattice([64,64],target="symmetric")
  Lattice([64,8],target="unit")
  Lattice([64,8],target="symmetric")
  MRu = Lattice([64 1;0 64],target="unit")
  MRs = Lattice([64 1;0 64],target="symmetric")
  Lattice([64 4;0 64],target="unit")
  Lattice([64 4;0 64],target="symmetric")
  Lattice([64 4;1 32],target="unit")
  Lattice([64 4;1 32],target="symmetric")

  @test MDu.m == 64^2
  @test MDs.m == 64^2
  @test MRu.m == 64^2
  @test MRs.m == 64^2

  @test MDu.target == "unit"
  @test MDs.target == "symmetric"
  @test MRu.target == "unit"
  @test MRs.target == "symmetric"

  @test MDu.d == 2
  @test MDs.d == 2
  @test MRu.d == 2
  @test MRs.d == 2

  @test MDu.size == (64,64)
  @test MDs.size == (64,64)
  @test MRu.size == (64^2,)
  @test MRs.size == (64^2,)

  @test MDu.rank == 2
  @test MDs.rank == 2
  @test MRu.rank == 1
  @test MRs.rank == 1

  @test MDu.samplingLatticeBasis ≈ [0.015625 0; 0 0.015625]
  @test MDs.samplingLatticeBasis ≈ [0.015625 0; 0 0.015625]
  @test MRu.samplingLatticeBasis ≈ [0.000244140625,0.984375]
  @test MRs.samplingLatticeBasis ≈ [0.000244140625,-0.015625]

  @test MDu.frequencyLatticeBasis ≈ [1 0; 0 1]
  @test MDs.frequencyLatticeBasis ≈ [1 0; 0 1]
  @test MRu.frequencyLatticeBasis ≈ [1,64]
  @test MRs.frequencyLatticeBasis ≈ [1,0]

  @test MDu.patternNormalForm == MDu.M
  @test MDs.patternNormalForm == MDs.M
  @test MRu.patternNormalForm == MRu.M
  @test MRs.patternNormalForm == MRs.M

  #check if we have a biorthogonal basis
  @test MDu.samplingLatticeBasis.'*MDu.frequencyLatticeBasis ≈
  [1/64 0.0;0.0 1/64]
  @test MDs.samplingLatticeBasis.'*MDs.frequencyLatticeBasis ≈
  [1/64 0.0;0.0 1/64]
  @test mod.(MRu.samplingLatticeBasis.'*MRu.frequencyLatticeBasis,1.0) ≈
  [1.0/64^2]
  @test mod.(MRs.samplingLatticeBasis.'*MRs.frequencyLatticeBasis,1.0) ≈
  [1.0/64^2]
end
