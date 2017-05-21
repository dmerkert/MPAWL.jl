using MPAWL
using Base.Test

@testset "Pattern" begin
  MDiag = PatternMatrix([64,64])
  MRankOne = PatternMatrix([64 1;0 64])

  @test modM([65;63.1],MDiag) ≈ [1.0,63.1]
  @test modM([65;63.1],MDiag,target="symmetric") ≈ [1.0,-0.9]
  @test modM([65;63.1],MRankOne) ≈ [1.0,63.1]
  @test modM([65;63.1],MRankOne,target="symmetric") ≈ [0.0,-0.9]

  @test getm(MDiag) == 64^2
  @test getm(MRankOne) == 64^2

  @test patternSize(MDiag) == [64,64]
  @test patternSize(MRankOne) == [64^2]

  @test patternDimension(MDiag) == 2
  @test patternDimension(MRankOne) == 1

  @test patternBasis(MDiag) ≈ [0.015625 0; 0 0.015625]
  @test patternBasis(MRankOne) ≈ [0.000244140625,0.984375]

  @test generatingSetBasis(MDiag) == [1 0;0 1]
  @test generatingSetBasis(MRankOne) == [1,63]

end

# write your own tests here
