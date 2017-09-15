using MPAWL
using Base.Test
using ProgressMeter

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

function getMatrix(
                   k1,k2,k3, a
                  )

  b = 32

  if a == 0
    return [b k1 k2;
            0 b k3;
            0 0 b]
  elseif a == 1
    return [b 0 k2;
            k1 b k3;
            0 0 b]
  elseif a == 2
    return [b k1 0;
            0 b k3;
            k2 0 b]
  elseif a == 3
    return [b k1 k2;
            0 b 0;
            0 k3 b]
  elseif a == 4
    return [b 0 0;
            k1 b k3;
            k2 0 b]
  elseif a == 5
    return [b 0 k2;
            k1 b 0;
            0 k3 b]
  elseif a == 6
    return [b k1 0;
            0 b 0;
            k2 k3 b]
  elseif a == 7
    return [b 0 0;
            k1 b 0;
            k2 k3 b]
  end
end

@testset "Lattice generation" begin
  aRange = 0:7
  k1Range = -32:1:32

  simulations = length(aRange)*length(k1Range)^3
  prog = Progress(simulations,dt=1.0, barglyphs=BarGlyphs("[=> ]"), barlen=50)


  for a in aRange
    for k1 in k1Range
      for k2 in k1Range
        for k3 in k1Range
          M = getMatrix(k1,k2,k3,a)
          next!(prog)
          if !(det(M) ≈ 0)
            L = Lattice(M)
          end
        end
      end
    end
  end
end
