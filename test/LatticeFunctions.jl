using MPAWL
using Base.Test

@testset "Lattice Functions" begin
  MDu = Lattice([64,64],target="unit")
  MDs = Lattice([64,64],target="symmetric")
  MRu = Lattice([64 1;0 64],target="unit")
  MRs = Lattice([64 1;0 64],target="symmetric")

  @test modM([65;63.1],MDu) ≈ [1.0,63.1]
  @test modM([65;63.1],MDs) ≈ [1.0,-0.9]
  @test modM([65;63.1],MRu) ≈ [1.0,63.1]
  @test modM([65;63.1],MRs) ≈ [0.0,-0.9]

  k = [66;5]
  v = frequencyLatticeBasisDecomp(k,MDu)
  @test modM(k,MDu.M',MDu.target) ==
  modM(MDu.frequencyLatticeBasis*v,MDu.M',MDu.target)

  v = frequencyLatticeBasisDecomp(k,MDs)
  @test modM(k,MDs.M',MDs.target) ==
  modM(MDu.frequencyLatticeBasis*v,MDs.M',MDs.target)

  v = frequencyLatticeBasisDecomp(k,MRu)
  @test modM(k,MRu.M',MRu.target) ==
  modM(MRu.frequencyLatticeBasis*v,MRu.M',MRu.target)

  v = frequencyLatticeBasisDecomp(k,MRs)
  @test modM(k,MRs.M',MRs.target) ==
  modM(MRu.frequencyLatticeBasis*v,MRs.M',MRs.target)


  for i in getSamplingIterator(MDu)
    @test all(0 .<= getSamplingPoint(MDu,i) .< 2.0pi)
    @test round.(getFrequencyPoint(MDu,i)) == getFrequencyPoint(MDu,i)
    @test all(-pi .<= getSamplingPoint(MDs,i) .< pi)
    @test round.(getFrequencyPoint(MDs,i)) == getFrequencyPoint(MDs,i)
  end

  for i in getFrequencyIterator(MRu)
    i = CartesianIndex((1))
    @test all(0 .<= getSamplingPoint(MRu,i) .< 2.0pi)
    @test round.(getFrequencyPoint(MRu,i)) == getFrequencyPoint(MRu,i)
    @test all(-pi .<= getSamplingPoint(MRs,i) .< pi)
    @test round.(getFrequencyPoint(MRs,i)) == getFrequencyPoint(MRs,i)
  end

  LUnit = Lattice([64 1;5 32],target="unit")
  LSym = Lattice([64 1;5 32],target="symmetric")

  @test getMaxDualLatticeIndex(LUnit) == [70;34]
  @test getMaxDualLatticeIndex(LSym) == [36;18]
  @test getMaxDualLatticeIndex(LUnit) == [70;34]
  @test getMaxDualLatticeIndex(LUnit,cubeSize=4.0) == [277;133]
  @test getMaxDualLatticeIndex(LSym,cubeSize=4.0) == [139;67]

  @test getUnitCell(Lattice([64 0;0 64])) ≈ [-1/32 -1/32;
                                              -1/32 1/32;
                                              1/32  1/32;
                                              1/32 -1/32]

  @test getUnitCell(Lattice([128 0;0 64])) ≈ [-1/64 -1/32;
                                              -1/64 1/32;
                                               1/64  1/32;
                                               1/64 -1/32]

  @test getMaxDualLatticeIndex(LUnit,[1.4;77.2]) == [476;2472]
  @test getMaxDualLatticeIndex(LSym,[1.4;77.2]) == [238;1236]
end
