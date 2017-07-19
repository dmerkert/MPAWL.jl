using MPAWL
using Base.Test

@testset "FFT" begin
  MDu = Lattice([64,64],target="unit")
  MRs = Lattice([64 1;0 64],target="symmetric")

  A = rand((64,64,3))
  AFrequency = FFT(A,MDu,[1;2])
  @test real(AFrequency[1,1,:])[:] ≈ sum(A,(1,2))[:]
  setFourierCoefficient!(AFrequency,MDu,[1.0;2.0;3.0],[1;2])
  A = IFFT(AFrequency,MDu,[1;2])
  @test sum(A,(1,2))[:] ≈ [1.0;2.0;3.0]

  A = rand((64^2,2))
  AFrequency = FFT(A,MRs,[1])
  @test real(AFrequency[1,:])[:] ≈ sum(A,(1))[:]
  setFourierCoefficient!(AFrequency,MRs,[3.1;4.1],[1])
  A = IFFT(AFrequency,MRs,[1])
  @test sum(A,(1))[:] ≈ [3.1;4.1]

  A = rand((64,3,64))
  AFrequency = FFT(A,MDu,[1;3])
  setFourierCoefficient!(AFrequency,MDu,[1.0;2.0;3.0],[1;3],[64;5])
  @test AFrequency[64,:,5] == [1.0;2.0;3.0]

  A = rand(Complex128,(64,64,3))
  AFrequency = FFT(A,MDu,[1;2])
  @test AFrequency[1,1,:][:] ≈ sum(A,(1,2))[:]
  setFourierCoefficient!(AFrequency,MDu,[1.0;2.0;3.0],[1;2])
  A = IFFT(AFrequency,MDu,[1;2])
  @test sum(A,(1,2))[:] ≈ [1.0;2.0;3.0]
end
