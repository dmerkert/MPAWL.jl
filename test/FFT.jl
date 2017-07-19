using MPAWL
using Base.Test

@testset "FFT" begin
  MDu = Lattice([64,64],target="unit")
  MRs = Lattice([64 1;0 64],target="symmetric")

  A = rand((64,64,3))
  AFrequency = FFT(A,MDu,FirstDimensionsFFT)
  @test real(AFrequency[1,1,:])[:] ≈ sum(A,(1,2))[:]
  setZerothFourierCoefficient!(AFrequency,MDu,[1.0;2.0;3.0],FirstDimensionsFFT)
  A = IFFT(AFrequency,MDu,FirstDimensionsFFT)
  @test sum(A,(1,2))[:] ≈ [1.0;2.0;3.0]

  A = rand((64^2,2))
  AFrequency = FFT(A,MRs,FirstDimensionsFFT)
  @test real(AFrequency[1,:])[:] ≈ sum(A,(1))[:]
  setZerothFourierCoefficient!(AFrequency,MRs,[3.1;4.1],FirstDimensionsFFT)
  A = IFFT(AFrequency,MRs,FirstDimensionsFFT)
  @test sum(A,(1))[:] ≈ [3.1;4.1]
end
