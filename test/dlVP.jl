using MPAWL
using Base.Test

@testset "pyramidFunction" begin
  alpha = [0.25;0.4]

  @test pyramidFunction(0.3,0.0) ≈ 1.0
  @test pyramidFunction(0.3,0.1) ≈ 1.0
  @test pyramidFunction(0.3,0.2) ≈ 1.0
  @test pyramidFunction(0.3,0.3) ≈ 0.83333333333333333
  @test pyramidFunction(0.3,0.4) ≈ 0.66666666666666666
  @test pyramidFunction(0.3,0.5) ≈ 0.5
  @test pyramidFunction(0.0,-0.5) ≈ 0.5
  @test pyramidFunction(0.0,-0.2) ≈ 1.0
  @test pyramidFunction(0.0,0.8) ≈ 0.0

  @test pyramidFunction(alpha, [0.1;0.0]) ≈ 1.0
  @test pyramidFunction(alpha, [0.1;0.1]) ≈ 1.0
  @test pyramidFunction(alpha, [0.1;0.2]) ≈ 0.8750
  @test pyramidFunction(alpha, [0.1;0.3]) ≈ 0.7500
  @test pyramidFunction(alpha, [0.1;0.4]) ≈ 0.6250
  @test pyramidFunction(alpha, [0.1;0.5]) ≈ 0.5
  @test pyramidFunction(alpha, [0.1;0.6]) ≈ 0.375
  @test pyramidFunction(alpha, [0.1;0.7]) ≈ 0.25
  @test pyramidFunction(alpha, [0.1;0.8]) ≈ 0.125
  @test pyramidFunction(alpha, [0.1;0.9]) ≈ 0.0



  @test pyramidFunction([0.2523423;0.42349],[0.123421;0.564534]) ≈ 0.423806937590026

end

@testset "dlVP" begin
  L = Lattice([2 0;0 4])
  g = [0.25;0.4]

  

  (ckphi,ckBSq) = delaValleePoussinMean(L,g,orthonormalize=false)
  (ckphiO,ckBSqO) = delaValleePoussinMean(L,g,orthonormalize=true)
  #= println("===================") =#
  #= (ckphiO,ckBSqO) = delaValleePoussinMean(L,g,orthonormalize=true) =#
  ckphim =
  [ 0 0 0 0 0 0 0 0 0 0 0
   0 0 0 0 0 0 0 0 0 0 0
   0 0 0.09375 0.2500 0.40625 0.5000 0.40625 0.2500 0.09375 0 0
   0 0 0.1875 0.5000 0.8125 1.0000 0.8125 0.5000 0.1875 0 0
   0 0 0.09375 0.2500 0.40625 0.5000 0.40625 0.2500 0.09375 0 0
   0 0 0 0 0 0 0 0 0 0 0
   0 0 0 0 0 0 0 0 0 0 0
  ]

  ckphimO =
 [
 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0
 0 0 0.158999682000954 0.500000000000000 0.688998622004134 0.707106781186547 0.688998622004134 0.500000000000000 0.158999682000954 0 0
 0 0 0.224859506698758 0.707106781186547 0.974391195694620 1.000000000000000 0.974391195694620 0.707106781186547 0.224859506698758 0 0
 0 0 0.158999682000954 0.500000000000000 0.688998622004134 0.707106781186547 0.688998622004134 0.500000000000000 0.158999682000954 0 0
 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0
]


@test ckphim ≈ real(fftshift(ckphi))
@test ckphimO ≈ real(fftshift(ckphiO))
@test norm(ckBSqO[:] - 1.0) ≈ 0.0 atol=1e-15


N = size(ckphi)
LUnit = Lattice(L.M,target="unit") :: Lattice{Int}
for coord in getFrequencyIterator(N)
  frequency = getFrequency(N,coord)
  h = frequencyLatticeBasisDecomp(frequency,LUnit)
  phiOrig = ckphi[coord]
  phiOrigO = ckphiO[coord]
  bSqOrig = ckBSq[(h+1)...]
  bSqOrigO = ckBSqO[(h+1)...]

  (phi,BSq) = delaValleePoussinMean(
                                    frequency,
                                    L,
                                    g,
                                    false
                                   )
  (phiO,BSqO) = delaValleePoussinMean(
                                      frequency,
                                    L,
                                    g,
                                    true
                                   )

  @test phi ≈ phiOrig
  @test phiO ≈ phiOrigO
  @test BSq ≈ bSqOrig
  @test BSqO ≈ bSqOrigO
end

#= frequency = CartesianIndex((1,0)) =#
#= for i in BracketSumIterator( =#
#=                             frequency, =#
#=                             CartesianRange( =#
#=                                            CartesianIndex(ntuple(i -> -1,2)), =#
#=                                            CartesianIndex(ntuple(i -> 1,2)) =#
#=                                           ), =#
#=                             L =#
#=                            ) =#
#=   @test frequencyLatticeBasisDecomp([frequency.I...],LUnit) == [frequency.I...] =#
#= end =#
end
