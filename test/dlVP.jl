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

  ckphim =
  ifftshift(
            [ 0 0 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0
             0 0 0.09375 0.2500 0.40625 0.5000 0.40625 0.2500 0.09375 0 0
             0 0 0.1875 0.5000 0.8125 1.0000 0.8125 0.5000 0.1875 0 0
             0 0 0.09375 0.2500 0.40625 0.5000 0.40625 0.2500 0.09375 0 0
             0 0 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0
            ]
           )

  ckphimO =
  ifftshift(
            [
             0 0 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0
             0 0 0.158999682000954 0.500000000000000 0.688998622004134 0.707106781186547 0.688998622004134 0.500000000000000 0.158999682000954 0 0
             0 0 0.224859506698758 0.707106781186547 0.974391195694620 1.000000000000000 0.974391195694620 0.707106781186547 0.224859506698758 0 0
             0 0 0.158999682000954 0.500000000000000 0.688998622004134 0.707106781186547 0.688998622004134 0.500000000000000 0.158999682000954 0 0
             0 0 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0
            ]
           )


  N = size(ckphim)
  for coord in getFrequencyIterator(N)
    frequency = getFrequency(N,coord)

    phi = delaValleePoussinMean(
                                      frequency,
                                      L,
                                      g,
                                      false
                                     )
    phiO = delaValleePoussinMean(
                                        frequency,
                                        L,
                                        g,
                                        true
                                       )

    @test phi ≈ ckphim[coord]
    @test phiO ≈ ckphimO[coord]
  end
end
