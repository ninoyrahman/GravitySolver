 !****************************************************************  
 ! Test_BoundaryValue.f90
 ! 
 ! Created on: Jul 26, 2014
 ! Author: ninoy
 !***************************************************************** 

 subroutine Test_BoundaryValue
    use MD_Parameter
    use MD_Helper
    use MD_BoundaryCondition
    implicit none

    real(kind=double) :: xcord,ycord,zcord,r2,d

    print *, ''
    print *, '***************************'
    print *, 'Boundary Value Test'
    print *, '***************************'
    print *, ''

    xcord = 0.1
    ycord = -0.25
    zcord = 1.75
    M=1.0_rp
    R = 0.5_rp
    r2 = xcord**2 + ycord**2 + zcord**2
    d = sqrt(r2)

    print *, 'h=',h
    print *, 'x=',xcord,ycord,zcord
    print *, 'M=',M
    print *, 'R=',R
    print *, 'r2=',r2
    print *, 'd=',d

    print *, 'VaccumSolution=',GetVaccumSolution()

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp

    print *, 'CenteredSphereSolution=',CalculatePotential(xcord,ycord,zcord)

    xc(0) = 0.1_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp

    print *, 'OffCenteredSphereSolution=',CalculatePotential(xcord,ycord,zcord)

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp

    print *, 'CenteredSolution=',GetCenteredSolution(xcord,ycord,zcord)

    print *, 'SinSolution=',GetSinSolution(xcord,ycord,zcord)

    print *, 'Poly6Solution=',GetPoly6Solution(xcord,ycord,zcord)

    print *, 'GaussianSolution=',GetGaussianSolution(xcord,ycord,zcord)

 end subroutine Test_BoundaryValue
