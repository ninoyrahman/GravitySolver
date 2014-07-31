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
    use MD_Quantity
    implicit none

    real(kind=double) :: x(0:2)

    print *, ''
    print *, '***************************'
    print *, 'Boundary Value Test'
    print *, '***************************'
    print *, ''

    call GetCordinate(2,7,3,x)
    print *, 'h=',h
    print *, 'x=',x(0),x(1),x(2)

    print *, 'VaccumSolution=',GetVaccumSolution()

    xc(0) = 0.0
    xc(1) = 0.0
    xc(2) = 0.0
    M=1.0

    print *, 'Potential=',CalculatePotential(x)

    xc(0) = 0.1
    xc(1) = 0.0
    xc(2) = 0.0
    M=1.0

    print *, 'Potential=',CalculatePotential(x)

    xc(0) = 0.0
    xc(1) = 0.0
    xc(2) = 0.0
    M=1.0

    print *, 'R2Solution=',GetR2Solution(x)

 end subroutine Test_BoundaryValue
