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

!    real(kind=double) :: xcord(0:2)

    print *, ''
    print *, '***************************'
    print *, 'Boundary Value Test'
    print *, '***************************'
    print *, ''

!    call GetCordinate(2,7,3,xcord)
    print *, 'h=',h
    print *, 'x=',x(2),y(7),z(3)

    print *, 'VaccumSolution=',GetVaccumSolution()

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M=1.0_rp

    print *, 'Potential=',CalculatePotential(x(2),y(7),z(3))

    xc(0) = 0.1_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M=1.0_rp

    print *, 'Potential=',CalculatePotential(x(2),y(7),z(3))

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M=1.0_rp

    print *, 'R2Solution=',GetR2Solution(x(2),y(7),z(3))

 end subroutine Test_BoundaryValue
