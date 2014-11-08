 !****************************************************************  
 ! MD_BoundaryCondition.f90
 ! 
 ! Created on: Jul 29, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_BoundaryCondition
    USE MD_Definition
    USE MD_Parameter
    USE MD_Quantity
    IMPLICIT NONE

CONTAINS

    real(kind=double) function GetVaccumSolution()
        implicit none
        GetVaccumSolution = 1.0_rp
    end function GetVaccumSolution

    real(kind=double) function CalculatePotential(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: radius
        radius = sqrt((xcord-xc(0))**2 + (ycord-xc(1))**2 + (zcord-xc(2))**2)
        CalculatePotential = -G*M/radius
    end function CalculatePotential

    real(kind=double) function GetR2Solution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: r2
        r2 = (xcord-xc(0))**2 + (ycord-xc(1))**2 + (zcord-xc(2))**2
        GetR2Solution = -2.0_rp*pi*G*(M/((4.0_rp/3.0_rp)*pi*R**3))*(R**2 - (1.0_rp/3.0_rp)*r2)
    end function GetR2Solution

    real(kind=double) function GetSinSolution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        GetSinSolution = -(4.0_rp*pi*G)*((2.0_rp*R/pi)**2)*(pi/(48.0_rp*R**3))* &
            (cos(pi*xcord/(2.0_rp*R)) + cos(pi*ycord/(2.0_rp*R)) + cos(pi*zcord/(2.0_rp*R)))
    end function GetSinSolution

    real(kind=double) function GetPoly6Solution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        GetPoly6Solution = xcord**6+ycord**6+zcord**6
    end function GetPoly6Solution

    real(kind=double) function GetGaussianSolution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: d, alpha
        d = sqrt(xcord**2 + ycord**2 + zcord**2)
        alpha = 64.0_rp
        GetGaussianSolution = -4.0_rp*pi*G*sqrt(pi)*erf(d*sqrt(alpha))/(4.0_rp*d*alpha**1.5)
    end function GetGaussianSolution

    real(kind=double) function GetBoundaryValue(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord

        if(scenerio .eq. 0) then
            GetBoundaryValue = GetVaccumSolution()
        else if(scenerio .eq. 1 .or. scenerio .eq. 2) then
            GetBoundaryValue = CalculatePotential(xcord,ycord,zcord)
        else if(scenerio .eq. 3) then
            GetBoundaryValue = GetR2Solution(xcord,ycord,zcord)
        else if(scenerio .eq. 4) then
            GetBoundaryValue =  GetSinSolution(xcord,ycord,zcord)
        else if(scenerio .eq. 5) then
            GetBoundaryValue = GetPoly6Solution(xcord,ycord,zcord)
        else if(scenerio .eq. 6) then
            GetBoundaryValue = GetGaussianSolution(xcord,ycord,zcord)
        end if

    end function GetBoundaryValue


END MODULE MD_BoundaryCondition
