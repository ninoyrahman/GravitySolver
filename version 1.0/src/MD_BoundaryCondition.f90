 !****************************************************************  
 ! MD_BoundaryCondition.f90
 ! 
 ! Created on: Jul 29, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_BoundaryCondition
    USE MD_Definition
    USE MD_Parameter
    USE MD_CordinateTransform
    IMPLICIT NONE

CONTAINS

    real(kind=double) function GetVaccumSolution()
        implicit none
        GetVaccumSolution = 1.0
    end function GetVaccumSolution

    real(kind=double) function CalculatePotential(x)
        implicit none
        real(kind=double), dimension(0:2), intent(in) :: x
        real(kind=double) :: r
        r = sqrt((x(0)-xc(0))**2 + (x(1)-xc(1))**2 + (x(2)-xc(2))**2)
        CalculatePotential = -G*M/r
    end function CalculatePotential

    real(kind=double) function GetR2Solution(x)
        implicit none
        real(kind=double), dimension(0:2), intent(in) :: x
        GetR2Solution = (x(0)-xc(0))**2 + (x(1)-xc(1))**2 + (x(2)-xc(2))**2
    end function GetR2Solution

    real(kind=double) function GetSinSolution(x)
        implicit none
        real(kind=double), dimension(0:2), intent(in) :: x
        GetSinSolution = sin(x(0))+sin(x(1))+sin(x(2))
    end function GetSinSolution

    real(kind=double) function GetPoly6Solution(x)
        implicit none
        real(kind=double), dimension(0:2), intent(in) :: x
        GetPoly6Solution = x(0)**6+x(1)**6+x(2)**6
    end function GetPoly6Solution

    real(kind=double) function GetBoundaryValue(x)
        implicit none
        real(kind=double), dimension(0:2), intent(in) :: x

        if(scenerio .eq. 0) then
            GetBoundaryValue = GetVaccumSolution()
        else if(scenerio .eq. 1 .or. scenerio .eq. 2) then
            GetBoundaryValue = CalculatePotential(x)
        else if(scenerio .eq. 3) then
            GetBoundaryValue = GetR2Solution(x)
        else if(scenerio .eq. 4) then
            GetBoundaryValue = GetSinSolution(x)
        else if(scenerio .eq. 5) then
            GetBoundaryValue = GetPoly6Solution(x)
        end if

    end function GetBoundaryValue


END MODULE MD_BoundaryCondition
