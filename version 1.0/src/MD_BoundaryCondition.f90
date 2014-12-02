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

    real(kind=double) function GetOCGaussianSolution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: left = 1.0_rp
        real(kind=double) :: right = 1.0_rp
        real(kind=double) :: alpha = 64.0_rp
        real(kind=double) :: rl, rr
        rl = sqrt((xcord + 0.5_rp)**2 + (ycord - 0.0_rp)**2 + (zcord - 0.0_rp)**2)
        rr = sqrt((xcord - 0.5_rp)**2 + (ycord - 0.0_rp)**2 + (zcord - 0.0_rp)**2)
        GetOCGaussianSolution = left*(-4.0_rp*pi*G)*sqrt(pi)*erf(rl*sqrt(alpha))/(4.0_rp*rl*alpha**1.5) + &
            right*(-4.0_rp*pi*G)*sqrt(pi)*erf(rr*sqrt(alpha))/(4.0_rp*rr*alpha**1.5)
    end function GetOCGaussianSolution

    real(kind=double) function GetR2Solution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: r2,unirho
        unirho = 1.0_rp/((4.0_rp/3.0_rp)*pi*R**3)
        r2 = (xcord-xc(0))**2 + (ycord-xc(1))**2 + (zcord-xc(2))**2
        if(r2 .lt. R*R) GetR2Solution = -2.0_rp*pi*G*unirho*(R**2 - (1.0_rp/3.0_rp)*r2)
        if(r2 .ge. R*R) GetR2Solution = -(4.0_rp/3.0_rp)*pi*G*unirho*(R**3/sqrt(r2))
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
        d = sqrt((xcord-xc(0))**2 + (ycord-xc(1))**2 + (zcord-xc(2))**2)
        alpha = 64.0_rp
        GetGaussianSolution = -4.0_rp*pi*G*sqrt(pi)*erf(d*sqrt(alpha))/(4.0_rp*d*alpha**1.5)
    end function GetGaussianSolution

    real(kind=double) function GetPloytrop1Solution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: alpha, beta, d
        alpha = 0.2_rp
        d = sqrt(xcord**2 + ycord**2 + zcord**2)
        beta = sin(d/alpha)/(d/alpha)
        if(d .eq. 0.0_rp) beta = 1.0
        if(d .le. R) GetPloytrop1Solution = -4.0_rp*pi*G*(alpha**2)*(1 + beta)
        if(d .gt. R) GetPloytrop1Solution = -4.0_rp*pi*G*(alpha**2)*(pi/(d/alpha))
    end function GetPloytrop1Solution

    real(kind=double) function GetPloytrop5Solution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: rc, r2
        rc = 0.2_rp
        r2 = xcord**2 + ycord**2 + zcord**2
        GetPloytrop5Solution = -(4.0_rp/3.0_rp)*pi*G*(rc**3)/sqrt(rc**2+r2)
    end function GetPloytrop5Solution

    real(kind=double) function GetCondensedSphereSolution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: rc, d, alpha, beta
        rc = 0.1_rp
        d = sqrt(xcord**2 + ycord**2 + zcord**2)
        alpha = d/rc
        if(alpha .eq. 0.0_rp) then
            beta = 0.0_rp
        else
            beta = atan(alpha)/alpha
        end if
        GetCondensedSphereSolution = 4.0_rp*pi*G*(rc**2)*(beta + &
            0.5_rp*log((1+alpha**2)/(1+(R/rc)**2)) - 1.0_rp)
    end function GetCondensedSphereSolution

    real(kind=double) function GetOCCondensedSolution(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord
        real(kind=double) :: left = 1.0_rp
        real(kind=double) :: right = 1.0_rp
        real(kind=double) :: rc = 0.1_rp
        real(kind=double) :: rl, rr, alphal, alphar, betal, betar
        rl = sqrt((xcord + 0.0_rp)**2 + (ycord + 0.5_rp)**2 + (zcord + 0.0_rp)**2)
        rr = sqrt((xcord - 0.0_rp)**2 + (ycord - 0.5_rp)**2 + (zcord - 0.0_rp)**2)
        alphal = rl/rc
        alphar = rr/rc
        if(alphal .eq. 0.0_rp) then
            betal = 0.0_rp
        else
            betal = atan(alphal)/alphal
        end if
        if(alphar .eq. 0.0_rp) then
            betar = 0.0_rp
        else
            betar = atan(alphar)/alphar
        end if
        GetOCCondensedSolution = left*4.0_rp*pi*G*(rc**2)*(betal + &
            0.5_rp*log((1+alphal**2)/(1+(R/rc)**2)) - 1.0_rp) + &
            right*4.0_rp*pi*G*(rc**2)*(betar + &
            0.5_rp*log((1+alphar**2)/(1+(R/rc)**2)) - 1.0_rp)
    end function GetOCCondensedSolution

    real(kind=double) function GetBoundaryValue(xcord,ycord,zcord)
        implicit none
        real(kind=double), intent(in) :: xcord,ycord,zcord

        if(scenerio .eq. 0) then
            GetBoundaryValue = GetVaccumSolution()
        else if(scenerio .eq. 1) then
            GetBoundaryValue = GetR2Solution(xcord,ycord,zcord)
        else if(scenerio .eq. 2) then
            GetBoundaryValue = GetOCGaussianSolution(xcord,ycord,zcord)
        else if(scenerio .eq. 3) then
            GetBoundaryValue = GetR2Solution(xcord,ycord,zcord)
        else if(scenerio .eq. 4) then
            GetBoundaryValue =  GetSinSolution(xcord,ycord,zcord)
        else if(scenerio .eq. 5) then
            GetBoundaryValue = GetPoly6Solution(xcord,ycord,zcord)
        else if(scenerio .eq. 6) then
            GetBoundaryValue = GetGaussianSolution(xcord,ycord,zcord)
        else if(scenerio .eq. 7) then
            GetBoundaryValue = GetPloytrop1Solution(xcord,ycord,zcord)
        else if(scenerio .eq. 8) then
            GetBoundaryValue = GetPloytrop5Solution(xcord,ycord,zcord)
        else if(scenerio .eq. 9) then
            GetBoundaryValue = GetCondensedSphereSolution(xcord,ycord,zcord)
        else if(scenerio .eq. 10) then
            GetBoundaryValue = GetOCCondensedSolution(xcord,ycord,zcord)
        end if

    end function GetBoundaryValue


END MODULE MD_BoundaryCondition
