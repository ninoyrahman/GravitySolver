!****************************************************************
! CalculateStencil.f90
!
! Created on: Jul 26, 2014
! Author: ninoy
!*****************************************************************

subroutine CalculateStencil(i,j,k)
    use MD_Parameter
    use MD_Quantity
    use MD_GeometricQuantity
    implicit none

    integer , intent(in) :: i,j,k
    integer :: nghi,nghj,nghk

    print *, ''
    print *, '********************'
    print *, 'Calculating Stencil'
    print *, '********************'
    print *, ''


    do nghi = -swidth, +swidth
        do nghj = -swidth, +swidth
            do nghk = -swidth, +swidth
                S(nghi,nghj,nghk) = 0.0
            end do
        end do
    end do


    S(-2,0,0) = h*((-5*bxx(i-1,j,k))/144. - (9*byx(i,j-1,k))/256. + (9*byx(i,j,k))/256. - &
        (9*bzx(i,j,k-1))/256. + (9*bzx(i,j,k))/256.) + h**2*(-Gxy(i-1,j,k)/128. - Gxz(i-1,j,k)/128. - &
        Gyx(i,j-1,k)/256. + Gyx(i,j,k)/256. - Gzx(i,j,k-1)/256. + Gzx(i,j,k)/256.)

    S(-1,0,0) = h*((15*bxx(i-1,j,k))/16. + (5*bxx(i,j,k))/144. + (45*byx(i,j-1,k))/128. - &
        (45*byx(i,j,k))/128. - (3*byy(i,j-1,k))/64. - (3*byy(i,j,k))/64. + (45*bzx(i,j,k-1))/128. - &
        (45*bzx(i,j,k))/128. - (3*bzz(i,j,k-1))/64. - (3*bzz(i,j,k))/64.) + h**2*((9*Gxy(i-1,j,k))/128. + &
        Gxy(i,j,k)/128. + (9*Gxz(i-1,j,k))/128. + Gxz(i,j,k)/128. - Gyx(i,j-1,k)/32. + &
        Gyx(i,j,k)/32. - Gzx(i,j,k-1)/32. + Gzx(i,j,k)/32.)

    S(0,0,0) = h*((-15*bxx(i-1,j,k))/16. - (15*bxx(i,j,k))/16. - (15*byy(i,j-1,k))/16. - (15*byy(i,j,k))/16. - &
        (15*bzz(i,j,k-1))/16. - (15*bzz(i,j,k))/16.) + h**2*((9*Gxy(i-1,j,k))/128. - (9*Gxy(i,j,k))/128. + &
        (9*Gxz(i-1,j,k))/128. - (9*Gxz(i,j,k))/128. + (9*Gyx(i,j-1,k))/128. - (9*Gyx(i,j,k))/128. + &
        (9*Gyz(i,j-1,k))/128. - (9*Gyz(i,j,k))/128. + (9*Gzx(i,j,k-1))/128. - (9*Gzx(i,j,k))/128. + &
        (9*Gzy(i,j,k-1))/128. - (9*Gzy(i,j,k))/128.)

    S(1,0,0) = h*((5*bxx(i-1,j,k))/144. + (15*bxx(i,j,k))/16. - (45*byx(i,j-1,k))/128. + &
        (45*byx(i,j,k))/128. - (3*byy(i,j-1,k))/64. - (3*byy(i,j,k))/64. - (45*bzx(i,j,k-1))/128. + &
        (45*bzx(i,j,k))/128. - (3*bzz(i,j,k-1))/64. - (3*bzz(i,j,k))/64.) + h**2*(-Gxy(i-1,j,k)/128. - &
        (9*Gxy(i,j,k))/128. - Gxz(i-1,j,k)/128. - (9*Gxz(i,j,k))/128. - Gyx(i,j-1,k)/32. + &
        Gyx(i,j,k)/32. - Gzx(i,j,k-1)/32. + Gzx(i,j,k)/32.)

    S(2,0,0) = h*((-5*bxx(i,j,k))/144. + (9*byx(i,j-1,k))/256. - (9*byx(i,j,k))/256. + &
        (9*bzx(i,j,k-1))/256. - (9*bzx(i,j,k))/256.) + h**2*(Gxy(i,j,k)/128. + Gxz(i,j,k)/128. - &
        Gyx(i,j-1,k)/256. + Gyx(i,j,k)/256. - Gzx(i,j,k-1)/256. + Gzx(i,j,k)/256.)

    S(-2,1,0) = h*(-bxx(i-1,j,k)/576. + (5*bxy(i-1,j,k))/128. + byx(i,j-1,k)/256. + &
        (9*byx(i,j,k))/256.) + h**2*(Gxy(i-1,j,k)/288. + Gyx(i,j-1,k)/2304. + Gyx(i,j,k)/256.)

    S(-1,1,0) = h*((3*bxx(i-1,j,k))/64. + bxx(i,j,k)/576. - (45*bxy(i-1,j,k))/128. - &
        (5*bxy(i,j,k))/128. - (5*byx(i,j-1,k))/128. - (45*byx(i,j,k))/128. + byy(i,j-1,k)/576. + &
        (3*byy(i,j,k))/64.) + h**2*(-Gxy(i-1,j,k)/32. - Gxy(i,j,k)/288. + Gyx(i,j-1,k)/288. + &
        Gyx(i,j,k)/32.)

    S(0,1,0) = h*((-3*bxx(i-1,j,k))/64. - (3*bxx(i,j,k))/64. - (45*bxy(i-1,j,k))/128. + &
        (45*bxy(i,j,k))/128. + (5*byy(i,j-1,k))/144. + (15*byy(i,j,k))/16. - (45*bzy(i,j,k-1))/128. + &
        (45*bzy(i,j,k))/128. - (3*bzz(i,j,k-1))/64. - (3*bzz(i,j,k))/64.) + h**2*(-Gxy(i-1,j,k)/32. + &
        Gxy(i,j,k)/32. - Gyx(i,j-1,k)/128. - (9*Gyx(i,j,k))/128. - Gyz(i,j-1,k)/128. - &
        (9*Gyz(i,j,k))/128. - Gzy(i,j,k-1)/32. + Gzy(i,j,k)/32.)

    S(1,1,0) = h*(bxx(i-1,j,k)/576. + (3*bxx(i,j,k))/64. + (5*bxy(i-1,j,k))/128. + &
        (45*bxy(i,j,k))/128. + (5*byx(i,j-1,k))/128. + (45*byx(i,j,k))/128. + byy(i,j-1,k)/576. + &
        (3*byy(i,j,k))/64.) + h**2*(Gxy(i-1,j,k)/288. + Gxy(i,j,k)/32. + &
        Gyx(i,j-1,k)/288. + Gyx(i,j,k)/32.)

    S(2,1,0) = h*(-bxx(i,j,k)/576. - (5*bxy(i,j,k))/128. - byx(i,j-1,k)/256. - (9*byx(i,j,k))/256.) + &
        h**2*(-Gxy(i,j,k)/288. + Gyx(i,j-1,k)/2304. + Gyx(i,j,k)/256.)


    S(-2,2,0) = h*(-bxy(i-1,j,k)/256. - byx(i,j,k)/256.) + h**2*(Gxy(i-1,j,k)/2304. - Gyx(i,j,k)/2304.)

    S(-1,2,0) = h*((9*bxy(i-1,j,k))/256. + bxy(i,j,k)/256. + (5*byx(i,j,k))/128. - byy(i,j,k)/576.) + &
        h**2*(-Gxy(i-1,j,k)/256. - Gxy(i,j,k)/2304. - Gyx(i,j,k)/288.)

    S(0,2,0) = h*((9*bxy(i-1,j,k))/256. - (9*bxy(i,j,k))/256. - (5*byy(i,j,k))/144. + &
        (9*bzy(i,j,k-1))/256. - (9*bzy(i,j,k))/256.) + h**2*(-Gxy(i-1,j,k)/256. + &
        Gxy(i,j,k)/256. + Gyx(i,j,k)/128. + Gyz(i,j,k)/128. - &
        Gzy(i,j,k-1)/256. + Gzy(i,j,k)/256.)

    S(1,2,0) = h*(-bxy(i-1,j,k)/256. - (9*bxy(i,j,k))/256. - (5*byx(i,j,k))/128. - byy(i,j,k)/576.) + &
        h**2*(Gxy(i-1,j,k)/2304. + Gxy(i,j,k)/256. - Gyx(i,j,k)/288.)

    S(2,2,0) = h*(bxy(i,j,k)/256. + byx(i,j,k)/256.) + h**2*(-Gxy(i,j,k)/2304. - Gyx(i,j,k)/2304.)

    S(-2,-1,0) = h*(-bxx(i-1,j,k)/576. - (5*bxy(i-1,j,k))/128. - (9*byx(i,j-1,k))/256. - &
        byx(i,j,k)/256.) + h**2*(Gxy(i-1,j,k)/288. - Gyx(i,j-1,k)/256. - Gyx(i,j,k)/2304.)

    S(-1,-1,0) = h*((3*bxx(i-1,j,k))/64. + bxx(i,j,k)/576. + (45*bxy(i-1,j,k))/128. + &
        (5*bxy(i,j,k))/128. + (45*byx(i,j-1,k))/128. + (5*byx(i,j,k))/128. + (3*byy(i,j-1,k))/64. + &
        byy(i,j,k)/576.) + h**2*(-Gxy(i-1,j,k)/32. - Gxy(i,j,k)/288. - &
        Gyx(i,j-1,k)/32. - Gyx(i,j,k)/288.)

    S(0,-1,0) = h*((-3*bxx(i-1,j,k))/64. - (3*bxx(i,j,k))/64. + (45*bxy(i-1,j,k))/128. - &
        (45*bxy(i,j,k))/128. + (15*byy(i,j-1,k))/16. + (5*byy(i,j,k))/144. + (45*bzy(i,j,k-1))/128. - &
        (45*bzy(i,j,k))/128. - (3*bzz(i,j,k-1))/64. - (3*bzz(i,j,k))/64.) + h**2*(-Gxy(i-1,j,k)/32. + &
        Gxy(i,j,k)/32. + (9*Gyx(i,j-1,k))/128. + Gyx(i,j,k)/128. + (9*Gyz(i,j-1,k))/128. + &
        Gyz(i,j,k)/128. - Gzy(i,j,k-1)/32. + Gzy(i,j,k)/32.)

    S(1,-1,0) = h*(bxx(i-1,j,k)/576. + (3*bxx(i,j,k))/64. - (5*bxy(i-1,j,k))/128. - (45*bxy(i,j,k))/128. - &
        (45*byx(i,j-1,k))/128. - (5*byx(i,j,k))/128. + (3*byy(i,j-1,k))/64. + byy(i,j,k)/576.) + &
        h**2*(Gxy(i-1,j,k)/288. + Gxy(i,j,k)/32. - Gyx(i,j-1,k)/32. - Gyx(i,j,k)/288.)

    S(2,-1,0) = h*(-bxx(i,j,k)/576. + (5*bxy(i,j,k))/128. + (9*byx(i,j-1,k))/256. + byx(i,j,k)/256.) + &
        h**2*(-Gxy(i,j,k)/288. - Gyx(i,j-1,k)/256. - Gyx(i,j,k)/2304.)

    S(-2,-2,0) = h*(bxy(i-1,j,k)/256. + byx(i,j-1,k)/256.) + h**2*(Gxy(i-1,j,k)/2304. + Gyx(i,j-1,k)/2304.)

    S(-1,-2,0) = h*((-9*bxy(i-1,j,k))/256. - bxy(i,j,k)/256. - (5*byx(i,j-1,k))/128. - byy(i,j-1,k)/576.) + &
        h**2*(-Gxy(i-1,j,k)/256. - Gxy(i,j,k)/2304. + Gyx(i,j-1,k)/288.)

    S(0,-2,0) = h*((-9*bxy(i-1,j,k))/256. + (9*bxy(i,j,k))/256. - (5*byy(i,j-1,k))/144. - &
        (9*bzy(i,j,k-1))/256. + (9*bzy(i,j,k))/256.) + h**2*(-Gxy(i-1,j,k)/256. + Gxy(i,j,k)/256. - &
        Gyx(i,j-1,k)/128. - Gyz(i,j-1,k)/128. - Gzy(i,j,k-1)/256. + Gzy(i,j,k)/256.)

    S(1,-2,0) = h*(bxy(i-1,j,k)/256. + (9*bxy(i,j,k))/256. + (5*byx(i,j-1,k))/128. - &
        byy(i,j-1,k)/576.) + h**2*(Gxy(i-1,j,k)/2304. + Gxy(i,j,k)/256. + Gyx(i,j-1,k)/288.)

    S(2,-2,0) = h*(-bxy(i,j,k)/256. - byx(i,j-1,k)/256.) + h**2*(-Gxy(i,j,k)/2304. + Gyx(i,j-1,k)/2304.)

    S(-2,0,1) = h*(-bxx(i-1,j,k)/576. + (5*bxz(i-1,j,k))/128. + bzx(i,j,k-1)/256. + (9*bzx(i,j,k))/256.) + &
        h**2*(Gxz(i-1,j,k)/288. + Gzx(i,j,k-1)/2304. + Gzx(i,j,k)/256.)

    S(-1,0,1) = h*((3*bxx(i-1,j,k))/64. + bxx(i,j,k)/576. - (45*bxz(i-1,j,k))/128. - &
        (5*bxz(i,j,k))/128. - (5*bzx(i,j,k-1))/128. - (45*bzx(i,j,k))/128. + bzz(i,j,k-1)/576. + &
        (3*bzz(i,j,k))/64.) + h**2*(-Gxz(i-1,j,k)/32. - Gxz(i,j,k)/288. + Gzx(i,j,k-1)/288. + &
        Gzx(i,j,k)/32.)

    S(1,0,1) = h*(bxx(i-1,j,k)/576. + (3*bxx(i,j,k))/64. + (5*bxz(i-1,j,k))/128. + (45*bxz(i,j,k))/128. + &
        (5*bzx(i,j,k-1))/128. + (45*bzx(i,j,k))/128. + bzz(i,j,k-1)/576. + (3*bzz(i,j,k))/64.) + &
        h**2*(Gxz(i-1,j,k)/288. + Gxz(i,j,k)/32. + Gzx(i,j,k-1)/288. + Gzx(i,j,k)/32.)

    S(2,0,1) = h*(-bxx(i,j,k)/576. - (5*bxz(i,j,k))/128. - bzx(i,j,k-1)/256. - (9*bzx(i,j,k))/256.) + &
        h**2*(-Gxz(i,j,k)/288. + Gzx(i,j,k-1)/2304. + Gzx(i,j,k)/256.)

    S(-2,0,2) = h*(-bxz(i-1,j,k)/256. - bzx(i,j,k)/256.) + h**2*(Gxz(i-1,j,k)/2304. - Gzx(i,j,k)/2304.)

    S(-1,0,2) = h*((9*bxz(i-1,j,k))/256. + bxz(i,j,k)/256. + (5*bzx(i,j,k))/128. - bzz(i,j,k)/576.) + &
        h**2*(-Gxz(i-1,j,k)/256. - Gxz(i,j,k)/2304. - Gzx(i,j,k)/288.)

    S(1,0,2) = h*(-bxz(i-1,j,k)/256. - (9*bxz(i,j,k))/256. - (5*bzx(i,j,k))/128. - bzz(i,j,k)/576.) + &
        h**2*(Gxz(i-1,j,k)/2304. + Gxz(i,j,k)/256. - Gzx(i,j,k)/288.)

    S(2,0,2) = h*(bxz(i,j,k)/256. + bzx(i,j,k)/256.) + h**2*(-Gxz(i,j,k)/2304. - Gzx(i,j,k)/2304.)

    S(-2,0,-1) = h*(-bxx(i-1,j,k)/576. - (5*bxz(i-1,j,k))/128. - (9*bzx(i,j,k-1))/256. - bzx(i,j,k)/256.) + &
        h**2*(Gxz(i-1,j,k)/288. - Gzx(i,j,k-1)/256. - Gzx(i,j,k)/2304.)

    S(-1,0,-1) = h*((3*bxx(i-1,j,k))/64. + bxx(i,j,k)/576. + (45*bxz(i-1,j,k))/128. + (5*bxz(i,j,k))/128. + &
        (45*bzx(i,j,k-1))/128. + (5*bzx(i,j,k))/128. + (3*bzz(i,j,k-1))/64. + bzz(i,j,k)/576.) + &
        h**2*(-Gxz(i-1,j,k)/32. - Gxz(i,j,k)/288. - Gzx(i,j,k-1)/32. - Gzx(i,j,k)/288.)

    S(0,0,-1) = h*((-3*bxx(i-1,j,k))/64. - (3*bxx(i,j,k))/64. + (45*bxz(i-1,j,k))/128. - &
        (45*bxz(i,j,k))/128. - (3*byy(i,j-1,k))/64. - (3*byy(i,j,k))/64. + (45*byz(i,j-1,k))/128. - &
        (45*byz(i,j,k))/128. + (15*bzz(i,j,k-1))/16. + (5*bzz(i,j,k))/144.) + h**2*(-Gxz(i-1,j,k)/32. + &
        Gxz(i,j,k)/32. - Gyz(i,j-1,k)/32. + Gyz(i,j,k)/32. + (9*Gzx(i,j,k-1))/128. + &
        Gzx(i,j,k)/128. + (9*Gzy(i,j,k-1))/128. + Gzy(i,j,k)/128.)

    S(1,0,-1) = h*(bxx(i-1,j,k)/576. + (3*bxx(i,j,k))/64. - (5*bxz(i-1,j,k))/128. - (45*bxz(i,j,k))/128. - &
        (45*bzx(i,j,k-1))/128. - (5*bzx(i,j,k))/128. + (3*bzz(i,j,k-1))/64. + bzz(i,j,k)/576.) + &
        h**2*(Gxz(i-1,j,k)/288. + Gxz(i,j,k)/32. - Gzx(i,j,k-1)/32. - Gzx(i,j,k)/288.)

    S(2,0,-1) = h*(-bxx(i,j,k)/576. + (5*bxz(i,j,k))/128. + (9*bzx(i,j,k-1))/256. + bzx(i,j,k)/256.) + &
        h**2*(-Gxz(i,j,k)/288. - Gzx(i,j,k-1)/256. - Gzx(i,j,k)/2304.)

    S(-2,0,-2) = h*(bxz(i-1,j,k)/256. + bzx(i,j,k-1)/256.) + h**2*(Gxz(i-1,j,k)/2304. + Gzx(i,j,k-1)/2304.)

    S(-1,0,-2) = h*((-9*bxz(i-1,j,k))/256. - bxz(i,j,k)/256. - (5*bzx(i,j,k-1))/128. - bzz(i,j,k-1)/576.) + &
        h**2*(-Gxz(i-1,j,k)/256. - Gxz(i,j,k)/2304. + Gzx(i,j,k-1)/288.)

    S(0,0,-2) = h*((-9*bxz(i-1,j,k))/256. + (9*bxz(i,j,k))/256. - (9*byz(i,j-1,k))/256. + &
        (9*byz(i,j,k))/256. - (5*bzz(i,j,k-1))/144.) + h**2*(-Gxz(i-1,j,k)/256. + Gxz(i,j,k)/256. - &
        Gyz(i,j-1,k)/256. + Gyz(i,j,k)/256. - Gzx(i,j,k-1)/128. - Gzy(i,j,k-1)/128.)

    S(1,0,-2) = h*(bxz(i-1,j,k)/256. + (9*bxz(i,j,k))/256. + (5*bzx(i,j,k-1))/128. - bzz(i,j,k-1)/576.) + &
        h**2*(Gxz(i-1,j,k)/2304. + Gxz(i,j,k)/256. + Gzx(i,j,k-1)/288.)

    S(2,0,-2) = h*(-bxz(i,j,k)/256. - bzx(i,j,k-1)/256.) + h**2*(-Gxz(i,j,k)/2304. + Gzx(i,j,k-1)/2304.)

    S(0,-2,1) = h*(-byy(i,j-1,k)/576. + (5*byz(i,j-1,k))/128. + bzy(i,j,k-1)/256. + &
        (9*bzy(i,j,k))/256.) + h**2*(Gyz(i,j-1,k)/288. + Gzy(i,j,k-1)/2304. + Gzy(i,j,k)/256.)

    S(0,-1,1) = h*((3*byy(i,j-1,k))/64. + byy(i,j,k)/576. - (45*byz(i,j-1,k))/128. - (5*byz(i,j,k))/128. - &
        (5*bzy(i,j,k-1))/128. - (45*bzy(i,j,k))/128. + bzz(i,j,k-1)/576. + (3*bzz(i,j,k))/64.) + &
        h**2*(-Gyz(i,j-1,k)/32. - Gyz(i,j,k)/288. + Gzy(i,j,k-1)/288. + Gzy(i,j,k)/32.)

    S(0,0,1) = h*((-3*bxx(i-1,j,k))/64. - (3*bxx(i,j,k))/64. - (45*bxz(i-1,j,k))/128. + (45*bxz(i,j,k))/128. - &
        (3*byy(i,j-1,k))/64. - (3*byy(i,j,k))/64. - (45*byz(i,j-1,k))/128. + (45*byz(i,j,k))/128. + &
        (5*bzz(i,j,k-1))/144. + (15*bzz(i,j,k))/16.) + h**2*(-Gxz(i-1,j,k)/32. + Gxz(i,j,k)/32. - &
        Gyz(i,j-1,k)/32. + Gyz(i,j,k)/32. - Gzx(i,j,k-1)/128. - (9*Gzx(i,j,k))/128. - Gzy(i,j,k-1)/128. - &
        (9*Gzy(i,j,k))/128.)

    S(0,1,1) = h*(byy(i,j-1,k)/576. + (3*byy(i,j,k))/64. + (5*byz(i,j-1,k))/128. + (45*byz(i,j,k))/128. + &
        (5*bzy(i,j,k-1))/128. + (45*bzy(i,j,k))/128. + bzz(i,j,k-1)/576. + (3*bzz(i,j,k))/64.) + &
        h**2*(Gyz(i,j-1,k)/288. + Gyz(i,j,k)/32. + Gzy(i,j,k-1)/288. + Gzy(i,j,k)/32.)

    S(0,2,1) = h*(-byy(i,j,k)/576. - (5*byz(i,j,k))/128. - bzy(i,j,k-1)/256. - (9*bzy(i,j,k))/256.) + &
        h**2*(-Gyz(i,j,k)/288. + Gzy(i,j,k-1)/2304. + Gzy(i,j,k)/256.)

    S(0,-2,2) = h*(-byz(i,j-1,k)/256. - bzy(i,j,k)/256.) + h**2*(Gyz(i,j-1,k)/2304. - Gzy(i,j,k)/2304.)

    S(0,-1,2) = h*((9*byz(i,j-1,k))/256. + byz(i,j,k)/256. + (5*bzy(i,j,k))/128. - bzz(i,j,k)/576.) + &
        h**2*(-Gyz(i,j-1,k)/256. - Gyz(i,j,k)/2304. - Gzy(i,j,k)/288.)

    S(0,0,2) = h*((9*bxz(i-1,j,k))/256. - (9*bxz(i,j,k))/256. + (9*byz(i,j-1,k))/256. - (9*byz(i,j,k))/256. - &
        (5*bzz(i,j,k))/144.) + h**2*(-Gxz(i-1,j,k)/256. + Gxz(i,j,k)/256. - Gyz(i,j-1,k)/256. + &
        Gyz(i,j,k)/256. + Gzx(i,j,k)/128. + Gzy(i,j,k)/128.)

    S(0,1,2) = h*(-byz(i,j-1,k)/256. - (9*byz(i,j,k))/256. - (5*bzy(i,j,k))/128. - bzz(i,j,k)/576.) + &
        h**2*(Gyz(i,j-1,k)/2304. + Gyz(i,j,k)/256. - Gzy(i,j,k)/288.)

    S(0,2,2) = h*(byz(i,j,k)/256. + bzy(i,j,k)/256.) + h**2*(-Gyz(i,j,k)/2304. - Gzy(i,j,k)/2304.)

    S(0,-2,-1) = h*(-byy(i,j-1,k)/576. - (5*byz(i,j-1,k))/128. - (9*bzy(i,j,k-1))/256. - bzy(i,j,k)/256.) + &
        h**2*(Gyz(i,j-1,k)/288. - Gzy(i,j,k-1)/256. - Gzy(i,j,k)/2304.)

    S(0,-1,-1) = h*((3*byy(i,j-1,k))/64. + byy(i,j,k)/576. + (45*byz(i,j-1,k))/128. + (5*byz(i,j,k))/128. + &
        (45*bzy(i,j,k-1))/128. + (5*bzy(i,j,k))/128. + (3*bzz(i,j,k-1))/64. + bzz(i,j,k)/576.) + &
        h**2*(-Gyz(i,j-1,k)/32. - Gyz(i,j,k)/288. - Gzy(i,j,k-1)/32. - Gzy(i,j,k)/288.)

    S(0,1,-1) = h*(byy(i,j-1,k)/576. + (3*byy(i,j,k))/64. - (5*byz(i,j-1,k))/128. - (45*byz(i,j,k))/128. - &
        (45*bzy(i,j,k-1))/128. - (5*bzy(i,j,k))/128. + (3*bzz(i,j,k-1))/64. + bzz(i,j,k)/576.) + &
        h**2*(Gyz(i,j-1,k)/288. + Gyz(i,j,k)/32. - Gzy(i,j,k-1)/32. - Gzy(i,j,k)/288.)

    S(0,2,-1) = h*(-byy(i,j,k)/576. + (5*byz(i,j,k))/128. + (9*bzy(i,j,k-1))/256. + bzy(i,j,k)/256.) + &
        h**2*(-Gyz(i,j,k)/288. - Gzy(i,j,k-1)/256. - Gzy(i,j,k)/2304.)

    S(0,-2,-2) = h*(byz(i,j-1,k)/256. + bzy(i,j,k-1)/256.) + h**2*(Gyz(i,j-1,k)/2304. + Gzy(i,j,k-1)/2304.)

    S(0,-1,-2) = h*((-9*byz(i,j-1,k))/256. - byz(i,j,k)/256. - (5*bzy(i,j,k-1))/128. - bzz(i,j,k-1)/576.) + &
        h**2*(-Gyz(i,j-1,k)/256. - Gyz(i,j,k)/2304. + Gzy(i,j,k-1)/288.)

    S(0,1,-2) = h*(byz(i,j-1,k)/256. + (9*byz(i,j,k))/256. + (5*bzy(i,j,k-1))/128. - bzz(i,j,k-1)/576.) + &
        h**2*(Gyz(i,j-1,k)/2304. + Gyz(i,j,k)/256. + Gzy(i,j,k-1)/288.)

    S(0,2,-2) = h*(-byz(i,j,k)/256. - bzy(i,j,k-1)/256.) + h**2*(-Gyz(i,j,k)/2304. + Gzy(i,j,k-1)/2304.)

end subroutine CalculateStencil
