!****************************************************************
! ComputeStencil4.f90
!
! Created on: Jul 26, 2014
! Author: ninoy
!*****************************************************************

subroutine ComputeStencil4
    use MD_Parameter
    use MD_Quantity
    use MD_GeometricQuantity
    use MD_Helper
    implicit none

    integer :: i,j,k, row, counter
    integer :: nghi, nghj, nghk
    character(len=100) :: format
    real(kind=double), dimension(:,:,:), allocatable :: S

    if(rank .eq. master) then
        print *, ''
        print *, '*****************************'
        print *, 'Calculating Stencil 4th-Order'
        print *, '*****************************'
        print *, ''
    endif

    if(.not. allocated(S)) allocate(S(-swidth:swidth,-swidth:swidth,-swidth:swidth))

    counter = 0
    do nghk = -swidth, +swidth
        do nghj = -swidth, +swidth
            do nghi = -swidth, +swidth
                stencilindex(nghi,nghj,nghk) = 0
                if((nghi .eq. 0) .or. (nghj .eq. 0) .or. (nghk .eq. 0)) then
                    counter = counter + 1
                    stencilindex(nghi,nghj,nghk) = counter
                end if
            end do
        end do
    end do


    format = "(F14.7)"

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,nghi,nghj,nghk,row,S)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax

                row = localindex(i,j,k)

                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            S(nghi,nghj,nghk) = 0.0_rp
                        end do
                    end do
                end do

                S(-2,0,0) = h*((-5.0_rp*bxx(i-1,j,k))/144.0_rp - (9.0_rp*byx(i,j-1,k))/256.0_rp + &
                    (9.0_rp*byx(i,j,k))/256.0_rp - &
                    (9.0_rp*bzx(i,j,k-1))/256.0_rp + (9.0_rp*bzx(i,j,k))/256.0_rp) + &
                    h**2*(-Gxy(i-1,j,k)/128.0_rp - Gxz(i-1,j,k)/128.0_rp - &
                    Gyx(i,j-1,k)/256.0_rp + Gyx(i,j,k)/256.0_rp - Gzx(i,j,k-1)/256.0_rp + Gzx(i,j,k)/256.0_rp)

                S(-1,0,0) = h*((15.0_rp*bxx(i-1,j,k))/16.0_rp + (5.0_rp*bxx(i,j,k))/144.0_rp + &
                    (45.0_rp*byx(i,j-1,k))/128.0_rp - &
                    (45.0_rp*byx(i,j,k))/128.0_rp - (3.0_rp*byy(i,j-1,k))/64.0_rp - &
                    (3.0_rp*byy(i,j,k))/64.0_rp + (45.0_rp*bzx(i,j,k-1))/128.0_rp - &
                    (45.0_rp*bzx(i,j,k))/128.0_rp - (3.0_rp*bzz(i,j,k-1))/64.0_rp - &
                    (3.0_rp*bzz(i,j,k))/64.0_rp) + &
                    h**2*((9.0_rp*Gxy(i-1,j,k))/128.0_rp + &
                    Gxy(i,j,k)/128.0_rp + (9.0_rp*Gxz(i-1,j,k))/128.0_rp + &
                    Gxz(i,j,k)/128.0_rp - Gyx(i,j-1,k)/32.0_rp + &
                    Gyx(i,j,k)/32.0_rp - Gzx(i,j,k-1)/32.0_rp + Gzx(i,j,k)/32.0_rp)

                S(0,0,0) = h*((-15.0_rp*bxx(i-1,j,k))/16.0_rp - (15.0_rp*bxx(i,j,k))/16.0_rp - &
                    (15.0_rp*byy(i,j-1,k))/16.0_rp - (15.0_rp*byy(i,j,k))/16.0_rp - &
                    (15.0_rp*bzz(i,j,k-1))/16.0_rp - (15.0_rp*bzz(i,j,k))/16.0_rp) + &
                    h**2*((9.0_rp*Gxy(i-1,j,k))/128.0_rp - (9.0_rp*Gxy(i,j,k))/128.0_rp + &
                    (9.0_rp*Gxz(i-1,j,k))/128.0_rp - (9.0_rp*Gxz(i,j,k))/128.0_rp + &
                    (9.0_rp*Gyx(i,j-1,k))/128.0_rp - (9.0_rp*Gyx(i,j,k))/128.0_rp + &
                    (9.0_rp*Gyz(i,j-1,k))/128.0_rp - (9.0_rp*Gyz(i,j,k))/128.0_rp + &
                    (9.0_rp*Gzx(i,j,k-1))/128.0_rp - (9.0_rp*Gzx(i,j,k))/128.0_rp + &
                    (9.0_rp*Gzy(i,j,k-1))/128.0_rp - (9.0_rp*Gzy(i,j,k))/128.0_rp)

                S(1,0,0) = h*((5.0_rp*bxx(i-1,j,k))/144.0_rp + (15.0_rp*bxx(i,j,k))/16.0_rp - &
                    (45.0_rp*byx(i,j-1,k))/128.0_rp + &
                    (45.0_rp*byx(i,j,k))/128.0_rp - (3.0_rp*byy(i,j-1,k))/64.0_rp - &
                    (3.0_rp*byy(i,j,k))/64.0_rp - (45.0_rp*bzx(i,j,k-1))/128.0_rp + &
                    (45.0_rp*bzx(i,j,k))/128.0_rp - (3.0_rp*bzz(i,j,k-1))/64.0_rp - &
                    (3.0_rp*bzz(i,j,k))/64.0_rp) + h**2*(-Gxy(i-1,j,k)/128.0_rp - &
                    (9.0_rp*Gxy(i,j,k))/128.0_rp - Gxz(i-1,j,k)/128.0_rp - &
                    (9.0_rp*Gxz(i,j,k))/128.0_rp - Gyx(i,j-1,k)/32.0_rp + &
                    Gyx(i,j,k)/32.0_rp - Gzx(i,j,k-1)/32.0_rp + Gzx(i,j,k)/32.0_rp)


                S(2,0,0) = h*((-5.0_rp*bxx(i,j,k))/144.0_rp + (9.0_rp*byx(i,j-1,k))/256.0_rp - &
                    (9.0_rp*byx(i,j,k))/256.0_rp + &
                    (9.0_rp*bzx(i,j,k-1))/256.0_rp - (9.0_rp*bzx(i,j,k))/256.0_rp) + &
                    h**2*(Gxy(i,j,k)/128.0_rp + Gxz(i,j,k)/128.0_rp - &
                    Gyx(i,j-1,k)/256.0_rp + Gyx(i,j,k)/256.0_rp - &
                    Gzx(i,j,k-1)/256.0_rp + Gzx(i,j,k)/256.0_rp)

                S(-2,1,0) = h*(-bxx(i-1,j,k)/576.0_rp + (5.0_rp*bxy(i-1,j,k))/128.0_rp + &
                    byx(i,j-1,k)/256.0_rp + &
                    (9.0_rp*byx(i,j,k))/256.0_rp) + h**2*(Gxy(i-1,j,k)/288.0_rp + &
                    Gyx(i,j-1,k)/2304.0_rp + Gyx(i,j,k)/256.0_rp)

                S(-1,1,0) = h*((3.0_rp*bxx(i-1,j,k))/64.0_rp + bxx(i,j,k)/576.0_rp - &
                    (45.0_rp*bxy(i-1,j,k))/128.0_rp - &
                    (5.0_rp*bxy(i,j,k))/128.0_rp - (5.0_rp*byx(i,j-1,k))/128.0_rp - &
                    (45.0_rp*byx(i,j,k))/128.0_rp + byy(i,j-1,k)/576.0_rp + &
                    (3.0_rp*byy(i,j,k))/64.0_rp) + h**2*(-Gxy(i-1,j,k)/32.0_rp - &
                    Gxy(i,j,k)/288.0_rp + Gyx(i,j-1,k)/288.0_rp + &
                    Gyx(i,j,k)/32.0_rp)

                S(0,1,0) = h*((-3.0_rp*bxx(i-1,j,k))/64.0_rp - (3.0_rp*bxx(i,j,k))/64.0_rp - &
                    (45.0_rp*bxy(i-1,j,k))/128.0_rp + &
                    (45.0_rp*bxy(i,j,k))/128.0_rp + (5.0_rp*byy(i,j-1,k))/144.0_rp + &
                    (15.0_rp*byy(i,j,k))/16.0_rp - (45.0_rp*bzy(i,j,k-1))/128.0_rp + &
                    (45.0_rp*bzy(i,j,k))/128.0_rp - (3.0_rp*bzz(i,j,k-1))/64.0_rp - &
                    (3.0_rp*bzz(i,j,k))/64.0_rp) + h**2*(-Gxy(i-1,j,k)/32.0_rp + &
                    Gxy(i,j,k)/32.0_rp - Gyx(i,j-1,k)/128.0_rp - (9.0_rp*Gyx(i,j,k))/128.0_rp - &
                    Gyz(i,j-1,k)/128.0_rp - &
                    (9.0_rp*Gyz(i,j,k))/128.0_rp - Gzy(i,j,k-1)/32.0_rp + Gzy(i,j,k)/32.0_rp)

                S(1,1,0) = h*(bxx(i-1,j,k)/576.0_rp + (3.0_rp*bxx(i,j,k))/64.0_rp + &
                    (5.0_rp*bxy(i-1,j,k))/128.0_rp + &
                    (45.0_rp*bxy(i,j,k))/128.0_rp + (5.0_rp*byx(i,j-1,k))/128.0_rp + &
                    (45.0_rp*byx(i,j,k))/128.0_rp + byy(i,j-1,k)/576.0_rp + &
                    (3.0_rp*byy(i,j,k))/64.0_rp) + h**2*(Gxy(i-1,j,k)/288.0_rp + Gxy(i,j,k)/32.0_rp + &
                    Gyx(i,j-1,k)/288.0_rp + Gyx(i,j,k)/32.0_rp)

                S(2,1,0) = h*(-bxx(i,j,k)/576.0_rp - (5.0_rp*bxy(i,j,k))/128.0_rp - &
                    byx(i,j-1,k)/256.0_rp - (9.0_rp*byx(i,j,k))/256.0_rp) + &
                    h**2*(-Gxy(i,j,k)/288.0_rp + Gyx(i,j-1,k)/2304.0_rp + Gyx(i,j,k)/256.0_rp)


                S(-2,2,0) = h*(-bxy(i-1,j,k)/256.0_rp - byx(i,j,k)/256.0_rp) + &
                    h**2*(Gxy(i-1,j,k)/2304.0_rp - Gyx(i,j,k)/2304.0_rp)

                S(-1,2,0) = h*((9.0_rp*bxy(i-1,j,k))/256.0_rp + bxy(i,j,k)/256.0_rp + &
                    (5.0_rp*byx(i,j,k))/128.0_rp - byy(i,j,k)/576.0_rp) + &
                    h**2*(-Gxy(i-1,j,k)/256.0_rp - Gxy(i,j,k)/2304.0_rp - Gyx(i,j,k)/288.0_rp)

                S(0,2,0) = h*((9.0_rp*bxy(i-1,j,k))/256.0_rp - (9.0_rp*bxy(i,j,k))/256.0_rp - &
                    (5.0_rp*byy(i,j,k))/144.0_rp + &
                    (9.0_rp*bzy(i,j,k-1))/256.0_rp - (9.0_rp*bzy(i,j,k))/256.0_rp) + &
                    h**2*(-Gxy(i-1,j,k)/256.0_rp + &
                    Gxy(i,j,k)/256.0_rp + Gyx(i,j,k)/128.0_rp + Gyz(i,j,k)/128.0_rp - &
                    Gzy(i,j,k-1)/256.0_rp + Gzy(i,j,k)/256.0_rp)

                S(1,2,0) = h*(-bxy(i-1,j,k)/256.0_rp - (9.0_rp*bxy(i,j,k))/256.0_rp - &
                    (5.0_rp*byx(i,j,k))/128.0_rp - byy(i,j,k)/576.0_rp) + &
                    h**2*(Gxy(i-1,j,k)/2304.0_rp + Gxy(i,j,k)/256.0_rp - Gyx(i,j,k)/288.0_rp)

                S(2,2,0) = h*(bxy(i,j,k)/256.0_rp + byx(i,j,k)/256.0_rp) + &
                    h**2*(-Gxy(i,j,k)/2304.0_rp - Gyx(i,j,k)/2304.0_rp)

                S(-2,-1,0) = h*(-bxx(i-1,j,k)/576.0_rp - (5.0_rp*bxy(i-1,j,k))/128.0_rp - &
                    (9.0_rp*byx(i,j-1,k))/256.0_rp - &
                    byx(i,j,k)/256.0_rp) + h**2*(Gxy(i-1,j,k)/288.0_rp - &
                    Gyx(i,j-1,k)/256.0_rp - Gyx(i,j,k)/2304.0_rp)

                S(-1,-1,0) = h*((3.0_rp*bxx(i-1,j,k))/64.0_rp + bxx(i,j,k)/576.0_rp + &
                    (45.0_rp*bxy(i-1,j,k))/128.0_rp + &
                    (5.0_rp*bxy(i,j,k))/128.0_rp + (45.0_rp*byx(i,j-1,k))/128.0_rp + &
                    (5.0_rp*byx(i,j,k))/128.0_rp + (3.0_rp*byy(i,j-1,k))/64.0_rp + &
                    byy(i,j,k)/576.0_rp) + h**2*(-Gxy(i-1,j,k)/32.0_rp - Gxy(i,j,k)/288.0_rp - &
                    Gyx(i,j-1,k)/32.0_rp - Gyx(i,j,k)/288.0_rp)

                S(0,-1,0) = h*((-3.0_rp*bxx(i-1,j,k))/64.0_rp - (3.0_rp*bxx(i,j,k))/64.0_rp + &
                    (45.0_rp*bxy(i-1,j,k))/128.0_rp - &
                    (45.0_rp*bxy(i,j,k))/128.0_rp + (15.0_rp*byy(i,j-1,k))/16.0_rp + &
                    (5.0_rp*byy(i,j,k))/144.0_rp + (45.0_rp*bzy(i,j,k-1))/128.0_rp - &
                    (45.0_rp*bzy(i,j,k))/128.0_rp - (3.0_rp*bzz(i,j,k-1))/64.0_rp - &
                    (3.0_rp*bzz(i,j,k))/64.0_rp) + h**2*(-Gxy(i-1,j,k)/32.0_rp + &
                    Gxy(i,j,k)/32.0_rp + (9.0_rp*Gyx(i,j-1,k))/128.0_rp + Gyx(i,j,k)/128.0_rp + &
                    (9.0_rp*Gyz(i,j-1,k))/128.0_rp + &
                    Gyz(i,j,k)/128.0_rp - Gzy(i,j,k-1)/32.0_rp + Gzy(i,j,k)/32.0_rp)

                S(1,-1,0) = h*(bxx(i-1,j,k)/576.0_rp + (3.0_rp*bxx(i,j,k))/64.0_rp - &
                    (5.0_rp*bxy(i-1,j,k))/128.0_rp - (45.0_rp*bxy(i,j,k))/128.0_rp - &
                    (45.0_rp*byx(i,j-1,k))/128.0_rp - (5.0_rp*byx(i,j,k))/128.0_rp + &
                    (3.0_rp*byy(i,j-1,k))/64.0_rp + byy(i,j,k)/576.0_rp) + &
                    h**2*(Gxy(i-1,j,k)/288.0_rp + Gxy(i,j,k)/32.0_rp - Gyx(i,j-1,k)/32.0_rp - &
                    Gyx(i,j,k)/288.0_rp)

                S(2,-1,0) = h*(-bxx(i,j,k)/576.0_rp + (5.0_rp*bxy(i,j,k))/128.0_rp + &
                    (9.0_rp*byx(i,j-1,k))/256.0_rp + byx(i,j,k)/256.0_rp) + &
                    h**2*(-Gxy(i,j,k)/288.0_rp - Gyx(i,j-1,k)/256.0_rp - Gyx(i,j,k)/2304.0_rp)

                S(-2,-2,0) = h*(bxy(i-1,j,k)/256.0_rp + byx(i,j-1,k)/256.0_rp) + &
                    h**2*(Gxy(i-1,j,k)/2304.0_rp + Gyx(i,j-1,k)/2304.0_rp)

                S(-1,-2,0) = h*((-9.0_rp*bxy(i-1,j,k))/256.0_rp - bxy(i,j,k)/256.0_rp - &
                    (5.0_rp*byx(i,j-1,k))/128.0_rp - byy(i,j-1,k)/576.0_rp) + &
                    h**2*(-Gxy(i-1,j,k)/256.0_rp - Gxy(i,j,k)/2304.0_rp + Gyx(i,j-1,k)/288.0_rp)

                S(0,-2,0) = h*((-9.0_rp*bxy(i-1,j,k))/256.0_rp + (9.0_rp*bxy(i,j,k))/256.0_rp - &
                    (5.0_rp*byy(i,j-1,k))/144.0_rp - &
                    (9.0_rp*bzy(i,j,k-1))/256.0_rp + (9.0_rp*bzy(i,j,k))/256.0_rp) + &
                    h**2*(-Gxy(i-1,j,k)/256.0_rp + Gxy(i,j,k)/256.0_rp - &
                    Gyx(i,j-1,k)/128.0_rp - Gyz(i,j-1,k)/128.0_rp - Gzy(i,j,k-1)/256.0_rp + Gzy(i,j,k)/256.0_rp)

                S(1,-2,0) = h*(bxy(i-1,j,k)/256.0_rp + (9.0_rp*bxy(i,j,k))/256.0_rp + &
                    (5.0_rp*byx(i,j-1,k))/128.0_rp - &
                    byy(i,j-1,k)/576.0_rp) + h**2*(Gxy(i-1,j,k)/2304.0_rp + &
                    Gxy(i,j,k)/256.0_rp + Gyx(i,j-1,k)/288.0_rp)

                S(2,-2,0) = h*(-bxy(i,j,k)/256.0_rp - byx(i,j-1,k)/256.0_rp) + &
                    h**2*(-Gxy(i,j,k)/2304.0_rp + Gyx(i,j-1,k)/2304.0_rp)

                S(-2,0,1) = h*(-bxx(i-1,j,k)/576.0_rp + (5.0_rp*bxz(i-1,j,k))/128.0_rp + &
                    bzx(i,j,k-1)/256.0_rp + (9.0_rp*bzx(i,j,k))/256.0_rp) + &
                    h**2*(Gxz(i-1,j,k)/288.0_rp + Gzx(i,j,k-1)/2304.0_rp + Gzx(i,j,k)/256.0_rp)

                S(-1,0,1) = h*((3.0_rp*bxx(i-1,j,k))/64.0_rp + bxx(i,j,k)/576.0_rp - &
                    (45.0_rp*bxz(i-1,j,k))/128.0_rp - &
                    (5.0_rp*bxz(i,j,k))/128.0_rp - (5.0_rp*bzx(i,j,k-1))/128.0_rp - &
                    (45.0_rp*bzx(i,j,k))/128.0_rp + bzz(i,j,k-1)/576.0_rp + &
                    (3.0_rp*bzz(i,j,k))/64.0_rp) + h**2*(-Gxz(i-1,j,k)/32.0_rp - &
                    Gxz(i,j,k)/288.0_rp + Gzx(i,j,k-1)/288.0_rp + &
                    Gzx(i,j,k)/32.0_rp)

                S(1,0,1) = h*(bxx(i-1,j,k)/576.0_rp + (3.0_rp*bxx(i,j,k))/64.0_rp + &
                    (5.0_rp*bxz(i-1,j,k))/128.0_rp + (45.0_rp*bxz(i,j,k))/128.0_rp + &
                    (5.0_rp*bzx(i,j,k-1))/128.0_rp + (45.0_rp*bzx(i,j,k))/128.0_rp + &
                    bzz(i,j,k-1)/576.0_rp + (3.0_rp*bzz(i,j,k))/64.0_rp) + &
                    h**2*(Gxz(i-1,j,k)/288.0_rp + Gxz(i,j,k)/32.0_rp + &
                    Gzx(i,j,k-1)/288.0_rp + Gzx(i,j,k)/32.0_rp)

                S(2,0,1) = h*(-bxx(i,j,k)/576.0_rp - (5.0_rp*bxz(i,j,k))/128.0_rp - &
                    bzx(i,j,k-1)/256.0_rp - (9.0_rp*bzx(i,j,k))/256.0_rp) + &
                    h**2*(-Gxz(i,j,k)/288.0_rp + Gzx(i,j,k-1)/2304.0_rp + Gzx(i,j,k)/256.0_rp)

                S(-2,0,2) = h*(-bxz(i-1,j,k)/256.0_rp - bzx(i,j,k)/256.0_rp) + &
                    h**2*(Gxz(i-1,j,k)/2304.0_rp - Gzx(i,j,k)/2304.0_rp)

                S(-1,0,2) = h*((9.0_rp*bxz(i-1,j,k))/256.0_rp + bxz(i,j,k)/256.0_rp + &
                    (5.0_rp*bzx(i,j,k))/128.0_rp - bzz(i,j,k)/576.0_rp) + &
                    h**2*(-Gxz(i-1,j,k)/256.0_rp - Gxz(i,j,k)/2304.0_rp - Gzx(i,j,k)/288.0_rp)

                S(1,0,2) = h*(-bxz(i-1,j,k)/256.0_rp - (9.0_rp*bxz(i,j,k))/256.0_rp - &
                    (5.0_rp*bzx(i,j,k))/128.0_rp - bzz(i,j,k)/576.0_rp) + &
                    h**2*(Gxz(i-1,j,k)/2304.0_rp + Gxz(i,j,k)/256.0_rp - Gzx(i,j,k)/288.0_rp)

                S(2,0,2) = h*(bxz(i,j,k)/256.0_rp + bzx(i,j,k)/256.0_rp) + &
                    h**2*(-Gxz(i,j,k)/2304.0_rp - Gzx(i,j,k)/2304.0_rp)

                S(-2,0,-1) = h*(-bxx(i-1,j,k)/576.0_rp - (5.0_rp*bxz(i-1,j,k))/128.0_rp - &
                    (9.0_rp*bzx(i,j,k-1))/256.0_rp - bzx(i,j,k)/256.0_rp) + &
                    h**2*(Gxz(i-1,j,k)/288.0_rp - Gzx(i,j,k-1)/256.0_rp - Gzx(i,j,k)/2304.0_rp)

                S(-1,0,-1) = h*((3.0_rp*bxx(i-1,j,k))/64.0_rp + bxx(i,j,k)/576.0_rp + &
                    (45.0_rp*bxz(i-1,j,k))/128.0_rp + (5.0_rp*bxz(i,j,k))/128.0_rp + &
                    (45.0_rp*bzx(i,j,k-1))/128.0_rp + (5.0_rp*bzx(i,j,k))/128.0_rp + &
                    (3.0_rp*bzz(i,j,k-1))/64.0_rp + bzz(i,j,k)/576.0_rp) + &
                    h**2*(-Gxz(i-1,j,k)/32.0_rp - Gxz(i,j,k)/288.0_rp - &
                    Gzx(i,j,k-1)/32.0_rp - Gzx(i,j,k)/288.0_rp)

                S(0,0,-1) = h*((-3.0_rp*bxx(i-1,j,k))/64.0_rp - (3.0_rp*bxx(i,j,k))/64.0_rp + &
                    (45.0_rp*bxz(i-1,j,k))/128.0_rp - &
                    (45.0_rp*bxz(i,j,k))/128.0_rp - (3.0_rp*byy(i,j-1,k))/64.0_rp - &
                    (3.0_rp*byy(i,j,k))/64.0_rp + (45.0_rp*byz(i,j-1,k))/128.0_rp - &
                    (45.0_rp*byz(i,j,k))/128.0_rp + (15.0_rp*bzz(i,j,k-1))/16.0_rp + &
                    (5.0_rp*bzz(i,j,k))/144.0_rp) + h**2*(-Gxz(i-1,j,k)/32.0_rp + &
                    Gxz(i,j,k)/32.0_rp - Gyz(i,j-1,k)/32.0_rp + Gyz(i,j,k)/32.0_rp + &
                    (9.0_rp*Gzx(i,j,k-1))/128.0_rp + &
                    Gzx(i,j,k)/128.0_rp + (9.0_rp*Gzy(i,j,k-1))/128.0_rp + Gzy(i,j,k)/128.0_rp)

                S(1,0,-1) = h*(bxx(i-1,j,k)/576.0_rp + (3.0_rp*bxx(i,j,k))/64.0_rp - &
                    (5.0_rp*bxz(i-1,j,k))/128.0_rp - (45.0_rp*bxz(i,j,k))/128.0_rp - &
                    (45.0_rp*bzx(i,j,k-1))/128.0_rp - (5.0_rp*bzx(i,j,k))/128.0_rp + &
                    (3.0_rp*bzz(i,j,k-1))/64.0_rp + bzz(i,j,k)/576.0_rp) + &
                    h**2*(Gxz(i-1,j,k)/288.0_rp + Gxz(i,j,k)/32.0_rp - &
                    Gzx(i,j,k-1)/32.0_rp - Gzx(i,j,k)/288.0_rp)

                S(2,0,-1) = h*(-bxx(i,j,k)/576.0_rp + (5.0_rp*bxz(i,j,k))/128.0_rp + &
                    (9.0_rp*bzx(i,j,k-1))/256.0_rp + bzx(i,j,k)/256.0_rp) + &
                    h**2*(-Gxz(i,j,k)/288.0_rp - Gzx(i,j,k-1)/256.0_rp - Gzx(i,j,k)/2304.0_rp)

                S(-2,0,-2) = h*(bxz(i-1,j,k)/256.0_rp + bzx(i,j,k-1)/256.0_rp) + &
                    h**2*(Gxz(i-1,j,k)/2304.0_rp + Gzx(i,j,k-1)/2304.0_rp)

                S(-1,0,-2) = h*((-9.0_rp*bxz(i-1,j,k))/256.0_rp - bxz(i,j,k)/256.0_rp - &
                    (5.0_rp*bzx(i,j,k-1))/128.0_rp - bzz(i,j,k-1)/576.0_rp) + &
                    h**2*(-Gxz(i-1,j,k)/256.0_rp - Gxz(i,j,k)/2304.0_rp + Gzx(i,j,k-1)/288.0_rp)

                S(0,0,-2) = h*((-9.0_rp*bxz(i-1,j,k))/256.0_rp + (9.0_rp*bxz(i,j,k))/256.0_rp - &
                    (9.0_rp*byz(i,j-1,k))/256.0_rp + &
                    (9.0_rp*byz(i,j,k))/256.0_rp - (5.0_rp*bzz(i,j,k-1))/144.0_rp) + &
                    h**2*(-Gxz(i-1,j,k)/256.0_rp + Gxz(i,j,k)/256.0_rp - &
                    Gyz(i,j-1,k)/256.0_rp + Gyz(i,j,k)/256.0_rp - &
                    Gzx(i,j,k-1)/128.0_rp - Gzy(i,j,k-1)/128.0_rp)

                S(1,0,-2) = h*(bxz(i-1,j,k)/256.0_rp + (9.0_rp*bxz(i,j,k))/256.0_rp + &
                    (5.0_rp*bzx(i,j,k-1))/128.0_rp - bzz(i,j,k-1)/576.0_rp) + &
                    h**2*(Gxz(i-1,j,k)/2304.0_rp + Gxz(i,j,k)/256.0_rp + Gzx(i,j,k-1)/288.0_rp)

                S(2,0,-2) = h*(-bxz(i,j,k)/256.0_rp - bzx(i,j,k-1)/256.0_rp) + &
                    h**2*(-Gxz(i,j,k)/2304.0_rp + Gzx(i,j,k-1)/2304.0_rp)

                S(0,-2,1) = h*(-byy(i,j-1,k)/576.0_rp + (5.0_rp*byz(i,j-1,k))/128.0_rp + &
                    bzy(i,j,k-1)/256.0_rp + &
                    (9.0_rp*bzy(i,j,k))/256.0_rp) + h**2*(Gyz(i,j-1,k)/288.0_rp + &
                    Gzy(i,j,k-1)/2304.0_rp + Gzy(i,j,k)/256.0_rp)

                S(0,-1,1) = h*((3.0_rp*byy(i,j-1,k))/64.0_rp + byy(i,j,k)/576.0_rp - &
                    (45.0_rp*byz(i,j-1,k))/128.0_rp - (5.0_rp*byz(i,j,k))/128.0_rp - &
                    (5.0_rp*bzy(i,j,k-1))/128.0_rp - (45.0_rp*bzy(i,j,k))/128.0_rp + &
                    bzz(i,j,k-1)/576.0_rp + (3.0_rp*bzz(i,j,k))/64.0_rp) + &
                    h**2*(-Gyz(i,j-1,k)/32.0_rp - Gyz(i,j,k)/288.0_rp + &
                    Gzy(i,j,k-1)/288.0_rp + Gzy(i,j,k)/32.0_rp)

                S(0,0,1) = h*((-3.0_rp*bxx(i-1,j,k))/64.0_rp - (3.0_rp*bxx(i,j,k))/64.0_rp - &
                    (45.0_rp*bxz(i-1,j,k))/128.0_rp + (45.0_rp*bxz(i,j,k))/128.0_rp - &
                    (3.0_rp*byy(i,j-1,k))/64.0_rp - (3.0_rp*byy(i,j,k))/64.0_rp - &
                    (45.0_rp*byz(i,j-1,k))/128.0_rp + (45.0_rp*byz(i,j,k))/128.0_rp + &
                    (5.0_rp*bzz(i,j,k-1))/144.0_rp + (15.0_rp*bzz(i,j,k))/16.0_rp) + &
                    h**2*(-Gxz(i-1,j,k)/32.0_rp + Gxz(i,j,k)/32.0_rp - &
                    Gyz(i,j-1,k)/32.0_rp + Gyz(i,j,k)/32.0_rp - Gzx(i,j,k-1)/128.0_rp - &
                    (9.0_rp*Gzx(i,j,k))/128.0_rp - Gzy(i,j,k-1)/128.0_rp - &
                    (9.0_rp*Gzy(i,j,k))/128.0_rp)

                S(0,1,1) = h*(byy(i,j-1,k)/576.0_rp + (3.0_rp*byy(i,j,k))/64.0_rp + &
                    (5.0_rp*byz(i,j-1,k))/128.0_rp + (45.0_rp*byz(i,j,k))/128.0_rp + &
                    (5.0_rp*bzy(i,j,k-1))/128.0_rp + (45.0_rp*bzy(i,j,k))/128.0_rp + &
                    bzz(i,j,k-1)/576.0_rp + (3.0_rp*bzz(i,j,k))/64.0_rp) + &
                    h**2*(Gyz(i,j-1,k)/288.0_rp + Gyz(i,j,k)/32.0_rp + &
                    Gzy(i,j,k-1)/288.0_rp + Gzy(i,j,k)/32.0_rp)

                S(0,2,1) = h*(-byy(i,j,k)/576.0_rp - (5.0_rp*byz(i,j,k))/128.0_rp - &
                    bzy(i,j,k-1)/256.0_rp - (9.0_rp*bzy(i,j,k))/256.0_rp) + &
                    h**2*(-Gyz(i,j,k)/288.0_rp + Gzy(i,j,k-1)/2304.0_rp + Gzy(i,j,k)/256.0_rp)

                S(0,-2,2) = h*(-byz(i,j-1,k)/256.0_rp - bzy(i,j,k)/256.0_rp) + &
                    h**2*(Gyz(i,j-1,k)/2304.0_rp - Gzy(i,j,k)/2304.0_rp)

                S(0,-1,2) = h*((9.0_rp*byz(i,j-1,k))/256.0_rp + byz(i,j,k)/256.0_rp + &
                    (5.0_rp*bzy(i,j,k))/128.0_rp - bzz(i,j,k)/576.0_rp) + &
                    h**2*(-Gyz(i,j-1,k)/256.0_rp - Gyz(i,j,k)/2304.0_rp - Gzy(i,j,k)/288.0_rp)

                S(0,0,2) = h*((9.0_rp*bxz(i-1,j,k))/256.0_rp - (9.0_rp*bxz(i,j,k))/256.0_rp + &
                    (9.0_rp*byz(i,j-1,k))/256.0_rp - (9.0_rp*byz(i,j,k))/256.0_rp - &
                    (5.0_rp*bzz(i,j,k))/144.0_rp) + h**2*(-Gxz(i-1,j,k)/256.0_rp + &
                    Gxz(i,j,k)/256.0_rp - Gyz(i,j-1,k)/256.0_rp + &
                    Gyz(i,j,k)/256.0_rp + Gzx(i,j,k)/128.0_rp + Gzy(i,j,k)/128.0_rp)

                S(0,1,2) = h*(-byz(i,j-1,k)/256.0_rp - (9.0_rp*byz(i,j,k))/256.0_rp - &
                    (5.0_rp*bzy(i,j,k))/128.0_rp - bzz(i,j,k)/576.0_rp) + &
                    h**2*(Gyz(i,j-1,k)/2304.0_rp + Gyz(i,j,k)/256.0_rp - Gzy(i,j,k)/288.0_rp)

                S(0,2,2) = h*(byz(i,j,k)/256.0_rp + bzy(i,j,k)/256.0_rp) + &
                    h**2*(-Gyz(i,j,k)/2304.0_rp - Gzy(i,j,k)/2304.0_rp)

                S(0,-2,-1) = h*(-byy(i,j-1,k)/576.0_rp - (5.0_rp*byz(i,j-1,k))/128.0_rp - &
                    (9.0_rp*bzy(i,j,k-1))/256.0_rp - bzy(i,j,k)/256.0_rp) + &
                    h**2*(Gyz(i,j-1,k)/288.0_rp - Gzy(i,j,k-1)/256.0_rp - Gzy(i,j,k)/2304.0_rp)

                S(0,-1,-1) = h*((3.0_rp*byy(i,j-1,k))/64.0_rp + byy(i,j,k)/576.0_rp + &
                    (45.0_rp*byz(i,j-1,k))/128.0_rp + (5.0_rp*byz(i,j,k))/128.0_rp + &
                    (45.0_rp*bzy(i,j,k-1))/128.0_rp + (5.0_rp*bzy(i,j,k))/128.0_rp + &
                    (3.0_rp*bzz(i,j,k-1))/64.0_rp + bzz(i,j,k)/576.0_rp) + &
                    h**2*(-Gyz(i,j-1,k)/32.0_rp - Gyz(i,j,k)/288.0_rp - &
                    Gzy(i,j,k-1)/32.0_rp - Gzy(i,j,k)/288.0_rp)

                S(0,1,-1) = h*(byy(i,j-1,k)/576.0_rp + (3.0_rp*byy(i,j,k))/64.0_rp - &
                    (5.0_rp*byz(i,j-1,k))/128.0_rp - (45.0_rp*byz(i,j,k))/128.0_rp - &
                    (45.0_rp*bzy(i,j,k-1))/128.0_rp - (5.0_rp*bzy(i,j,k))/128.0_rp + &
                    (3.0_rp*bzz(i,j,k-1))/64.0_rp + bzz(i,j,k)/576.0_rp) + &
                    h**2*(Gyz(i,j-1,k)/288.0_rp + Gyz(i,j,k)/32.0_rp - &
                    Gzy(i,j,k-1)/32.0_rp - Gzy(i,j,k)/288.0_rp)

                S(0,2,-1) = h*(-byy(i,j,k)/576.0_rp + (5.0_rp*byz(i,j,k))/128.0_rp + &
                    (9.0_rp*bzy(i,j,k-1))/256.0_rp + bzy(i,j,k)/256.0_rp) + &
                    h**2*(-Gyz(i,j,k)/288.0_rp - Gzy(i,j,k-1)/256.0_rp - Gzy(i,j,k)/2304.0_rp)

                S(0,-2,-2) = h*(byz(i,j-1,k)/256.0_rp + bzy(i,j,k-1)/256.0_rp) + &
                    h**2*(Gyz(i,j-1,k)/2304.0_rp + Gzy(i,j,k-1)/2304.0_rp)

                S(0,-1,-2) = h*((-9.0_rp*byz(i,j-1,k))/256.0_rp - byz(i,j,k)/256.0_rp - &
                    (5.0_rp*bzy(i,j,k-1))/128.0_rp - bzz(i,j,k-1)/576.0_rp) + &
                    h**2*(-Gyz(i,j-1,k)/256.0_rp - Gyz(i,j,k)/2304.0_rp + Gzy(i,j,k-1)/288.0_rp)

                S(0,1,-2) = h*(byz(i,j-1,k)/256.0_rp + (9.0_rp*byz(i,j,k))/256.0_rp + &
                    (5.0_rp*bzy(i,j,k-1))/128.0_rp - bzz(i,j,k-1)/576.0_rp) + &
                    h**2*(Gyz(i,j-1,k)/2304.0_rp + Gyz(i,j,k)/256.0_rp + Gzy(i,j,k-1)/288.0_rp)

                S(0,2,-2) = h*(-byz(i,j,k)/256.0_rp - bzy(i,j,k-1)/256.0_rp) + &
                    h**2*(-Gyz(i,j,k)/2304.0_rp + Gzy(i,j,k-1)/2304.0_rp)

                Stencil(row,0) = 0
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            Stencil(row,stencilindex(nghi,nghj,nghk)) = S(nghi,nghj,nghk)
                        end do
                    end do
                end do

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    if(allocated(S)) deallocate(S)

end subroutine ComputeStencil4
