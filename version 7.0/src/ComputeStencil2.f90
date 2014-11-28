 !****************************************************************  
 ! ComputeStencil2.f90
 ! 
 ! Created on: Sep 21, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine ComputeStencil2
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
        print *, 'Calculating Stencil 2nd-Order'
        print *, '*****************************'
        print *, ''
    end if

    if(.not. allocated(S)) allocate(S(-swidth:swidth,-swidth:swidth,-swidth:swidth))

    counter = 0
    do nghk = -swidth, +swidth
        do nghj = -swidth, +swidth
            do nghi = -swidth, +swidth
                stencilindex(nghi,nghj,nghk) = 0
                if((nghi .eq. 0 .and. nghk .eq. 0) .or. (nghj .eq. 0 .and. nghk .eq. 0) &
                    .or. (nghi .eq. 0 .and. nghj .eq. 0)) then
                    counter = counter + 1
                    stencilindex(nghi,nghj,nghk) = counter
                end if
            end do
        end do
    end do

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

                S(-2,0,0) = h*(-bxx(-0.5 + i,j,k)/24. - (3*byx(i,-0.5 + j,k))/64. + &
                    (3*byx(i,0.5 + j,k))/64. - (3*bzx(i,j,-0.5 + k))/64. + (3*bzx(i,j,0.5 + k))/64.)

                S(-1,0,0) = h*((9*bxx(-0.5 + i,j,k))/8. + bxx(0.5 + i,j,k)/24. + (3*byx(i,-0.5 + j,k))/8. - &
                    (3*byx(i,0.5 + j,k))/8. + (3*bzx(i,j,-0.5 + k))/8. - (3*bzx(i,j,0.5 + k))/8.)

                S(0,0,0) = h*((-9*bxx(-0.5 + i,j,k))/8. - (9*bxx(0.5 + i,j,k))/8. - (9*byy(i,-0.5 + j,k))/8. - &
                    (9*byy(i,0.5 + j,k))/8. - (9*bzz(i,j,-0.5 + k))/8. - (9*bzz(i,j,0.5 + k))/8.)

                S(1,0,0) = h*(bxx(-0.5 + i,j,k)/24. + (9*bxx(0.5 + i,j,k))/8. - (3*byx(i,-0.5 + j,k))/8. + &
                    (3*byx(i,0.5 + j,k))/8. - (3*bzx(i,j,-0.5 + k))/8. + (3*bzx(i,j,0.5 + k))/8.)

                S(2,0,0) = h*(-bxx(0.5 + i,j,k)/24. + (3*byx(i,-0.5 + j,k))/64. - &
                    (3*byx(i,0.5 + j,k))/64. + (3*bzx(i,j,-0.5 + k))/64. - (3*bzx(i,j,0.5 + k))/64.)

                S(-2,1,0) = h*(bxy(-0.5 + i,j,k)/24. + byx(i,-0.5 + j,k)/192. + (3*byx(i,0.5 + j,k))/64.)

                S(-1,1,0) = h*((-3*bxy(-0.5 + i,j,k))/8. - bxy(0.5 + i,j,k)/24. - &
                    byx(i,-0.5 + j,k)/24. - (3*byx(i,0.5 + j,k))/8.)

                S(0,1,0) = h*((-3*bxy(-0.5 + i,j,k))/8. + (3*bxy(0.5 + i,j,k))/8. + byy(i,-0.5 + j,k)/24. + &
                    (9*byy(i,0.5 + j,k))/8. - (3*bzy(i,j,-0.5 + k))/8. + (3*bzy(i,j,0.5 + k))/8.)

                S(1,1,0) = h*(bxy(-0.5 + i,j,k)/24. + (3*bxy(0.5 + i,j,k))/8. + &
                    byx(i,-0.5 + j,k)/24. + (3*byx(i,0.5 + j,k))/8.)

                S(2,1,0) = h*(-bxy(0.5 + i,j,k)/24. - byx(i,-0.5 + j,k)/192. - &
                    (3*byx(i,0.5 + j,k))/64.)

                S(-2,2,0) = h*(-bxy(-0.5 + i,j,k)/192. - byx(i,0.5 + j,k)/192.)

                S(-1,2,0) = h*((3*bxy(-0.5 + i,j,k))/64. + bxy(0.5 + i,j,k)/192. + byx(i,0.5 + j,k)/24.)

                S(0,2,0) = h*((3*bxy(-0.5 + i,j,k))/64. - (3*bxy(0.5 + i,j,k))/64. - &
                    byy(i,0.5 + j,k)/24. + (3*bzy(i,j,-0.5 + k))/64. - (3*bzy(i,j,0.5 + k))/64.)

                S(1,2,0) = h*(-bxy(-0.5 + i,j,k)/192. - (3*bxy(0.5 + i,j,k))/64. - &
                    byx(i,0.5 + j,k)/24.)

                S(2,2,0) = h*(bxy(0.5 + i,j,k)/192. + byx(i,0.5 + j,k)/192.)

                S(-2,-1,0) = h*(-bxy(-0.5 + i,j,k)/24. - (3*byx(i,-0.5 + j,k))/64. - byx(i,0.5 + j,k)/192.)

                S(-1,-1,0) = h*((3*bxy(-0.5 + i,j,k))/8. + bxy(0.5 + i,j,k)/24. + &
                    (3*byx(i,-0.5 + j,k))/8. + byx(i,0.5 + j,k)/24.)

                S(0,-1,0) = h*((3*bxy(-0.5 + i,j,k))/8. - (3*bxy(0.5 + i,j,k))/8. + (9*byy(i,-0.5 + j,k))/8. + &
                    byy(i,0.5 + j,k)/24. + (3*bzy(i,j,-0.5 + k))/8. - (3*bzy(i,j,0.5 + k))/8.)

                S(1,-1,0) = h*(-bxy(-0.5 + i,j,k)/24. - (3*bxy(0.5 + i,j,k))/8. - &
                    (3*byx(i,-0.5 + j,k))/8. - byx(i,0.5 + j,k)/24.)

                S(2,-1,0) = h*(bxy(0.5 + i,j,k)/24. + (3*byx(i,-0.5 + j,k))/64. + byx(i,0.5 + j,k)/192.)

                S(-2,-2,0) = h*(bxy(-0.5 + i,j,k)/192. + byx(i,-0.5 + j,k)/192.)

                S(-1,-2,0) = h*((-3*bxy(-0.5 + i,j,k))/64. - bxy(0.5 + i,j,k)/192. - byx(i,-0.5 + j,k)/24.)

                S(0,-2,0) = h*((-3*bxy(-0.5 + i,j,k))/64. + (3*bxy(0.5 + i,j,k))/64. - &
                    byy(i,-0.5 + j,k)/24. - (3*bzy(i,j,-0.5 + k))/64. + (3*bzy(i,j,0.5 + k))/64.)

                S(1,-2,0) = h*(bxy(-0.5 + i,j,k)/192. + (3*bxy(0.5 + i,j,k))/64. + byx(i,-0.5 + j,k)/24.)

                S(2,-2,0) = h*(-bxy(0.5 + i,j,k)/192. - byx(i,-0.5 + j,k)/192.)

                S(-2,0,1) = h*(bxz(-0.5 + i,j,k)/24. + bzx(i,j,-0.5 + k)/192. + (3*bzx(i,j,0.5 + k))/64.)

                S(-1,0,1) = h*((-3*bxz(-0.5 + i,j,k))/8. - bxz(0.5 + i,j,k)/24. - &
                    bzx(i,j,-0.5 + k)/24. - (3*bzx(i,j,0.5 + k))/8.)

                S(1,0,1) = h*(bxz(-0.5 + i,j,k)/24. + (3*bxz(0.5 + i,j,k))/8. + &
                    bzx(i,j,-0.5 + k)/24. + (3*bzx(i,j,0.5 + k))/8.)

                S(2,0,1) = h*(-bxz(0.5 + i,j,k)/24. - bzx(i,j,-0.5 + k)/192. - (3*bzx(i,j,0.5 + k))/64.)

                S(-2,0,2) = h*(-bxz(-0.5 + i,j,k)/192. - bzx(i,j,0.5 + k)/192.)

                S(-1,0,2) = h*((3*bxz(-0.5 + i,j,k))/64. + bxz(0.5 + i,j,k)/192. + bzx(i,j,0.5 + k)/24.)

                S(1,0,2) = h*(-bxz(-0.5 + i,j,k)/192. - (3*bxz(0.5 + i,j,k))/64. - bzx(i,j,0.5 + k)/24.)

                S(2,0,2) = h*(bxz(0.5 + i,j,k)/192. + bzx(i,j,0.5 + k)/192.)

                S(-2,0,-1) = h*(-bxz(-0.5 + i,j,k)/24. - (3*bzx(i,j,-0.5 + k))/64. - bzx(i,j,0.5 + k)/192.)

                S(-1,0,-1) = h*((3*bxz(-0.5 + i,j,k))/8. + bxz(0.5 + i,j,k)/24. + &
                    (3*bzx(i,j,-0.5 + k))/8. + bzx(i,j,0.5 + k)/24.)

                S(0,0,-1) = h*((3*bxz(-0.5 + i,j,k))/8. - (3*bxz(0.5 + i,j,k))/8. + (3*byz(i,-0.5 + j,k))/8. - &
                    (3*byz(i,0.5 + j,k))/8. + (9*bzz(i,j,-0.5 + k))/8. + bzz(i,j,0.5 + k)/24.)

                S(1,0,-1) = h*(-bxz(-0.5 + i,j,k)/24. - (3*bxz(0.5 + i,j,k))/8. - &
                    (3*bzx(i,j,-0.5 + k))/8. - bzx(i,j,0.5 + k)/24.)

                S(2,0,-1) = h*(bxz(0.5 + i,j,k)/24. + (3*bzx(i,j,-0.5 + k))/64. + bzx(i,j,0.5 + k)/192.)

                S(-2,0,-2) = h*(bxz(-0.5 + i,j,k)/192. + bzx(i,j,-0.5 + k)/192.)

                S(-1,0,-2) = h*((-3*bxz(-0.5 + i,j,k))/64. - bxz(0.5 + i,j,k)/192. - bzx(i,j,-0.5 + k)/24.)

                S(0,0,-2) = h*((-3*bxz(-0.5 + i,j,k))/64. + (3*bxz(0.5 + i,j,k))/64. - &
                    (3*byz(i,-0.5 + j,k))/64. + (3*byz(i,0.5 + j,k))/64. - bzz(i,j,-0.5 + k)/24.)

                S(1,0,-2) = h*(bxz(-0.5 + i,j,k)/192. + (3*bxz(0.5 + i,j,k))/64. + bzx(i,j,-0.5 + k)/24.)

                S(2,0,-2) = h*(-bxz(0.5 + i,j,k)/192. - bzx(i,j,-0.5 + k)/192.)

                S(0,-2,1) = h*(byz(i,-0.5 + j,k)/24. + bzy(i,j,-0.5 + k)/192. + (3*bzy(i,j,0.5 + k))/64.)

                S(0,-1,1) = h*((-3*byz(i,-0.5 + j,k))/8. - byz(i,0.5 + j,k)/24. - &
                    bzy(i,j,-0.5 + k)/24. - (3*bzy(i,j,0.5 + k))/8.)

                S(0,0,1) = h*((-3*bxz(-0.5 + i,j,k))/8. + (3*bxz(0.5 + i,j,k))/8. - (3*byz(i,-0.5 + j,k))/8. + &
                    (3*byz(i,0.5 + j,k))/8. + bzz(i,j,-0.5 + k)/24. + (9*bzz(i,j,0.5 + k))/8.)

                S(0,1,1) = h*(byz(i,-0.5 + j,k)/24. + (3*byz(i,0.5 + j,k))/8. + bzy(i,j,-0.5 + k)/24. + &
                    (3*bzy(i,j,0.5 + k))/8.)

                S(0,2,1) = h*(-byz(i,0.5 + j,k)/24. - bzy(i,j,-0.5 + k)/192. - (3*bzy(i,j,0.5 + k))/64.)

                S(0,-2,2) = h*(-byz(i,-0.5 + j,k)/192. - bzy(i,j,0.5 + k)/192.)

                S(0,-1,2) = h*((3*byz(i,-0.5 + j,k))/64. + byz(i,0.5 + j,k)/192. + bzy(i,j,0.5 + k)/24.)

                S(0,0,2) = h*((3*bxz(-0.5 + i,j,k))/64. - (3*bxz(0.5 + i,j,k))/64. + (3*byz(i,-0.5 + j,k))/64. - &
                    (3*byz(i,0.5 + j,k))/64. - bzz(i,j,0.5 + k)/24.)

                S(0,1,2) = h*(-byz(i,-0.5 + j,k)/192. - (3*byz(i,0.5 + j,k))/64. - bzy(i,j,0.5 + k)/24.)

                S(0,2,2) = h*(byz(i,0.5 + j,k)/192. + bzy(i,j,0.5 + k)/192.)

                S(0,-2,-1) = h*(-byz(i,-0.5 + j,k)/24. - (3*bzy(i,j,-0.5 + k))/64. - bzy(i,j,0.5 + k)/192.)

                S(0,-1,-1) = h*((3*byz(i,-0.5 + j,k))/8. + byz(i,0.5 + j,k)/24. + (3*bzy(i,j,-0.5 + k))/8. + &
                    bzy(i,j,0.5 + k)/24.)

                S(0,1,-1) = h*(-byz(i,-0.5 + j,k)/24. - (3*byz(i,0.5 + j,k))/8. - (3*bzy(i,j,-0.5 + k))/8. - &
                    bzy(i,j,0.5 + k)/24.)

                S(0,2,-1) = h*(byz(i,0.5 + j,k)/24. + (3*bzy(i,j,-0.5 + k))/64. + bzy(i,j,0.5 + k)/192.)

                S(0,-2,-2) = h*(byz(i,-0.5 + j,k)/192. + bzy(i,j,-0.5 + k)/192.)

                S(0,-1,-2) = h*((-3*byz(i,-0.5 + j,k))/64. - byz(i,0.5 + j,k)/192. - bzy(i,j,-0.5 + k)/24.)

                S(0,1,-2) = h*(byz(i,-0.5 + j,k)/192. + (3*byz(i,0.5 + j,k))/64. + bzy(i,j,-0.5 + k)/24.)

                S(0,2,-2) = h*(-byz(i,0.5 + j,k)/192. - bzy(i,j,-0.5 + k)/192.)

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

end subroutine ComputeStencil2
