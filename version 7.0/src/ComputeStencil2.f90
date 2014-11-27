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

                S(0,0,0) = -162.0_rp*h/24.0_rp

                S(0,0,-2) = -1.0_rp*h/24.0_rp
                S(0,0,-1) = +28.0_rp*h/24.0_rp
                S(0,0,+1) = +28.0_rp*h/24.0_rp
                S(0,0,+2) = -1.0_rp*h/24.0_rp

                S(0,-2,0) = -1.0_rp*h/24.0_rp
                S(0,-1,0) = +28.0_rp*h/24.0_rp
                S(0,+1,0) = +28.0_rp*h/24.0_rp
                S(0,+2,0) = -1.0_rp*h/24.0_rp

                S(-2,0,0) = -1.0_rp*h/24.0_rp
                S(-1,0,0) = +28.0_rp*h/24.0_rp
                S(+1,0,0) = +28.0_rp*h/24.0_rp
                S(+2,0,0) = -1.0_rp*h/24.0_rp


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
