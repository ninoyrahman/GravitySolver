 !****************************************************************  
 ! CalculateStencil2.f90
 ! 
 ! Created on: Sep 21, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine CalculateStencil
    use MD_Parameter
    use MD_Quantity
    use MD_GeometricQuantity
    use MD_Helper
    implicit none

    integer :: i,j,k, row
    integer :: nghi,nghj,nghk

    print *, ''
    print *, '********************'
    print *, 'Calculating Stencil'
    print *, '********************'
    print *, ''

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax

                row = GetIndex(i,j,k)

                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            S(row,nghi,nghj,nghk) = 0.0_rp
                        end do
                    end do
                end do

                S(row, 0,0,0) = -162.0_rp*h/24.0_rp

                S(row, 0,0,-2) = -1.0_rp*h/24.0_rp
                S(row, 0,0,-1) = +28.0_rp*h/24.0_rp
                S(row, 0,0,+1) = +28.0_rp*h/24.0_rp
                S(row, 0,0,+2) = -1.0_rp*h/24.0_rp

                S(row, 0,-2,0) = -1.0_rp*h/24.0_rp
                S(row, 0,-1,0) = +28.0_rp*h/24.0_rp
                S(row, 0,+1,0) = +28.0_rp*h/24.0_rp
                S(row, 0,+2,0) = -1.0_rp*h/24.0_rp

                S(row, -2,0,0) = -1.0_rp*h/24.0_rp
                S(row, -1,0,0) = +28.0_rp*h/24.0_rp
                S(row, +1,0,0) = +28.0_rp*h/24.0_rp
                S(row, +2,0,0) = -1.0_rp*h/24.0_rp

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

end subroutine CalculateStencil
