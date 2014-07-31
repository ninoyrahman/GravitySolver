 !****************************************************************
 ! ComputeTransverseGradient.f90
 !
 ! Created on: Jul 28, 2014
 ! Author: ninoy
 !*****************************************************************

 subroutine ComputeTransverseGradient
    use MD_Parameter
    use MD_GeometricQuantity
    implicit none

    integer :: i,j,k

    print *, ''
    print *, '********************************'
    print *, 'Calculating Transverse Gradient'
    print *, '********************************'
    print *, ''


    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do i = imin-1, imax
        do j = jmin-1, jmax
            do k = kmin-1, kmax

                bxx(i,j,k) = 1.0
                bxy(i,j,k) = 0.0
                bxz(i,j,k) = 0.0
                byx(i,j,k) = 0.0
                byy(i,j,k) = 1.0
                byz(i,j,k) = 0.0
                bzx(i,j,k) = 0.0
                bzy(i,j,k) = 0.0
                bzz(i,j,k) = 1.0

                Gxy(i,j,k) = 0.0
                Gxz(i,j,k) = 0.0
                Gyx(i,j,k) = 0.0
                Gyz(i,j,k) = 0.0
                Gzx(i,j,k) = 0.0
                Gzy(i,j,k) = 0.0

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

 end subroutine ComputeTransverseGradient
