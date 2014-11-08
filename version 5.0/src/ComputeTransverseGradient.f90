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

    if(rank .eq. master) then
        print *, ''
        print *, '********************************'
        print *, 'Calculating Transverse Gradient'
        print *, '********************************'
        print *, ''
    endif

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin-1, imax

                Gxy(i,j,k) = (1.0_rp/(2.0_rp*h))*(bxy(i,j+1,k)-bxy(i,j-1,k))
                Gxz(i,j,k) = (1.0_rp/(2.0_rp*h))*(bxz(i,j,k+1)-bxz(i,j,k-1))

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin-1, jmax
            do i = imin, imax

                Gyx(i,j,k) = (1.0_rp/(2.0_rp*h))*(byx(i+1,j,k)-byx(i-1,j,k))
                Gyz(i,j,k) = (1.0_rp/(2.0_rp*h))*(byz(i,j,k+1)-byz(i,j,k-1))

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax
        do j = jmin, jmax
            do i = imin, imax

                Gzx(i,j,k) = (1.0_rp/(2.0_rp*h))*(bzx(i+1,j,k)-bzx(i-1,j,k))
                Gzy(i,j,k) = (1.0_rp/(2.0_rp*h))*(bzy(i,j+1,k)-bzy(i,j-1,k))

            end do
        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL


end subroutine ComputeTransverseGradient
