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
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin-1, imax

                bxx(i,j,k) = 1.0_rp
                bxy(i,j,k) = 0.0_rp
                bxz(i,j,k) = 0.0_rp

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

                byx(i,j,k) = 0.0_rp
                byy(i,j,k) = 1.0_rp
                byz(i,j,k) = 0.0_rp

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

                bzx(i,j,k) = 0.0_rp
                bzy(i,j,k) = 0.0_rp
                bzz(i,j,k) = 1.0_rp

            end do
        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin-1, imax

                if(j .eq. jmin) then
                    Gxy(i,j,k) = -(1.0_rp/(2.0_rp*h))*(bxy(i,j+2,k)-4.0_rp*bxy(i,j+1,k)+3.0_rp*bxy(i,j,k))
                else if(j .eq. jmax) then
                    Gxy(i,j,k) = (1.0_rp/(2.0_rp*h))*(bxy(i,j-2,k)-4.0_rp*bxy(i,j-1,k)+3.0_rp*bxy(i,j,k))
                else
                    Gxy(i,j,k) = (1.0_rp/(2.0_rp*h))*(bxy(i,j+1,k)-bxy(i,j-1,k))
                end if

                if(k .eq. kmin) then
                    Gxz(i,j,k) = -(1.0_rp/(2.0_rp*h))*(bxz(i,j,k+2)-4.0_rp*bxz(i,j,k+1)+3.0_rp*bxz(i,j,k))
                else if(k .eq. kmax) then
                    Gxz(i,j,k) = (1.0_rp/(2.0_rp*h))*(bxz(i,j,k-2)-4.0_rp*bxz(i,j,k-1)+3.0_rp*bxz(i,j,k))
                else
                    Gxz(i,j,k) = (1.0_rp/(2.0_rp*h))*(bxz(i,j,k+1)-bxz(i,j,k-1))
                end if

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

                if(i .eq. imin) then
                    Gyx(i,j,k) = -(1.0_rp/(2.0_rp*h))*(byx(i+2,j,k)-4.0_rp*byx(i+1,j,k)+3.0_rp*byx(i,j,k))
                else if(i .eq. imax) then
                    Gyx(i,j,k) = (1.0_rp/(2.0_rp*h))*(byx(i-2,j,k)-4.0_rp*byx(i-1,j,k)+3.0_rp*byx(i,j,k))
                else
                    Gyx(i,j,k) = (1.0_rp/(2.0_rp*h))*(byx(i+1,j,k)-byx(i-1,j,k))
                end if

                if(k .eq. kmin) then
                    Gyz(i,j,k) = -(1.0_rp/(2.0_rp*h))*(byz(i,j,k+2)-4.0_rp*byz(i,j,k+1)+3.0_rp*byz(i,j,k))
                else if(k .eq. kmax) then
                    Gyz(i,j,k) = (1.0_rp/(2.0_rp*h))*(byz(i,j,k-2)-4.0_rp*byz(i,j,k-1)+3.0_rp*byz(i,j,k))
                else
                    Gyz(i,j,k) = (1.0_rp/(2.0_rp*h))*(byz(i,j,k+1)-byz(i,j,k-1))
                end if

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

                if(i .eq. imin) then
                    Gzx(i,j,k) = -(1.0_rp/(2.0_rp*h))*(bzx(i+2,j,k)-4.0_rp*bzx(i+1,j,k)+3.0_rp*bzx(i,j,k))
                else if(i .eq. imax) then
                    Gzx(i,j,k) = (1.0_rp/(2.0_rp*h))*(bzx(i-2,j,k)-4.0_rp*bzx(i-1,j,k)+3.0_rp*bzx(i,j,k))
                else
                    Gzx(i,j,k) = (1.0_rp/(2.0_rp*h))*(bzx(i+1,j,k)-bzx(i-1,j,k))
                end if

                if(j .eq. jmin) then
                    Gzy(i,j,k) = -(1.0_rp/(2.0_rp*h))*(bzy(i,j+2,k)-4.0_rp*bzy(i,j+1,k)+3.0_rp*bzy(i,j,k))
                else if(j .eq. jmax) then
                    Gzy(i,j,k) = (1.0_rp/(2.0_rp*h))*(bzy(i,j-2,k)-4.0_rp*bzy(i,j-1,k)+3.0_rp*bzy(i,j,k))
                else
                    Gzy(i,j,k) = (1.0_rp/(2.0_rp*h))*(bzy(i,j+1,k)-bzy(i,j-1,k))
                end if

            end do
        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

end subroutine ComputeTransverseGradient
