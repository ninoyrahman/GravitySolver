 !****************************************************************  
 ! ComputeCartesianMapping.f90
 ! 
 ! Created on: Oct 31, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine ComputeCartesianMapping
    use MD_Parameter
    use MD_GeometricQuantity
    implicit none

    integer :: i,j,k

    if(rank .eq. master) then
        print *, ''
        print *, '******************************'
        print *, 'Calculating Cartesian Mapping'
        print *, '******************************'
        print *, ''
    endif


    ! looping over inner cell and calculating bxx
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin-1, imax
                bxx(i,j,k) = 1.0_rp
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    ! looping over inner cell and calculating bxy
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin-1, jmax+1
            do i = imin-1, imax
                bxy(i,j,k) = 0.0_rp
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    ! looping over inner cell and calculating bxz
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin, jmax
            do i = imin-1, imax
                bxz(i,j,k) = 0.0_rp
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    ! looping over inner cell and calculating byx
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin-1, jmax
            do i = imin-1, imax+1
                byx(i,j,k) = 0.0_rp
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell and calculating byy
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin-1, jmax
            do i = imin, imax
                byy(i,j,k) = 1.0_rp
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell and calculating byz
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax
            do i = imin, imax
                byz(i,j,k) = 0.0_rp
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell and calculating bzx
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax
        do j = jmin, jmax
            do i = imin-1, imax+1
                bzx(i,j,k) = 0.0_rp
            end do
        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    ! looping over inner cell and calculating bzy
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax
        do j = jmin-1, jmax+1
            do i = imin, imax
                bzy(i,j,k) = 0.0_rp
            end do
        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    ! looping over inner cell and calculating bzz
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax
        do j = jmin, jmax
            do i = imin, imax
                bzz(i,j,k) = 1.0_rp
            end do
        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

end subroutine ComputeCartesianMapping
