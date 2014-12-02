 !****************************************************************  
 ! IC_Ploytrop1.f90
 ! 
 ! Created on: Nov 16, 2014
 ! Author: ninoy
 !***************************************************************** 
subroutine IC_Ploytrop1
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k,counter
    real(kind=double) :: d,buffer,alpha
    integer :: ierr

    print *, ''
    print *, '***********************************'
    print *, 'Initializing Scenerio: Polytrop n=1'
    print *, '***********************************'
    print *, ''

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
    alpha = 0.2_rp
    R = pi*alpha

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,d)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                d = sqrt(x(i)**2 + y(j)**2 + z(k)**2)
                rho(i,j,k) = sin(d/alpha)/(d/alpha)
                if(d .eq. 0.0_rp) rho(i,j,k) = 1.0_rp
                if(rho(i,j,k) .lt. 0.0_rp) rho(i,j,k) = 1.0d-17
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    M = M*h**3

    print *, 'M=', M
    print *, 'R=', R
    print *, 'alpha=', alpha


end subroutine IC_Ploytrop1
