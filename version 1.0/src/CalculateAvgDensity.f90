 !****************************************************************  
 ! CalculateAvgDensity.f90
 ! 
 ! Created on: Aug 13, 2014
 ! Author: ninoy
 !***************************************************************** 

 subroutine CalculateAvgDensity
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k
    real(kind=double) :: Lx,Ly,Lz,Mtmp

    print *, ''
    print *, '***************************'
    print *, 'Calculating Average Density'
    print *, '***************************'
    print *, ''


    Mtmp = 0.0

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,Lx,Ly,Lz)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:Mtmp) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax

                Lx = (rho(i-1,j,k) - 2.0_rp*rho(i,j,k) + rho(i+1,j,k))
                Ly = (rho(i,j-1,k) - 2.0_rp*rho(i,j,k) + rho(i,j+1,k))
                Lz = (rho(i,j,k-1) - 2.0_rp*rho(i,j,k) + rho(i,j,k+1))
                rho4(i,j,k) = rho(i,j,k) + ((Lx + Ly + Lz)/24.0_rp)
                Mtmp = Mtmp + rho4(i,j,k)

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    Mtmp = Mtmp*h**3

    print *, M, Mtmp

 end subroutine CalculateAvgDensity
