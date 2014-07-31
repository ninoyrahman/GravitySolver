 !****************************************************************  
 ! IC_CentredSphere.f90
 ! 
 ! Created on: Jul 26, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine IC_CentredSphere
    use MD_Parameter
    use MD_Quantity
    use MD_CordinateTransform
    implicit none

    integer :: i,j,k,counter
    real(kind=double) :: r,r2
    real(kind=double), dimension(0:2) :: x

    print *, ''
    print *, '**************************************'
    print *, 'Initializing Scenerio: Centred Sphere'
    print *, '**************************************'
    print *, ''

    xc(0) = 0.0
    xc(1) = 0.0
    xc(2) = 0.0
    M = 0.0
    r = 0.1
    counter = 0

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,x,r2)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:counter)
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                call GetCordinate(i,j,k,x)
                r2 = x(0)**2 + x(1)**2 + x(2)**2
                if(r2 .le. r*r) then
                    counter = counter + 1
                end if
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,x,r2)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M)
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                call GetCordinate(i,j,k,x)
                r2 = x(0)**2 + x(1)**2 + x(2)**2
                if(r2 .le. r*r) then
                    rho(i,j,k) = 1.0/counter
                else
                    rho(i,j,k) = 0.0
                end if
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


end subroutine IC_CentredSphere
