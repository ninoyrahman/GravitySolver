 !****************************************************************  
 ! IC_OffCentredSphere.f90
 ! 
 ! Created on: Jul 27, 2014
 ! Author: ninoy
 !***************************************************************** 
subroutine IC_OffCentredSphere
    use MD_Parameter
    use MD_Quantity
    use MD_CordinateTransform
    implicit none

    integer :: i,j,k,counter
    real(kind=double) :: r, r2
    real(kind=double), dimension(0:2) :: x
    real(kind=double), dimension(0:2) :: xtmp

    print *, ''
    print *, '******************************************'
    print *, 'Initializing Scenerio: Off-Centred Sphere'
    print *, '******************************************'
    print *, ''

    xc(0) = 0.1
    xc(1) = 0.0
    xc(2) = 0.0
    M = 0.0
    r = 0.1
    counter = 0
    xtmp(0) = 0.0
    xtmp(1) = 0.0
    xtmp(2) = 0.0

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,x,r2)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3) REDUCTION(+:counter)
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                call GetCordinate(i,j,k,x)
                r2 = (x(0)-xc(0))*(x(0)-xc(0))+(x(1)-xc(1))*(x(1)-xc(1))+(x(2)-xc(2))*(x(2)-xc(2))
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
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3) REDUCTION(+:M)
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                call GetCordinate(i,j,k,x)
                r2 = (x(0)-xc(0))*(x(0)-xc(0))+(x(1)-xc(1))*(x(1)-xc(1))+(x(2)-xc(2))*(x(2)-xc(2))
                if(r2 .le. r*r) then
                    rho(i,j,k) = 1.0/counter
                    xtmp(0) = xtmp(0) + x(0)*rho(i,j,k)
                    xtmp(1) = xtmp(1) + x(1)*rho(i,j,k)
                    xtmp(2) = xtmp(2) + x(2)*rho(i,j,k)
                else
                    rho(i,j,k) = 0.0
                end if
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    xc(0) = xtmp(0)/M
    xc(1) = xtmp(1)/M
    xc(2) = xtmp(2)/M


end subroutine IC_OffCentredSphere
