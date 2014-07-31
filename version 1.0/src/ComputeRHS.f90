 !****************************************************************  
 ! ComputeRHS.f90
 ! 
 ! Created on: Jul 26, 2014
 ! Author: ninoy
 !***************************************************************** 

SUBROUTINE ComputeRHS
#include <finclude/petscdef.h>
    use petscksp
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    use MD_Helper
    use MD_CordinateTransform
    use MD_BoundaryCondition
    implicit none

    PetscErrorCode ierr
    PetscScalar bi(1)
    PetscOffset offset
    PetscInt i,j,k
    PetscInt nghi,nghj,nghk
    PetscInt row
    real(kind=double) ,dimension(0:2) :: x
    PetscViewer viewer

    print *, ''
    print *, '****************'
    print *, 'Calculating RHS'
    print *, '****************'
    print *, ''


    call VecZeroEntries(b,ierr)
    call VecGetArray(b,bi,offset,ierr)

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                row = GetIndex(i,j,k)
                bi(offset + 1 + row) = 4*pi*G*rho(i,j,k)*h**3
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=imin
    i=imin
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j = jmin, jmax
        do k = kmin, kmax

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=imin+1
    i=imin+1
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j = jmin, jmax
        do k = kmin, kmax

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    ! cell adjacent to boundary i=imax
    i=imax
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j = jmin, jmax
        do k = kmin, kmax

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=imax-1
    i=imax-1
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j = jmin, jmax
        do k = kmin, kmax

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    ! cell adjacent to boundary j=jmin
    j=jmin
    !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i = imin+2, imax-2
        do k = kmin, kmax

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary j=jmin+1
    j=jmin+1
    !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i = imin+2, imax-2
        do k = kmin, kmax

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    ! cell adjacent to boundary j=jmax
    j=jmax
    !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i = imin+2, imax-2
        do k = kmin, kmax

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary j=jmax-1
    j=jmax-1
    !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i = imin+2, imax-2
        do k = kmin, kmax

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary k=kmin
    k=kmin
    !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i = imin+2, imax-2
        do j = jmin+2, jmax-2

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary k=kmin+1
    k=kmin+1
    !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i = imin+2, imax-2
        do j = jmin+2, jmax-2

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary k=kmax
    k=kmax
    !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i = imin+2, imax-2
        do j = jmin+2, jmax-2

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary k=kmax
    k=kmax-1
    !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i = imin+2, imax-2
        do j = jmin+2, jmax-2

            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(nghi,nghj,nghk) .ne. 0.0) then
                            row = GetIndex(i,j,k)
                            call GetCordinate(i+nghi,j+nghj,k+nghk,x)
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(nghi,nghj,nghk)*GetBoundaryValue(x)
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call VecRestoreArray(b,bi,offset,ierr)
    call VecAssemblyBegin(b,ierr)
    call VecAssemblyEnd(b,ierr)




END SUBROUTINE ComputeRHS
