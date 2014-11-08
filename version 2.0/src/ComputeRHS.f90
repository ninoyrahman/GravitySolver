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
    use MD_BoundaryCondition
    implicit none

    PetscErrorCode ierr
    PetscScalar bi(1)
    PetscOffset offset
    PetscInt i,j,k
    PetscInt nghi,nghj,nghk
    PetscInt row
    PetscViewer viewer
    real(kind=double) :: ue = 0.0_rp

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
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetIndex(i,j,k)
                bi(offset + 1 + row) = 4.0_rp*pi*G*(ue*rho(i,j,k) + (1.0_rp - ue)*rho4(i,j,k))*h**3
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=imin
    i=imin
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do j = jmin, jmax
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do j = jmin, jmax
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do j = jmin, jmax
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do j = jmin, jmax
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do i = imin+2, imax-2
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do i = imin+2, imax-2
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do i = imin+2, imax-2
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do i = imin+2, imax-2
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j = jmin+2, jmax-2
        do i = imin+2, imax-2
            row = GetIndex(i,j,k)
            do nghi = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghk = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j = jmin+2, jmax-2
        do i = imin+2, imax-2
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j = jmin+2, jmax-2
        do i = imin+2, imax-2
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j = jmin+2, jmax-2
        do i = imin+2, imax-2
            row = GetIndex(i,j,k)
            do nghk = -swidth, +swidth
                do nghj = -swidth, +swidth
                    do nghi = -swidth, +swidth
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                            bi(offset + 1 + row) = bi(offset + 1 + row) - &
                                S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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

    call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

END SUBROUTINE ComputeRHS
