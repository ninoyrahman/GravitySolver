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
    use MD_PETScQuantity
    use MD_Helper
    use MD_BoundaryCondition
    implicit none

    PetscErrorCode ierr, ierror
    PetscInt i,j,k
    PetscInt nghi,nghj,nghk
    PetscInt row
    PetscInt istart, iend, jstart, jend
    real(kind=double) :: ue = 0.0_rp
    character(len=10) :: rankstr
    character(len=10) :: format
    PetscScalar, dimension(:), allocatable :: bi
    PetscInt, dimension(:), allocatable :: rowi

    if(.not. allocated(bi)) allocate(bi(0:ilnc*jlnc*klnc-1))
    if(.not. allocated(rowi)) allocate(rowi(0:ilnc*jlnc*klnc-1))

    if(rank .eq. master) then
        print *, ''
        print *, '****************'
        print *, 'Calculating RHS'
        print *, '****************'
        print *, ''
    endif


    call VecZeroEntries(b,ierr)

    if(igmin .eq. ngh) then
        istart = imin + 2
    else if(igmin .eq. ngh+1) then
        istart = imin + 1
    else
        istart = imin
    end if

    if(igmax .eq. ngh+inc-1) then
        iend = imax - 2
    else if(igmax .eq. ngh+inc-2) then
        iend = imax - 1
    else
        iend = imax
    end if

    if(jgmin .eq. ngh) then
        jstart = jmin + 2
    else if(jgmin .eq. ngh+1) then
        jstart = jmin + 1
    else
        jstart = jmin
    end if

    if(jgmax .eq. ngh+jnc-1) then
        jend = jmax - 2
    else if(jgmax .eq. ngh+jnc-2) then
        jend = jmax - 1
    else
        jend = jmax
    end if

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetLocalIndex(i,j,k)
                bi(row) = 4.0_rp*pi*G*(ue*rho(i,j,k) + (1.0_rp - ue)*rho4(i,j,k))*h**3
                rowi(row) = GetIndex(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=ibottom
    if(igmin .eq. ngh) then
        i=imin
        !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do k = kmin, kmax
            do j = jmin, jmax
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary i=ibottom+1
    if(igmin .le. ngh+1 .and. igmax .ge. ngh+1) then
        i=imin+1
        !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do k = kmin, kmax
            do j = jmin, jmax
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary i=itop
    if(igmax .eq. ngh+inc-1) then
        i=imax
        !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do k = kmin, kmax
            do j = jmin, jmax
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary i=itop-1
    if(igmin .le. ngh+inc-2 .and. igmax .ge. ngh+inc-2) then
        i=imax-1
        !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do k = kmin, kmax
            do j = jmin, jmax
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary j=jbottom
    if(jgmin .eq. ngh) then
        j=jmin
        !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do k = kmin, kmax
            do i = istart, iend
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if


    ! cell adjacent to boundary j=jbottom+1
    if(jgmin .le. ngh+1 .and. jgmax .ge. ngh+1) then
        j=jmin+1
        !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do k = kmin, kmax
            do i = istart, iend
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary j=jtop
    if(jgmax .eq. ngh+jnc-1) then
        j=jmax
        !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do k = kmin, kmax
            do i = istart, iend
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary j=jtop-1
    if(jgmin .le. ngh+jnc-2 .and. jgmax .ge. ngh+jnc-2) then
        j=jmax-1
        !$OMP PARALLEL PRIVATE(i,k,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do k = kmin, kmax
            do i = istart, iend
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary k=kbottom
    if(kgmin .eq. ngh) then
        k=kmin
        !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do j = jstart, jend
            do i = istart, iend
                row = GetLocalIndex(i,j,k)
                do nghi = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghk = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary k=kbottom+1
    if(kgmin .le. ngh+1 .and. kgmax .ge. ngh+1) then
        k=kmin+1
        !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do j = jstart, jend
            do i = istart, iend
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary k=ktop
    if(kgmax .eq. ngh+knc-1) then
        k=kmax
        !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do j = jstart, jend
            do i = istart, iend
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    ! cell adjacent to boundary k=kmax
    if(kgmin .le. ngh+knc-2 .and. kgmax .ge. ngh+knc-2) then
        k=kmax-1
        !$OMP PARALLEL PRIVATE(i,j,nghi,nghj,nghk,row)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
        do j = jstart, jend
            do i = istart, iend
                row = GetLocalIndex(i,j,k)
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then
                                bi(row) = bi(row) - &
                                    S(row,nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                            end if
                        end do
                    end do
                end do

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    call VecSetValues(b,ilnc*jlnc*klnc,rowi,bi,INSERT_VALUES,ierr)

    if(allocated(bi)) deallocate(bi)
    if(allocated(rowi)) deallocate(rowi)

    call VecAssemblyBegin(b,ierr)
    call VecAssemblyEnd(b,ierr)
!    call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)

END SUBROUTINE ComputeRHS
