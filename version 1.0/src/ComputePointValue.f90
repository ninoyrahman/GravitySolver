 !****************************************************************  
 ! ComputePointValue.f90
 ! 
 ! Created on: Nov 6, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine ComputePointValue(cv, pv)
#include <finclude/petscdef.h>
    use petscksp
    use MD_Definition
    use MD_Quantity
    use MD_BoundaryCondition
    use MD_Helper
    implicit none

    Vec, intent(in) :: cv
    Vec, intent(in) :: pv
    Vec rhs
    Mat sysmat
    KSP kry
    PC prec
    PetscErrorCode ierr
    KSPConvergedReason reason
    PetscInt iteration
    PetscScalar, dimension(0:0,0:7) :: sysmata
    PetscInt, dimension(0:7) :: icolumn
    PetscScalar rhsi(1)
    PetscOffset offset
    integer :: i,j,k
    integer :: nghi,nghj,nghk
    integer :: row, column
    integer :: sw
    parameter (sw = 1)
    real(kind=double) :: diag = 0.75_rp
    real(kind=double) :: offdiag = 1.0_rp/24.0_rp
    real(kind=double), dimension(-sw:sw,-sw:sw,-sw:sw) :: St
    character(len=100) :: format

    print *, ''
    print *, '*********************'
    print *, 'Computing Point Value'
    print *, '*********************'
    print *, ''

    do nghk = -sw, +sw
        do nghj = -sw, +sw
            do nghi = -sw, +sw
                St(nghi,nghj,nghk) = 0.0_rp
            end do
        end do
    end do

    St(0,0,0) = diag
    St(-1,0,0) = offdiag
    St(+1,0,0) = offdiag
    St(0,-1,0) = offdiag
    St(0,+1,0) = offdiag
    St(0,0,-1) = offdiag
    St(0,0,+1) = offdiag

    call MatCreate(PETSC_COMM_WORLD,sysmat,ierr)
    call MatSetSizes(sysmat,PETSC_DECIDE,PETSC_DECIDE,inc*jnc*knc,inc*jnc*knc,ierr)
    call MatSetType(sysmat,MATSEQAIJ,ierr)
    call MatSeqAIJSetPreallocation(sysmat,7,PETSC_NULL_INTEGER,ierr)
    call MatSetUp(sysmat,ierr)

    ! looping over inner cell
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax

                row = GetIndex(i,j,k)
                column = 0

                ! looping over inner cell in stencil
                do nghk = -sw, +sw
                    do nghj = -sw, +sw
                        do nghi = -sw, +sw

                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 0 .and. St(nghi,nghj,nghk) .ne. 0.0_rp) then

                                icolumn(column) = GetIndex(i+nghi,j+nghj,k+nghk)
                                sysmata(0,column) = St(nghi,nghj,nghk)
                                column = column + 1

                            end if
                        end do
                    end do
                end do

                call MatSetValues(sysmat,1,row,column,icolumn,sysmata,INSERT_VALUES,ierr)

            end do
        end do
    end do

    call MatAssemblyBegin(sysmat, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(sysmat, MAT_FINAL_ASSEMBLY, ierr)

!    call MatView(sysmat,PETSC_VIEWER_STDOUT_WORLD,ierr)

    call KSPCreate(PETSC_COMM_WORLD,kry,ierr)
    call KSPSetOperators(kry,sysmat,sysmat,SAME_PRECONDITIONER,ierr)
    call KSPSetType(kry,KSPCG,ierr)
    call KSPSetInitialGuessNonzero(kry,PETSC_TRUE,ierr)
    call KSPGetPC(kry,prec,ierr)
    call PCSetType(prec,PCSOR,ierr)
    call KSPSetTolerances(kry,1.0d-17,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION, &
        PETSC_DEFAULT_INTEGER,ierr)

    call VecDuplicate(cv,rhs,ierr)
    call VecCopy(cv,rhs,ierr)
    call VecGetArray(rhs,rhsi,offset,ierr)

    ! cell adjacent to boundary i=imin
    i=imin
    !$OMP PARALLEL PRIVATE(j,k,nghi,nghj,nghk,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do k = kmin, kmax
        do j = jmin, jmax
            row = GetIndex(i,j,k)
            do nghk = -sw, +sw
                do nghj = -sw, +sw
                    do nghi = -sw, +sw
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. St(nghi,nghj,nghk) .ne. 0.0_rp) then
                            rhsi(offset + 1 + row) = rhsi(offset + 1 + row) - &
                                St(nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
            do nghk = -sw, +sw
                do nghj = -sw, +sw
                    do nghi = -sw, +sw
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. St(nghi,nghj,nghk) .ne. 0.0_rp) then
                            rhsi(offset + 1 + row) = rhsi(offset + 1 + row) - &
                                St(nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
        do i = imin+1, imax-1
            row = GetIndex(i,j,k)
            do nghk = -sw, +sw
                do nghj = -sw, +sw
                    do nghi = -sw, +sw
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. St(nghi,nghj,nghk) .ne. 0.0_rp) then
                            rhsi(offset + 1 + row) = rhsi(offset + 1 + row) - &
                                St(nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
        do i = imin+1, imax-1
            row = GetIndex(i,j,k)
            do nghk = -sw, +sw
                do nghj = -sw, +sw
                    do nghi = -sw, +sw
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. St(nghi,nghj,nghk) .ne. 0.0_rp) then
                            rhsi(offset + 1 + row) = rhsi(offset + 1 + row) - &
                                St(nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    do j = jmin+1, jmax-1
        do i = imin+1, imax-1
            row = GetIndex(i,j,k)
            do nghi = -sw, +sw
                do nghj = -sw, +sw
                    do nghk = -sw, +sw
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. St(nghi,nghj,nghk) .ne. 0.0_rp) then
                            rhsi(offset + 1 + row) = rhsi(offset + 1 + row) - &
                                St(nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
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
    do j = jmin+1, jmax-1
        do i = imin+1, imax-1
            row = GetIndex(i,j,k)
            do nghk = -sw, +sw
                do nghj = -sw, +sw
                    do nghi = -sw, +sw
                        if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 1 .and. St(nghi,nghj,nghk) .ne. 0.0_rp) then
                            rhsi(offset + 1 + row) = rhsi(offset + 1 + row) - &
                                St(nghi,nghj,nghk)*GetBoundaryValue(x(i+nghi),y(j+nghj),z(k+nghk))
                        end if
                    end do
                end do
            end do

        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call VecRestoreArray(rhs,rhsi,offset,ierr)
    call VecAssemblyBegin(rhs,ierr)
    call VecAssemblyEnd(rhs,ierr)

    call KSPSolve(kry,rhs,pv,ierr)

    format = "(A20, I3)"
    call KSPView(kry, PETSC_VIEWER_STDOUT_WORLD,ierr)

    call KSPGetConvergedReason(kry,reason,ierr)
    if(reason .lt. 0) then
        print *,'Diverged!!!, ', reason
    else
        call KSPGetIterationNumber(kry,iteration,ierr)
        print format,'iteration number=', iteration
    end if

    call MatDestroy(sysmat,ierr)
    call VecDestroy(rhs,ierr)
    call KSPDestroy(kry,ierr)

end subroutine ComputePointValue
