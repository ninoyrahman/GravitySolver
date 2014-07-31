 !****************************************************************  
 ! ComputeRHS.f90
 ! 
 ! Created on: Jul 24, 2014
 ! Author: ninoy
 !***************************************************************** 

SUBROUTINE ComputeRHS(b)
#include <finclude/petscdef.h>
    use petscksp
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    use MD_Helper
    implicit none

    Vec b
    PetscErrorCode ierr
    PetscScalar bi(1)
    PetscOffset offset
    PetscInt i,j,k
    PetscInt row

    call VecZeroEntries(b,ierr)
    call VecGetArray(b,bi,offset,ierr)

    ! cell adjacent to boundary k=minin
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=minin, maxin
        do j = minin, maxin
            row = GetIndex(i,j,minin)
            bi(offset + 1 + row) = bi(offset + 1 + row) - 16.0 + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary k=maxin
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=minin, maxin
        do j = minin, maxin
            row = GetIndex(i,j,maxin)
            bi(offset + 1 + row) = bi(offset + 1 + row) - 16.0 + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary j=minin
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=minin, maxin
        do k = minin, maxin
            row = GetIndex(i,minin,k)
            bi(offset + 1 + row) = bi(offset + 1 + row) - 16.0 + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary j=maxin
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=minin, maxin
        do k = minin, maxin
            row = GetIndex(i,maxin,k)
            bi(offset + 1 + row) = bi(offset + 1 + row) - 16.0 + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=minin
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j=minin, maxin
        do k = minin, maxin
            row = GetIndex(minin,j,k)
            bi(offset + 1 + row) = bi(offset + 1 + row) - 16.0 + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=maxin
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j=minin, maxin
        do k = minin, maxin
            row = GetIndex(maxin,j,k)
            bi(offset + 1 + row) = bi(offset + 1 + row) - 16.0 + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    ! cell adjacent to boundary k=minin+1
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=minin, maxin
        do j = minin, maxin
            row = GetIndex(i,j,minin+1)
            bi(offset + 1 + row) = bi(offset + 1 + row) + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary k=maxin-1
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=minin, maxin
        do j = minin, maxin
            row = GetIndex(i,j,maxin-1)
            bi(offset + 1 + row) = bi(offset + 1 + row) + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary j=minin+1
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=minin, maxin
        do k = minin, maxin
            row = GetIndex(i,minin+1,k)
            bi(offset + 1 + row) = bi(offset + 1 + row) + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary j=maxin-1
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do i=minin, maxin
        do k = minin, maxin
            row = GetIndex(i,maxin-1,k)
            bi(offset + 1 + row) = bi(offset + 1 + row) + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=minin+1
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j=minin, maxin
        do k = minin, maxin
            row = GetIndex(minin+1,j,k)
            bi(offset + 1 + row) = bi(offset + 1 + row) + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! cell adjacent to boundary i=maxin-1
    !$OMP PARALLEL PRIVATE(i,j,k,row,ierr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    do j=minin, maxin
        do k = minin, maxin
            row = GetIndex(maxin-1,j,k)
            bi(offset + 1 + row) = bi(offset + 1 + row) + 1.0
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call VecRestoreArray(b,bi,offset,ierr)

    call VecAssemblyBegin(b,ierr)
    call VecAssemblyEnd(b,ierr)

END SUBROUTINE ComputeRHS
