!****************************************************************
! ComputeMatrix.f90
!
! Created on: Jul 24, 2014
! Author: ninoy
!*****************************************************************

SUBROUTINE ComputeMatrix(A)
#include <finclude/petscdef.h>
    use petscksp
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    use MD_Helper
    implicit none

    Mat A
    PetscErrorCode ierr
    PetscInt i,j,k
    PetscInt nghi,nghj,nghk
    PetscInt row, column
    PetscScalar Aij

    ! main diagonal column = row
    do i=minin,maxin
        do j = minin, maxin
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row
                Aij = -90.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row - 1
    do i=minin,maxin
        do j = minin, maxin
            do k = minin+1, maxin
                row = GetIndex(i,j,k)
                column = row - 1
                Aij = 16.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row + 1
    do i=minin,maxin
        do j = minin, maxin
            do k = minin, maxin-1
                row = GetIndex(i,j,k)
                column = row + 1
                Aij = 16.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row - 2
    do i=minin,maxin
        do j = minin, maxin
            do k = minin+2, maxin
                row = GetIndex(i,j,k)
                column = row - 2
                Aij = -1.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row + 2
    do i=minin,maxin
        do j = minin, maxin
            do k = minin, maxin-2
                row = GetIndex(i,j,k)
                column = row + 2
                Aij = -1.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row - nin
    do i=minin,maxin
        do j = minin+1, maxin
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row - nin
                Aij = 16.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row + nin
    do i=minin, maxin
        do j = minin, maxin-1
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row + nin
                Aij = 16.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

        ! off diagonal column = row - 2*nin
    do i=minin,maxin
        do j = minin+2, maxin
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row - 2*nin
                Aij = -1.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row + 2*nin
    do i=minin, maxin
        do j = minin, maxin-2
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row + 2*nin
                Aij = -1.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row - nin*nin
    do i=minin+1,maxin
        do j = minin, maxin
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row - nin*nin
                Aij = 16.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row + nin*nin
    do i=minin, maxin-1
        do j = minin, maxin
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row + nin*nin
                Aij = 16.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row - 2*nin*nin
    do i=minin+2,maxin
        do j = minin, maxin
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row - 2*nin*nin
                Aij = -1.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    ! off diagonal column = row + 2*nin*nin
    do i=minin, maxin-2
        do j = minin, maxin
            do k = minin, maxin
                row = GetIndex(i,j,k)
                column = row + 2*nin*nin
                Aij = -1.0
                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)
            end do
        end do
    end do

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

END SUBROUTINE ComputeMatrix
