 !****************************************************************  
 ! ComputeMatrix.f90
 ! 
 ! Created on: Jul 26, 2014
 ! Author: ninoy
 !***************************************************************** 

SUBROUTINE ComputeMatrix
#include <finclude/petscdef.h>
    use petscksp
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    use MD_Helper
    implicit none

    PetscErrorCode ierr
    PetscInt i,j,k
    PetscInt nghi,nghj,nghk
    PetscInt row, column
    PetscScalar Aij

    print *, ''
    print *, '***************************'
    print *, 'Calculating System Matrix'
    print *, '***************************'
    print *, ''


    ! looping over inner cell
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetIndex(i,j,k)

                ! looping over inner cell in stencil
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth

                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 0 .and. S(row,nghi,nghj,nghk) .ne. 0.0_rp) then

                                Aij = S(row,nghi,nghj,nghk)!*576.0_rp/h
                                column = GetIndex(i+nghi,j+nghj,k+nghk)
                                call MatSetValue(A,row,column,Aij,INSERT_VALUES,ierr)

                            end if
                        end do
                    end do
                end do

            end do
        end do
    end do

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

!    call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

END SUBROUTINE ComputeMatrix
