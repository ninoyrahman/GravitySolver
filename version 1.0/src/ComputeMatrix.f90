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
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                row = GetIndex(i,j,k)

                ! looping over inner cell in stencil
                do nghi = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghk = -swidth, +swidth
                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 0 .and. S(nghi,nghj,nghk) .ne. 0.0) then

                                Aij = S(nghi,nghj,nghk)
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

END SUBROUTINE ComputeMatrix
