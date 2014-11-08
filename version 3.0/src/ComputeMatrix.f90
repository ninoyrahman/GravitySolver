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
    use MD_PETScQuantity
    use MD_Helper
    implicit none

    PetscErrorCode ierr
    PetscInt i,j,k
    PetscInt nghi,nghj,nghk
    PetscInt row,lrow,column
    PetscScalar, dimension(0:0,0:61) :: aa
    PetscInt, dimension(0:61) :: icolumn

    if(rank .eq. master) then
        print *, ''
        print *, '***************************'
        print *, 'Calculating System Matrix'
        print *, '***************************'
        print *, ''
    end if

    ! looping over inner cell
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax

                row = GetIndex(i,j,k)
                lrow = GetLocalIndex(i,j,k)
                column = 0

                ! looping over inner cell in stencil
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth

                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 0 .and. S(lrow,nghi,nghj,nghk) .ne. 0.0_rp) then

                                icolumn(column) = GetIndex(i+nghi,j+nghj,k+nghk)
                                aa(0,column) = S(lrow,nghi,nghj,nghk)
                                column = column + 1

                            end if
                        end do
                    end do
                end do

                call MatSetValues(A,1,row,column,icolumn,aa,INSERT_VALUES,ierr)

            end do
        end do
    end do

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

END SUBROUTINE ComputeMatrix
