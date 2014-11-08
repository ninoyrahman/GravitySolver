PROGRAM main
#include <finclude/petscdef.h>
    USE petscksp
    USE MD_Definition
    USE MD_Parameter
    USE MD_Quantity
    USE MD_Helper
    USE MD_IO
    IMPLICIT NONE

    integer :: iteration
    PetscErrorCode ierr

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    call ReadParameterFile

    if(test .eq. 0) then
        call Init
        call Solver(iteration)
        call Output_3D(0)
    else
        call RunTest
    end if

    call Finalize

    call PetscFinalize(ierr)

END PROGRAM main
