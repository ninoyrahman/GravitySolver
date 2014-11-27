!****************************************************************
! ParallelGravitySolver.f90
!
! Created on: Sep 29, 2014
! Author: ninoy
!*****************************************************************

program ParallelGravitySolver
#include <finclude/petscdef.h>
    USE petscksp
    USE MD_Parameter
    USE MD_IO
    implicit none

    integer :: time = 0
    integer :: iteration
    integer :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, nproc, ierr)

    if(rank .eq. master) then
        print *, 'Parallel Gravity Solver'
    end if

    call ReadParameterFile

    if(test .eq. 0) then
        call Init
        call Solver(iteration)
        call Output_3D(time)
        time = time + 1
    else
        call RunTest
    end if

    call Finalize
    call PetscFinalize(ierr)

end program ParallelGravitySolver
