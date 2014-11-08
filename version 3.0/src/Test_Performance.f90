!****************************************************************
! Test_Performance.f90
!
! Created on: Nov 3, 2014
! Author: ninoy
!*****************************************************************
subroutine Test_Performance
    use MD_Parameter
    implicit none
#include <f90papi.h>

    integer(kind=8) :: start_time, end_time
    real(kind=double) :: init_time, comp_time, out_time
    integer :: ierr
    integer :: iteration
    character(len=100) :: format

    print *, ''
    print *, '*****************'
    print *, 'Performance Test'
    print *, '*****************'
    print *, ''

    coordinate = 0
    scenerio = 4

    call PAPIF_library_init(ierr)
    if(ierr .eq. PAPI_VER_CURRENT) then
        print *, 'PAPI library initialized.'
    else if(ierr .gt. 0) then
        print *, 'error: PAPI library version mis-match.'
        call exit(1)
    else
        print *, 'error: PAPI initialization error.'
        call exit(1)
    end if

    call PAPIF_get_real_usec(start_time)
    call Init
    call PAPIF_get_real_usec(end_time)
    init_time = dble(end_time - start_time)*1.0d-6

    call PAPIF_get_real_usec(start_time)
    call Solver(iteration)
    call PAPIF_get_real_usec(end_time)
    comp_time = dble(end_time - start_time)*1.0d-6

    call PAPIF_get_real_usec(start_time)
    call Output_3D(0)
    call PAPIF_get_real_usec(end_time)
    out_time = dble(end_time - start_time)*1.0d-6

    call MPI_Barrier(cart_comm, ierr)
    call MPI_Allreduce(init_time, init_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, cart_comm, ierr)
    call MPI_Allreduce(comp_time, comp_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, cart_comm, ierr)
    call MPI_Allreduce(out_time, out_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, cart_comm, ierr)

    if(rank .eq. master) then

        format = "(A20, I5, I5, I5)"
        write (*,format) 'proc distribution=',inp,jnp,knp
        write (*,format) 'resolution=',inc,jnc,knc
        format = "(A12, E10.3)"
        write (*,format) 'init time=',init_time
        write (*,format) 'comp time=',comp_time
        write (*,format) 'out time=',out_time

    end if


end subroutine Test_Performance
