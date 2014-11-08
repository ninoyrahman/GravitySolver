 !****************************************************************  
 ! MD_GridParameter.f90
 ! 
 ! Created on: Jul 23, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_Parameter
    USE MD_Definition
    IMPLICIT NONE
    SAVE

    !********************************************
    !       computational grid parameter        !
    !********************************************

    ! mpi parameter
    integer :: nproc, rank, master
    parameter (master = 0)
    integer :: inp, jnp, knp
    integer, dimension(0:2) :: coord

    ! global inner cell number
    integer :: inc, jnc, knc
    ! local inner cell number
    integer :: ilnc, jlnc, klnc

    ! ghost cell number
    integer, parameter :: ngh = 4

    ! grid size
    real(KIND=double) :: h

    ! starting and ending index
    integer :: imin, imax, jmin
    integer :: jmax, kmin, kmax
    integer :: igmin, jgmin, kgmin
    integer :: igmax, jgmax, kgmax

    ! stencil width
    integer, parameter :: swidth = 2

    !********************************************
    !           setting parameter               !
    !********************************************

    integer :: test = 0
    integer :: coordinate = 0
    integer :: scenerio = 1
    real(kind=double) :: M
    real(kind=double) :: R
    real(kind=double) :: xc(0:2)


END MODULE MD_Parameter
