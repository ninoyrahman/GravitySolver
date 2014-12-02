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

    ! inner cell number
    integer, parameter :: inc = 80
    integer, parameter :: jnc = 80
    integer, parameter :: knc = 80

    ! ghost cell number
    integer, parameter :: ngh = 4

    ! grid size
    real(KIND=double), parameter :: h = 2.0_rp/dble(inc)

    ! starting and ending index
    integer, parameter :: imin = ngh
    integer, parameter :: imax = inc + imin - 1
    integer, parameter :: jmin = ngh
    integer, parameter :: jmax = jnc + jmin - 1
    integer, parameter :: kmin = ngh
    integer, parameter :: kmax = knc + kmin - 1

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
