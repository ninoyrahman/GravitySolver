!****************************************************************
! RunTest.f90
!
! Created on: Jul 24, 2014
! Author: ninoy
!*****************************************************************

subroutine RunTest
    use MD_Parameter
    implicit none

    if(test .eq. 1) then
        call Test_Indexing
    else if(test .eq. 2) then
        call Test_BoundaryValue
    else if(test .eq. 3) then
        call Test_Accuracy
    else if(test .eq. 4) then
        call Test_Performance
    else if(test .eq. 5) then
        call Test_SolverConvergence
    else
        print *, 'Wrong Test Number!'
    end if

end subroutine RunTest
