!****************************************************************
! Test_Indexing.f90
!
! Created on: Jul 25, 2014
! Author: ninoy
!*****************************************************************

SUBROUTINE Test_Indexing
    use MD_Parameter
    use MD_Helper
    implicit none

    integer :: i,j,k
    integer :: row,column,n

    print *, ''
    print *, '**************'
    print *, 'Indexing Test'
    print *, '**************'
    print *, ''

    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                print *, 'i=',i,'j=',j,'k=',k,'index=',GetIndex(i,j,k)
            end do
        end do
    end do

    n = inc*jnc*knc-1

    do row = 0, n
        do column = 0, n
            print *, 'row=',row,'column=',column,'index=',GetLinearIndex(row,column)
        end do
    end do

END SUBROUTINE Test_Indexing
