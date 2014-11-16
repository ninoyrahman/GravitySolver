!****************************************************************
! Test_Indexing.f90
!
! Created on: Jul 25, 2014
! Author: ninoy
!*****************************************************************

SUBROUTINE Test_Indexing
    use MD_Parameter
    use MD_Helper
    use MD_CalculateCordinate
    implicit none

    integer :: i,j,k
    integer :: row,column,n
    integer :: whichIndexing, whichRank

    if(rank .eq. master) then
        print *, ''
        print *, '**************'
        print *, 'Indexing Test'
        print *, '**************'
        print *, ''
    endif

    h = 2.0_rp/dble(inc)
    whichIndexing = 1   !0=global normal, 1=global PETSc, 2=local, 3=linear
    whichRank = 0

    call DomainDecomposition
    call CalculateCartesianCordinate

    print *, 'rank: ', rank, 'igmin=', igmin, 'igmax=', igmax
    print *, 'rank: ', rank, 'jgmin=', jgmin, 'jgmax=', jgmax
    print *, 'rank: ', rank, 'kgmin=', kgmin, 'kgmax=', kgmax

    if(rank .eq. whichRank .and. whichIndexing .eq. 0) then
        do k = kmin, kmax
            do j = jmin, jmax
                do i = imin, imax
                    print *, 'rank=', rank, 'i=',i,'j=',j,'k=',k,'index=',GetIndex(i,j,k)
                end do
            end do
        end do
    end if

    if(rank .eq. whichRank .and. whichIndexing .eq. 1) then
        do k = kmin, kmax
            do j = jmin, jmax
                do i = imin, imax
                    print *, 'rank=', rank, 'i=',i,'j=',j,'k=',k,'index=',GetPETScIndex(i,j,k)
                end do
            end do
        end do
    end if

    if(rank .eq. whichRank .and. whichIndexing .eq. 2) then
        do k = kmin, kmax
            do j = jmin, jmax
                do i = imin, imax
                    print *, 'rank=', rank, 'i=',i,'j=',j,'k=',k,'index=',GetLocalIndex(i,j,k)
                end do
            end do
        end do
    end if

    n = inc*jnc*knc-1

    if(rank .eq. whichRank .and. whichIndexing .eq. 3) then
        do row = 0, n
            do column = 0, n
                print *, 'row=',row,'column=',column,'index=',GetLinearIndex(row,column)
            end do
        end do
    end if


END SUBROUTINE Test_Indexing
