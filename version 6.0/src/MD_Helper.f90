 !****************************************************************  
 ! MD_Helper.f90
 ! 
 ! Created on: Jul 24, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_Helper
#include <finclude/petscdef.h>
    USE petscksp
    USE MD_Parameter
    USE MD_Quantity
    IMPLICIT NONE

CONTAINS

    !***************************************************************
    ! returns 1 if cell(i,j,k) is a boundary cell
    ! returns 0 if cell(i,j,k) is a inner cell
    !***************************************************************
    integer function IsBoundary(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        integer :: ig,jg,kg
        ig = i-ngh + igmin
        jg = j-ngh + jgmin
        kg = k-ngh + kgmin
        if(ig .lt. ngh .or. ig .gt. ngh+inc-1 .or. &
            jg .lt. ngh .or. jg .gt. ngh+jnc-1 .or. &
            kg .lt. ngh .or. kg .gt. ngh+knc-1) then
            IsBoundary = 1
        else
            IsBoundary = 0
        end if
    end function IsBoundary

    !***************************************************************
    ! this funtion convert multidimensional local index(i,j,k)
    ! of inner cell to linear normal global index starting from zero.
    ! cell index: (2,2,3)-->1
    ! with 2 boundary cell in each dimension
    !***************************************************************
    integer function GetIndex(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        GetIndex = (i-ngh+igmin-ngh) + inc*(j-ngh+jgmin-ngh) + inc*jnc*(k-ngh+kgmin-ngh)
    end function GetIndex

    !***************************************************************
    ! this funtion convert multidimensional local index(i,j,k)
    ! of inner cell to linear PETSc global index starting from zero.
    ! cell index: (2,2,3)-->1
    ! with 2 boundary cell in each dimension
    !***************************************************************
    integer function GetPETScIndex(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        GetPETScIndex = (igmin-ngh)*jlnc*klnc + (jgmin-ngh)*inc*klnc + (kgmin-ngh)*inc*jnc + &
            (i-imin) + ilnc*(j-jmin) + ilnc*jlnc*(k-kmin)
    end function GetPETScIndex

    !***************************************************************
    ! this funtion convert multidimensional local index(i,j,k)
    ! of inner cell to linear local index starting from zero.
    !***************************************************************
    integer function GetLocalIndex(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        GetLocalIndex = (i-imin) + ilnc*(j-jmin) + ilnc*jlnc*(k-kmin)
    end function GetLocalIndex

    !***************************************************************
    ! this funtion convert multidimensional index(row,column)
    ! of matrix entry to linear index starting from zero.
    ! cell index: (1,1)-->3
    ! with 2 column in each row
    !***************************************************************
    integer function GetLinearIndex(row,column)
        implicit none
        integer, intent(in) :: row,column
        GetLinearIndex = column + row*inc*jnc*knc
    end function GetLinearIndex

    !***************************************************************
    ! this function computes local and global index of cell
    !***************************************************************
    subroutine ComputeIndex
        implicit none
        integer :: i,j,k
        if(rank .eq. master) then
            print *, ''
            print *, '***************'
            print *, 'Computing Index'
            print *, '***************'
            print *, ''
        endif

        !$OMP PARALLEL PRIVATE(i,j,k)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
        do k = kmin, kmax
            do j = jmin, jmax
                do i = imin, imax
                    localindex(i,j,k) = GetLocalIndex(i,j,k)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL PRIVATE(i,j,k)
        !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
        do k = kmin-swidth, kmax+swidth
            do j = jmin-swidth, jmax+swidth
                do i = imin-swidth, imax+swidth
                    globalindex(i,j,k) = GetIndex(i,j,k)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine ComputeIndex

    !***************************************************************
    ! returns random double precision real number
    !***************************************************************
    real(kind=double) function ran2(idum)
        integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        real(kind=double) :: AM,EPS,RNMX
        parameter(IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
            IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
            NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        integer :: idum2,j,k,iv(NTAB),iy
        save iv,iy,idum2
        data idum2/123456789/, iv/NTAB*0/, iy/0/
        if (idum.le.0) then
            idum=max(-idum,1)
            idum2=idum
            do 11 j=NTAB+8,1,-1

                k=idum/IQ1
                idum=IA1*(idum-k*IQ1)-k*IR1
                if (idum.lt.0) idum=idum+IM1
                if (j.le.NTAB) iv(j)=idum
11          continue
            iy=iv(1)
        endif
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idum
        if(iy.lt.1)iy=iy+IMM1
        ran2=min(AM*iy,RNMX)
    end function ran2

END MODULE MD_Helper

