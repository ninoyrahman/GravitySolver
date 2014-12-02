 !****************************************************************  
 ! Output_3D.f90
 ! 
 ! Created on: Jul 23, 2014
 ! Author: ninoy
 !***************************************************************** 

SUBROUTINE Output_3D(time)
#include <finclude/petscdef.h>
    use petscksp
    use MD_Parameter
    use MD_Quantity
    use MD_Helper
    implicit none

    PetscScalar phii(1)
    PetscOffset offset
    PetscErrorCode ierr
    integer, intent(in) :: time
    integer :: i,j,k,row,ierror
!    real(kind=double), dimension(0:2) :: xcord
    character(len=200) :: format
    format = "(A10, I3)"

    print *, ''
    print *, '************'
    print *, 'writing vtk'
    print *, '************'
    print *, ''
    print format, 'time=',time

    call VecGetArray(phi,phii,offset,ierr)

    open(unit=0, file='data/GravitationalPotential.vtk',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)

100 FORMAT(F8.4,F8.4,F8.4)

    write(0,'(a)')"# vtk DataFile Version 2.0"
    write(0,'(a)')' Generated by Ninoy Rahman'
    write(0,'(a)')'ASCII'
    write(0,*)''
    write(0,'(a)')'DATASET UNSTRUCTURED_GRID'
    write(0,*)'POINTS ', inc*jnc*knc ,' double'

    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
!                call GetCordinate(i,j,k,xcord)
                write(0,100) x(i),y(j),z(k)
            end do
        end do
    end do

    write(0,*) ''
    write(0,*) 'POINT_DATA ', inc*jnc*knc
    write(0,'(a)') 'SCALARS Potential double 1'
    write(0,'(a)') 'LOOKUP_TABLE DEFAULT'

    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetIndex(i,j,k)
                write(0,*) PetscRealPart(phii(offset + 1 + row))
            end do
        end do
    end do

    write(0,*) ''
    write(0,'(a)') 'SCALARS Density double 1'
    write(0,'(a)') 'LOOKUP_TABLE DEFAULT'

    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                write(0,*) rho(i,j,k)
            end do
        end do
    end do

    call VecRestoreArray(phi,phii,offset,ierr)

    close(unit=0)

END SUBROUTINE Output_3D
