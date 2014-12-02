 !****************************************************************  
 ! Output_1D.f90
 ! 
 ! Created on: Aug 13, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine Output_1D
#include <finclude/petscdef.h>
    use MD_Parameter
    use MD_Quantity
    use MD_Helper
    use MD_BoundaryCondition
    implicit none

    integer :: i,j,k,row
    integer :: ihalf, jhalf, khalf, ierror
    PetscScalar phii(1)
    PetscScalar bi(1)
    PetscOffset offset,offset1
    PetscErrorCode ierr
    character(len=100) :: format

    ihalf = (imin + imax)/2
    jhalf = (jmin + jmax)/2
    khalf = (kmin + kmax)/2

    open(unit=1000, file='data/rho_x.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1001, file='data/rho_y.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1002, file='data/rho_z.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1003, file='data/phi_x.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1004, file='data/phi_y.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1005, file='data/phi_z.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1006, file='data/error_x.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1007, file='data/error_y.dat',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)

    open(unit=1008, file='data/error_z.dat',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)

    open(unit=1009, file='data/exact_x.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1010, file='data/exact_y.dat',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)

    open(unit=1011, file='data/exact_z.dat',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)

    open(unit=1012, file='data/rhs_x.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1013, file='data/rhs_y.dat',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)

    open(unit=1014, file='data/rhs_z.dat',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)

    call VecGetArray(phi,phii,offset,ierr)
    call VecGetArray(b,bi,offset1,ierr)

    format = "(E20.3E5, A4, E20.3E5)"

    j = jhalf
    k = khalf
    do i = imin, imax
        row = GetIndex(i,j,k)
        write(1000,*) x(i), '    ', rho(i,j,k)
        write(1003,*) x(i), '    ', PetscRealPart(phii(offset + 1 + row))
        write(1006,*) x(i), '    ', error(i,j,k)
        write(1009,*) x(i), '    ', GetBoundaryValue(x(i),y(j),z(k))
        write(1012,*) x(i), '    ', PetscRealPart(bi(offset1 + 1 + row))
!        print *, x(i),y(j),z(k),sqrt(x(i)**2 + y(j)**2 + z(k)**2),GetBoundaryValue(x(i),y(j),z(k))
    end do

    i = ihalf
    k = khalf
    do j = jmin, jmax
        row = GetIndex(i,j,k)
        write(1001,*) y(j), '    ', rho(i,j,k)
        write(1004,*) y(j), '    ', PetscRealPart(phii(offset + 1 + row))
        write(1007,*) y(j), '    ', error(i,j,k)
        write(1010,*) y(j), '    ', GetBoundaryValue(x(i),y(j),z(k))
        write(1013,*) y(j), '    ', PetscRealPart(bi(offset + 1 + row))
    end do

    i = ihalf
    j = jhalf
    do k = kmin, kmax
        row = GetIndex(i,j,k)
        write(1002,*) z(k), '    ', rho(i,j,k)
        write(1005,*) z(k), '    ', PetscRealPart(phii(offset + 1 + row))
        write(1008,*) z(k), '    ', error(i,j,k)
        write(1011,*) z(k), '    ', GetBoundaryValue(x(i),y(j),z(k))
        write(1014,*) z(k), '    ', PetscRealPart(bi(offset + 1 + row))
    end do

    call VecRestoreArray(phi,phii,offset,ierr)
    call VecRestoreArray(b,bi,offset1,ierr)

    do i = 1000, 1014
        close(unit=i)
    end do

end subroutine Output_1D
