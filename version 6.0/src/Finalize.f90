!****************************************************************
! Finalize.f90
!
! Created on: Jul 26, 2014
! Author: ninoy
!*****************************************************************

subroutine Finalize
#include <finclude/petscdef.h>
    use petscksp
    use MD_Quantity
    use MD_PETScQuantity
    use MD_GeometricQuantity
    implicit none

    PetscErrorCode ierr

    if(rank .eq. master) then
        print *, ''
        print *, '***********'
        print *, 'Finalizing'
        print *, '***********'
        print *, ''
    endif

    call MatDestroy(A,ierr)
    call VecDestroy(b,ierr)
    call VecDestroy(phi,ierr)
    call KSPDestroy(ksp,ierr)

    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)
    if(allocated(z)) deallocate(z)

    if(allocated(rho)) deallocate(rho)
    if(allocated(rho4)) deallocate(rho4)
    if(allocated(globalindex)) deallocate(globalindex)
    if(allocated(localindex)) deallocate(localindex)

    if(allocated(S)) deallocate(S)

    if(allocated(bxx)) deallocate(bxx)
    if(allocated(bxy)) deallocate(bxy)
    if(allocated(bxz)) deallocate(bxz)

    if(allocated(byx)) deallocate(byx)
    if(allocated(byy)) deallocate(byy)
    if(allocated(byz)) deallocate(byz)

    if(allocated(bzx)) deallocate(bzx)
    if(allocated(bzy)) deallocate(bzy)
    if(allocated(bzz)) deallocate(bzz)

    if(allocated(Gxy)) deallocate(Gxy)
    if(allocated(Gxz)) deallocate(Gxz)

    if(allocated(Gyx)) deallocate(Gyx)
    if(allocated(Gyz)) deallocate(Gyz)

    if(allocated(Gzx)) deallocate(Gzx)
    if(allocated(Gzy)) deallocate(Gzy)

end subroutine Finalize
