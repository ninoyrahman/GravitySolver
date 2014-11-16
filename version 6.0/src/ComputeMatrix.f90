 !****************************************************************  
 ! ComputeMatrix.f90
 ! 
 ! Created on: Jul 26, 2014
 ! Author: ninoy
 !***************************************************************** 

SUBROUTINE ComputeMatrix
#include <finclude/petscdef.h>
    use petscksp
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    use MD_PETScQuantity
    use MD_Helper
    implicit none

    Vec dnnz
    Vec onnz
    PetscErrorCode ierr
    PetscInt i,j,k
    PetscInt proc,whichproc,buffer
    PetscInt nghi,nghj,nghk
    PetscInt row,lrow,column
    PetscScalar, dimension(:,:), allocatable :: aa
    PetscInt, dimension(:,:), allocatable :: icolumn
    PetscInt, dimension(:), allocatable :: ncolumn
    PetscInt, dimension(:), allocatable :: low, high
    PetscScalar, dimension(:), allocatable :: dnnzi, onnzi
    PetscInt, dimension(:), allocatable :: dnnzii, onnzii
    PetscInt, dimension(:), allocatable :: irow
    PetscScalar, pointer :: di(:)
    PetscScalar, pointer :: oi(:)

    if(rank .eq. master) then
        print *, ''
        print *, '***************************'
        print *, 'Calculating System Matrix'
        print *, '***************************'
        print *, ''
    end if

    if(.not. allocated(aa)) allocate(aa(0:ilnc*jlnc*klnc-1,0:61))
    if(.not. allocated(icolumn)) allocate(icolumn(0:ilnc*jlnc*klnc-1,0:61))
    if(.not. allocated(ncolumn)) allocate(ncolumn(0:ilnc*jlnc*klnc-1))
    if(.not. allocated(dnnzi)) allocate(dnnzi(0:ilnc*jlnc*klnc-1))
    if(.not. allocated(onnzi)) allocate(onnzi(0:ilnc*jlnc*klnc-1))
    if(.not. allocated(dnnzii)) allocate(dnnzii(0:ilnc*jlnc*klnc-1))
    if(.not. allocated(onnzii)) allocate(onnzii(0:ilnc*jlnc*klnc-1))
    if(.not. allocated(irow)) allocate(irow(0:ilnc*jlnc*klnc-1))
    if(.not. allocated(low)) allocate(low(0:nproc))
    if(.not. allocated(high)) allocate(high(0:nproc))

    call VecDuplicate(b, dnnz, ierr)
    call VecDuplicate(b, onnz, ierr)
    call VecGetOwnershipRange(dnnz, low(rank), high(rank), ierr)
    high(rank) = high(rank) - 1

    call MPI_Barrier(cart_comm, ierr)
    buffer = low(rank)
    call MPI_Allgather(buffer, 1, MPI_INT, low, 1, MPI_INT, cart_comm, ierr)
    buffer = high(rank)
    call MPI_Allgather(buffer, 1, MPI_INT, high, 1, MPI_INT, cart_comm, ierr)


    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,nghi,nghj,nghk,row,lrow,column,proc,whichproc)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax

                row = globalindex(i,j,k)
                lrow = localindex(i,j,k)
                column = 0

                do proc = 0, nproc-1
                    if(row .ge. low(proc) .and. row .le. high(proc)) whichproc = proc
                end do

                irow(lrow) = row
                dnnzi(lrow) = 0
                onnzi(lrow) = 0

                ! looping over inner cell in stencil
                do nghk = -swidth, +swidth
                    do nghj = -swidth, +swidth
                        do nghi = -swidth, +swidth

                            if(IsBoundary(i+nghi,j+nghj,k+nghk) .eq. 0 .and. S(lrow,nghi,nghj,nghk) .ne. 0.0_rp) then

                                icolumn(lrow,column) = globalindex(i+nghi,j+nghj,k+nghk)
                                aa(lrow,column) = S(lrow,nghi,nghj,nghk)!*576.0_rp/h


                                if(icolumn(lrow,column) .ge. low(whichproc) .and. &
                                    icolumn(lrow,column) .le. high(whichproc)) then
                                    dnnzi(lrow) = dnnzi(lrow) + 1
                                else
                                    onnzi(lrow) = onnzi(lrow) + 1
                                end if

                                column = column + 1
                            end if
                        end do
                    end do
                end do

                ncolumn(lrow) = column

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call VecSetValues(dnnz,ilnc*jlnc*klnc,irow,dnnzi,INSERT_VALUES,ierr)
    call VecAssemblyBegin(dnnz,ierr)
    call VecAssemblyEnd(dnnz,ierr)
    call VecSetValues(onnz,ilnc*jlnc*klnc,irow,onnzi,INSERT_VALUES,ierr)
    call VecAssemblyBegin(onnz,ierr)
    call VecAssemblyEnd(onnz,ierr)

    call VecGetArrayF90(dnnz, di, ierr)
    call VecGetArrayF90(onnz, oi, ierr)

    !$OMP PARALLEL PRIVATE(lrow)
    !$OMP DO SCHEDULE(STATIC)
    do lrow = 0, ilnc*jlnc*klnc-1
        dnnzii(lrow) = int(PetscRealPart(di(lrow+1)))
        onnzii(lrow) = int(PetscRealPart(oi(lrow+1)))
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call MatCreate(cart_comm,A,ierr)
    call MatSetSizes(A,ilnc*jlnc*klnc,ilnc*jlnc*klnc,inc*jnc*knc,inc*jnc*knc,ierr)
    call MatSetType(A,MATMPIAIJ,ierr)
    call MatMPIAIJSetPreallocation(A,0,dnnzii,0,onnzii,ierr)
    call MatSetUp(A,ierr)

    call VecRestoreArrayF90(dnnz, di, ierr)
    call VecRestoreArrayF90(onnz, oi, ierr)

    ! looping over inner cell
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax

                row = globalindex(i,j,k)
                lrow = localindex(i,j,k)

                call MatSetValues(A,1,row,ncolumn(lrow),icolumn(lrow,:),aa(lrow,:),INSERT_VALUES,ierr)

            end do
        end do
    end do

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

    if(allocated(aa)) deallocate(aa)
    if(allocated(icolumn)) deallocate(icolumn)
    if(allocated(ncolumn)) deallocate(ncolumn)
    if(allocated(dnnzi)) deallocate(dnnzi)
    if(allocated(onnzi)) deallocate(onnzi)
    if(allocated(dnnzii)) deallocate(dnnzii)
    if(allocated(onnzii)) deallocate(onnzii)
    if(allocated(irow)) deallocate(irow)
    if(allocated(low)) deallocate(low)
    if(allocated(high)) deallocate(high)

!    call MatView(A, PETSC_VIEWER_MATLAB_WORLD, ierr)

END SUBROUTINE ComputeMatrix
