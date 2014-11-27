!****************************************************************
! DomainDecomposition.f90
!
! Created on: Oct 4, 2014
! Author: ninoy
!*****************************************************************

subroutine DomainDecomposition
    use petscksp
    use MD_Parameter
    use MD_Quantity
    use MD_PETScQuantity
    use MD_GeometricQuantity
    use MD_Helper
    implicit none

    integer :: iquotient, jquotient, kquotient
    integer :: iremainder, jremainder, kremainder
    logical, dimension(0:2) :: period
    integer, dimension(0:2) :: dims
    logical :: reorder
    integer :: ierr

    if(rank .eq. master) then
        print *, ''
        print *, '********************'
        print *, 'Domain Decomposition'
        print *, '********************'
        print *, ''
    endif

    if(nproc .ne. inp*jnp*knp) then
        write(*,*) nproc, inp*jnp*knp
        call abort
    endif

    dims(2) = inp
    dims(1) = jnp
    dims(0) = knp

    period(0) = .false.
    period(1) = .false.
    period(2) = .false.

    reorder = .false.

    call MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, reorder, cart_comm, ierr)
    call MPI_Cart_coords(cart_comm, rank, 3, coord, ierr)

    iquotient = inc/inp
    jquotient = jnc/jnp
    kquotient = knc/knp

    iremainder = mod(inc,inp)
    jremainder = mod(jnc,jnp)
    kremainder = mod(knc,knp)

    if(coord(2) .lt. iremainder) then
        ilnc = iquotient + 1
    else
        ilnc = iquotient
    end if

    if(coord(1) .lt. jremainder) then
        jlnc = jquotient + 1
    else
        jlnc = jquotient
    end if

    if(coord(0) .lt. kremainder) then
        klnc = kquotient + 1
    else
        klnc = kquotient
    end if

    imin = ngh
    imax = ilnc + imin - 1
    jmin = ngh
    jmax = jlnc + jmin - 1
    kmin = ngh
    kmax = klnc + kmin - 1

    igmin = ngh + coord(2)*iquotient
    jgmin = ngh + coord(1)*jquotient
    kgmin = ngh + coord(0)*kquotient

    if(coord(2) .lt. iremainder) then
        igmin = igmin + coord(2)
    else
        igmin = igmin + iremainder
    endif

    if(coord(1) .lt. jremainder) then
        jgmin = jgmin + coord(1)
    else
        jgmin = jgmin + jremainder
    endif

    if(coord(0) .lt. kremainder) then
        kgmin = kgmin + coord(0)
    else
        kgmin = kgmin + kremainder
    endif

    igmax = igmin + ilnc - 1
    jgmax = jgmin + jlnc - 1
    kgmax = kgmin + klnc - 1

    if(.not. allocated(x)) allocate(x(imin-ngh:imax+ngh))
    if(.not. allocated(y)) allocate(y(jmin-ngh:jmax+ngh))
    if(.not. allocated(z)) allocate(z(kmin-ngh:kmax+ngh))

    if(.not. allocated(rho)) allocate(rho(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
    if(.not. allocated(rho4)) allocate(rho4(imin:imax,jmin:jmax,kmin:kmax))

    if(.not. allocated(globalindex)) allocate(globalindex(imin-swidth:imax+swidth, &
        jmin-swidth:jmax+swidth,kmin-swidth:kmax+swidth))
    if(.not. allocated(localindex)) allocate(localindex(imin:imax,jmin:jmax,kmin:kmax))
    if(.not. allocated(stencilindex)) allocate(stencilindex(-swidth:+swidth,-swidth:+swidth,-swidth:+swidth))


    if(.not. allocated(Stencil)) allocate(Stencil(0:ilnc*jlnc*klnc-1,0:61))

    if(.not. allocated(bxx)) allocate(bxx(imin-1:imax,jmin:jmax,kmin:kmax))
    if(.not. allocated(bxy)) allocate(bxy(imin-1:imax,jmin-1:jmax+1,kmin:kmax))
    if(.not. allocated(bxz)) allocate(bxz(imin-1:imax,jmin:jmax,kmin-1:kmax+1))

    if(.not. allocated(byx)) allocate(byx(imin-1:imax+1,jmin-1:jmax,kmin:kmax))
    if(.not. allocated(byy)) allocate(byy(imin:imax,jmin-1:jmax,kmin:kmax))
    if(.not. allocated(byz)) allocate(byz(imin:imax,jmin-1:jmax,kmin-1:kmax+1))

    if(.not. allocated(bzx)) allocate(bzx(imin-1:imax+1,jmin:jmax,kmin-1:kmax))
    if(.not. allocated(bzy)) allocate(bzy(imin:imax,jmin-1:jmax+1,kmin-1:kmax))
    if(.not. allocated(bzz)) allocate(bzz(imin:imax,jmin:jmax,kmin-1:kmax))

    if(.not. allocated(Gxy)) allocate(Gxy(imin-1:imax,jmin:jmax,kmin:kmax))
    if(.not. allocated(Gxz)) allocate(Gxz(imin-1:imax,jmin:jmax,kmin:kmax))

    if(.not. allocated(Gyx)) allocate(Gyx(imin:imax,jmin-1:jmax,kmin:kmax))
    if(.not. allocated(Gyz)) allocate(Gyz(imin:imax,jmin-1:jmax,kmin:kmax))

    if(.not. allocated(Gzx)) allocate(Gzx(imin:imax,jmin:jmax,kmin-1:kmax))
    if(.not. allocated(Gzy)) allocate(Gzy(imin:imax,jmin:jmax,kmin-1:kmax))

end subroutine DomainDecomposition
