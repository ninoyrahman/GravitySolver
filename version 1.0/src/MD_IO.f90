 !****************************************************************  
 ! MD_IO.f90
 ! 
 ! Created on: Jul 27, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_IO
    USE MD_Parameter
    USE MD_CordinateTransform
    USE MD_Quantity
    IMPLICIT NONE

CONTAINS

    ! this routine is written by pedro montero
    SUBROUTINE Read_parameter(name,value)
        IMPLICIT NONE
        CHARACTER(len=*),INTENT(IN) :: name
        CHARACTER(len=*),INTENT(OUT) :: value
        INTEGER :: i,j,unit,ierror,m,n
        CHARACTER(len=40) :: nametmp
        CHARACTER(len=100) :: line
        LOGICAL ::  r


        i = index(name,'=')


        IF (i.GT.0) THEN
            nametmp = name(1:i-1)
        ELSE
            nametmp = name
        ENDIF



        INQUIRE(unit=15,opened=r)

        IF(.NOT.r) THEN
            unit = 15
        ELSE
            WRITE(*,*) 'the unit 15 is already opened'
            stop
        ENDIF

        OPEN(unit=unit,file='Parameter.dat',status='old')

10      format(A)

15  continue

    READ(unit,10,IOSTAT=ierror) line


    IF (ierror.GT.0) THEN
        WRITE(*,*) 'Error reading the parameter file'
        stop
    ELSE IF(ierror.EQ.-1) THEN
        GOTO 20
    ELSE
    continue
END IF

i = index(line,'=')
j = index(line,'#')

!*********************************************************************
!    Comment line or no equal sign, then go on to the next line      *
!*********************************************************************
if(i .eq. 0 .or. j .eq. 1) goto 15

!***************************************************************
!    Equal sign in a comment, then go on to the next line      *
!***************************************************************
if(j .gt. 0 .and. j .lt. i) goto 15

!***************************************************
!     There is a comment in the line, make sure    *
!     we do not include it as part of the value    *
!***************************************************


m = len(line)

if(j.gt.0) m = j-1

n=0
17 continue
   n=n+1
   if(line(n:n) .eq. ' ') goto 17
   if(nametmp .eq. line(n:i-1))then
       value=line(i+1:m)
       goto 20
   endif
   goto 15
20 continue
   close(unit)
   !      WRITE(*,*) nametmp,value
   END SUBROUTINE Read_parameter

   SUBROUTINE Read_real_parm(nam,val)
       USE MD_Definition
       IMPLICIT NONE
       character(len=*),INTENT(IN):: nam
       real(KIND=double),INTENT(OUT):: val
       character(len=100) str_val
       call Read_parameter(nam,str_val)

       if(index(str_val,'.') .eq. 0)then
           print *,'need decimal point in ',nam
           stop
       endif
       read(str_val,*) val

   END SUBROUTINE Read_real_parm



   SUBROUTINE Read_integer_parm(nam,val)
       IMPLICIT NONE
       character(len=*),INTENT(IN):: nam
       integer,INTENT(OUT):: val
       character(len=20) str_val

       call Read_parameter(nam,str_val)

       read(str_val,*) val

   END SUBROUTINE Read_integer_parm


   subroutine ReadParameterFile
       implicit none

       print *, ''
       print *, '********************'
       print *, 'Read Parameter File'
       print *, '********************'
       print *, ''

       call Read_integer_parm('coordinate=',coordinate)
       call Read_integer_parm('scenerio=',scenerio)
       call Read_integer_parm('test=',test)

   end subroutine ReadParameterFile

   subroutine ViewInfo
       use MD_Parameter
       use MD_CordinateTransform
       implicit none

       real(kind=double), dimension(0:2) :: x
       character(len=100) :: format
       integer :: nghi,nghj,nghk

       open(UNIT=1, FILE='data/Param.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='REWIND')

       print *, ''
       print *, '********'
       print *, 'Setting'
       print *, '********'
       print *, ''

       if(coordinate .eq. 0) write(*,*) 'coordinate: cartesian'
       if(coordinate .eq. 0) write(1,*) 'coordinate: cartesian'
       if(coordinate .eq. 1) write(*,*) 'coordinate: nonuniformly scaled cartesian'
       if(coordinate .eq. 1) write(1,*) 'coordinate: nonuniformly scaled cartesian'

       format = "(A12, I5, I5, I5)"
       write (*,format) 'resolution=',inc,jnc,knc
       write (1,format) 'resolution=',inc,jnc,knc

       format = "(A6, F10.3, F10.3, F10.3)"
       call GetCordinate(imin,jmin,kmin,x)
       write (*,format) 'xmin=',x(0),x(1),x(2)
       write (1,format) 'xmin=',x(0),x(1),x(2)

       call GetCordinate(imax,jmax,kmax,x)
       write (*,format) 'xmax=',x(0),x(1),x(2)
       write (1,format) 'xmax=',x(0),x(1),x(2)
       write(1,*) ''


       if(test .ne. 0) then
           write(*,*) 'Running Test'
       else
           if(scenerio .eq. 0) write(*,*) 'scenerio: vaccum'
           if(scenerio .eq. 0) write(1,*) 'scenerio: vaccum'
           if(scenerio .eq. 1) write(*,*) 'scenerio: centred sphere'
           if(scenerio .eq. 1) write(1,*) 'scenerio: centred sphere'
           if(scenerio .eq. 2) write(*,*) 'scenerio: off-centred sphere'
           if(scenerio .eq. 2) write(1,*) 'scenerio: off-centred sphere'

       end if

       format = "(A10, F10.3)"
       write (*,format), 'mass, M =', M
       write (1,format), 'mass, M =', M

       format = "(A21, F10.3, F10.3, F10.3)"
       write (*,format) 'centre of mass, xc =', xc(0),xc(1),xc(2)
       write (1,format) 'centre of mass, xc =', xc(0),xc(1),xc(2)
       write(1,*) ''

       print *, ''
       print *, '********'
       print *, 'Stencil'
       print *, '********'
       print *, ''

       format = "(F10.3)"


       write(1,*) 'stencil'
       do nghi = -swidth, +swidth
           do nghj = -swidth, +swidth
               do nghk = -swidth, +swidth
                   write(*,format,advance='no') S(nghi,nghj,nghk)
                   write(1,format,advance='no') S(nghi,nghj,nghk)
               end do
               write(*,*) ''
               write(1,*) ''
           end do
           write(*,*) ''
           write(1,*) ''
       end do

       close(unit=1)

   end subroutine ViewInfo

   END MODULE MD_IO
