 !****************************************************************  
 ! MD_Definition.f90
 ! 
 ! Created on: Jul 23, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_Definition
#define SUCCESS 1

    IMPLICIT NONE
    SAVE

    integer, parameter :: double = kind(1.d0)
    real(KIND=double), parameter :: G = 1.0
    real(KIND=double), parameter :: pi = 3.1415927

CONTAINS

END MODULE MD_Definition

