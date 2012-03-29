!
!---->---1---------2---------3---------4---------5---------6---------7--<
!
!    NeuroFEM license:
!    =================
!    Copyright (c) 2007 by 
!    Dr.Carsten Wolters, Dr.Alfred Anwander,  
!    Prof.Dr.H.Buchner, Prof.Dr.G.Knoll, Dr.Adrian Rienaecker, 
!    Rainer Beckmann. 
!
!    Permission is hereby granted, free of charge, to any person 
!    obtaining a copy of the NeuroFEM software and associated 
!    documentation files (the "Software"), to deal in the Software 
!    without restrictions, including without limitations the rights to 
!    use, copy, modify, merge, publish, distribute, sublicense, 
!    and/or sell copies of the Software, and to permit persons to 
!    whom the Software is furnished to do so, subject to the following 
!    conditions:
!
!    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
!    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
!    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
!    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
!    OR OTHER DEALINGS IN THE SOFTWARE. THE ENTIRE RISK AS TO THE 
!    QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH YOU.
!
!    The above copyright notice, permission and warranty notices 
!    shall be included in all copies of the Software. A copyright 
!    notice of the contributing authors and the above permission 
!    and warranty notices shall be included in all source code 
!    files of the Software. 
!
!---->---1---------2---------3---------4---------5---------6---------7--<
      SUBROUTINE QFDEL (FILNAM, IERR)
C     QFDEL loescht den File mit dem Namen FILNAM
C
C     ***********************
C     * RECHNER - ABHAENGIG *
C     ***********************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C FILNAM    C   I *(*)            Name des Files
C IERR      I   O           =  1  File ist offen, kann nicht geloescht werden
C                           = -1  File existiert nicht
C
      CHARACTER*(*) FILNAM
      CHARACTER*83 COMMAND
      INTEGER IERR, LEFI, QCLEN
      LOGICAL LEX, OP
C
      IERR = 0
      LEFI=QCLEN(FILNAM)
      INQUIRE (FILE=FILNAM(1:LEFI), OPENED=OP, EXIST=LEX)
      IF (.NOT.LEX) THEN
         IERR=-1
         R E T U R N
      ENDIF
      IF (OP) THEN
         IERR=1
         R E T U R N
      ENDIF
C
      COMMAND= 'rm '//FILNAM
      CALL SYSTEM (COMMAND(1:LEFI+3))
C
      R E T U R N
      END
