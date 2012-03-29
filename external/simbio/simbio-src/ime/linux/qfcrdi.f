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
      SUBROUTINE QFCRDI (IUNIT, FILNAM, FM, IR, IERR)
C
C     Erstellen eines neuen Direct-Access-Files. Falls er bereits existiert,
C     wird er vor dem Oeffnen geloescht.
C
C     ****************
C     * STANDARD F77 *
C     ****************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C IUNIT     I   I                 FORTRAN-UNIT, der der File zugeordn. wird
C FILNAM    C   I *(*)            Name des Files
C FM        C   I *(*)            Formatierung des Files
C                                 'UN formatted'
C                                 'FO rmatted'
C IR        I   I                 Recordlaenge, die in der RECL-Option des
C                                 OPEN-Statements stehen muss
C IERR      I   O                 FEHLERCODE =1: File existiert schon
C                                                und ist offen
C
      INCLUDE 'imelibh'
C
      CHARACTER*(*)  FILNAM
      CHARACTER*(*)  FM
      INTEGER IR, IUNIT, IERR
C
C  Abfrage, ob der File existiert
C
      CALL QFDEL (FILNAM, IERR)
      IF (IERR .GT. 0) RETURN
C
C
      IF (IERR .EQ. -1)  IERR = 0
C
         CALL QFOPDI (IUNIT, FILNAM, 'NEW', FM, IR, IERR)
C
      R E T U R N
      END
