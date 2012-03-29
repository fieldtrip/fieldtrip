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
      SUBROUTINE QSDEFO (NTTE, FILNA, IERR)
C     Setzt die Ausgabe entweder auf Terminal, oder auf File.
C
C     ****************
C     * STANDARD F77 *
C     ****************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C NTTE      I   I                 UNIT-Nummer fuer die Ausgabe
C FILNA     C   I  *(*)           Name des Ausgabe-Files
C                                 bei FILNA=CHAR(0) wird die Ausgabe
C                                 aufs Terminal gesetzt
C IERR      I   O                 Fehler beim Oeffnen des OUTPUT-Files
C
C
      INCLUDE 'imelibh'
C
      INTEGER NTTE, IERR
      CHARACTER*(*) FILNA
C
      IERR=0
C
C     Wenn vorher Ausgabe ueber Datei, dann geoeffnete
C     Ausgabe-Datei schliessen.
C
      IF (.NOT.TERM) CALL QFCLOS (NTTE, 0)
C
C     Wenn Ausgabe aufs Terminal
      IF (FILNA.EQ.CHAR(0)) THEN
        TERM=.TRUE.
        NTTERM=ISTOUT
        NTTE=NTTERM
      ELSE
C
C     Ausgabe auf File
        TERM=.FALSE.
        NTTERM=NTTE
        CALL QFCRSE (NTTE,FILNA,'FO',IERR)
      END IF
C
      R E T U R N
      END
