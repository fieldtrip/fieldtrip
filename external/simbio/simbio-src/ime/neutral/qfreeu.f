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
      INTEGER FUNCTION QFREEU()
C     Bestimmt die naechste freie Fortran-Unit fuer Files
C
C     ****************
C     * STANDARD F77 *
C     ****************
C
C
      INCLUDE 'imelibh'
C
      INTEGER IUNIT
C
C      WRITE  (*,40) IFIMIP
C 40   FORMAT (1X,'IFIMIP          ',I2)
C      WRITE  (*,50) MANFIP
C 50   FORMAT (1X,'MANFIP          ',I2)
      DO 10 IUNIT=IFIMIP, MANFIP
C      	 WRITE  (*,30) IUNIT
C 30   	 FORMAT (1X,'IUNIT          ',I2)
C         WRITE  (*,60) IFIKEN(IUNIT)
C 60   	 FORMAT (1X,'IFIKEN(IUNIT)          ',I2)
         IF (IFIKEN(IUNIT).EQ.2) GOTO 20
 10   CONTINUE
C
C     Hier sollte er normalerweise nicht hinkommen:
C
      STOP '>> Keine freie FORTRAN-UNIT verfuegbar <<'
C
 20   QFREEU=IUNIT
      IFIKEN(IUNIT)=1
      R E T U R N
      END
