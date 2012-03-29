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
      SUBROUTINE QFCLAL
C     Schliesst alle geoeffneten Files, ausser Tastatur und Terminal
C
C     *******************************
C     * WAHRSCHEINLICH STANDARD F77 *
C     *******************************
C
C
      INCLUDE 'imelibh'
C
      INTEGER NATERM, NATAST, IUN
C
      IF(.NOT.TAST)THEN
         NATAST=ABS(NTTAST)
         IFIKEN(NATAST)=3
      END IF
      IF(.NOT.TERM)THEN
         NATERM=ABS(NTTERM)
         IFIKEN(NATERM)=3
      END IF
      DO 10 IUN=IFIMIP, MANFIP
         IF (IFIKEN(IUN).EQ.0) THEN
            CLOSE(IUN)
            IFIKEN(IUN)=2
         END IF
 10   CONTINUE
      IF (.NOT.TAST) IFIKEN(NATAST)=0
      IF (.NOT.TERM) IFIKEN(NATERM)=0
      R E T U R N
      END
