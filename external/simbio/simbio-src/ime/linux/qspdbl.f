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
      DOUBLE PRECISION FUNCTION QSPDBL(I)
C
C     Liefert die Rechnergenauigkeit fuer DOUBLE PRECISION -
C     Variablen.
C     In Kommentaren sind fuer verschiedene Rechner die Genauigkeiten
C     angegeben.
C     Die Funktion liefert fuer I = 1, 2, 3 verschiedene Maschinenparameter
C     Es gilt fuer Rechner mit B hoch T Bits und mit dem kleinsten und
C     groessten Exponenten EMIN und EMAX
C
C         DPMPAR(1) = B**(1 - T), Rechnergenauigkeit
C
C         QSPDBL(2) = B**(EMIN - 1), kleinste Zahl
C
C         QSPDBL(3) = B**EMAX*(1 - B**(-T)), groesste Zahl
C
C     Diese Funktion geht urspruenglich zurueck auf:
C     Argonne National Laboratory. Minpack Project. June 1983.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C
C
      INTEGER I
      DOUBLE PRECISION  DMACH(3)
C
C
C.....IBM RS/6000 (AIX)
      DMACH(1) =1.D-16
      DMACH(2) =0.2225074D-307
      DMACH(3) =0.1797693D+308
C
      QSPDBL = DMACH(I)
      RETURN
      END
