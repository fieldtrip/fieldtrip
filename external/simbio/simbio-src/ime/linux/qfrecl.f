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
      INTEGER FUNCTION QFRECL (NI, ND, NR)
C
C     Berechnung der Recordlaenge eines Direct-Access-Files, bei dem
C     pro Record NI Integer-, NR Real- und ND Double-Zahlen
C     gespeichert werden sollen.
C
C     ***********************
C     * RECHNER - ABHAENGIG *
C     ***********************
C
C
C     ERLAEUTERUNG IBM3090 UNTER VM/XA
C
C     Die Einheit der Recordlaenge im OPEN-Statement das BYTE
C     Die Wortlaenge von INTEGER und REAL ist bei
C     Standard-F77 4 Byte, also 32 Bit
C     -> 1 INTEGER = Recordlaenge 4
C     Die Wortlaenge von DOUBLE ist bei Standard-F77 8 Byte, also 64 Bit
C     -> 1 DOUBLE = Recordlaenge 8
C
C     ---> QFRECL = 4*INTEGER + 4*REAL + 8*DOUBLE
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C NI        I   I                 Anzahl der INTEGER - Zahlen im Record
C ND        I   I                 Anzahl der DOUBLE  - Zahlen im Record
C NR        I   I                 Anzahl der REAL    - Zahlen im Record
C
C
      INCLUDE 'imelibh'
C
      INTEGER NI, ND, NR
C
      QFRECL = 4*(NI+2*ND+NR)
C
      R E T U R N
      END
