      SUBROUTINE GETADJ(IELTYP,ADJ1,ADJ2,ADJ3,IERR)                     NEU
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
C
C#
C
C GETADJ:
C =======
C
C HOLT DIE BEI DER BERECHNUNG DER KINEMATISCH AEQUIVALENTEN KNOTENKRAEFTE
C AUS DEN KNOTENDRUECKEN BENOETIGTEN ADJUNKTEN AUS DEM COMMON-BLOCK UND
C UEBERGIBT SIE ANS RUFENDE PROGRAMM
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C IELTYP    I   I                 TYP DES ELEMENTES ABCDE
C                                 A: ELEMENTKENNUNG
C                                    0 - NORMALES ELEMENT
C                                    1 - RANDELEMENT
C                                 B: =0, FUER SPAETERE VERWENDUNG
C                                 C: GLOBALE DIMENSION DES KOORDINATENRAUMES
C                                 D: LAUFENDE NUMMER DES TYPS
C                                 E: LOKALE DIMENSION DES ELEMENTES
C ADJ1
C ADJ2      D   O                 ADJUNKTEN
C ADJ3  
C IERR      I   O                 FEHLERPARAMETER
C                                 0: ALLES KLAR
C                                 1: FALSCHER ELEMENTTYP
C
C#
C
      INCLUDE 'fofulib.inc'
C
      CHARACTER*6 CHTYP
C
C +++++
C ELEMENTTYP DEKODIEREN
C +++++
      WRITE(CHTYP,'(I5)')IELTYP
      READ(CHTYP,'(BN,2X,3I1)')IGLDIM,IELNUM,ILODIM
      IERR=0
      IF(IGLDIM.EQ.2)THEN
C +++++
C 2D-RECHNUNG
C +++++
        ADJ1=ADJ21
        ADJ2=ADJ22
      ELSE IF(IGLDIM.EQ.3)THEN
C +++++
C 3D-RECHNUNG
C +++++
        ADJ1=ADJ31
        ADJ2=ADJ32
        ADJ3=ADJ33
      ELSE
C +++++
C FALSCHER ELEMENTTYP
C +++++
        IERR=1
        R E T U R N
      END IF
      R E T U R N
      END
