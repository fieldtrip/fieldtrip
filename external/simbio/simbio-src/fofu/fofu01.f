      SUBROUTINE FOFU01(IBEWER,IBEABL,R     ,LKNEL ,FWE   ,FAL   )      NEU
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
C FOFU01:
C =======
C
C BERECHNUNG DER WERTE UND DER LOKALEN ABLEITUNGEN DER FORMFUNKTIONEN
C DES ISOPARAMETRISCHEN STABES (ELEMENTTYP XXX01)
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C IBEWER    I   I                 1: WERTE DER FORMF. BERECHNEN
C IBEABL    I   I                 1: WERTE DER ABLEITUNGEN DER FORMF. BER.
C R         D   I                 1. LOKALE KOORDINATE DES AUSWERTEPUNKTES
C LKNEL     I   I  (MKNPEP)       KNOTEN-ELEMENT-ZUORDNUNG UNCOMPRESSED
C FWE       D   O  (MKNPEP)       WERTE DER FORMFUNKTIONEN
C FAL       D   O  (3,MKNPEP)     WERTE DER LOK. ABL. DER FORMF.
C
C#
C
      INCLUDE 'fofulib.inc'
C
      DIMENSION LKNEL (  MKNPEP),FWE   (  MKNPEP),FAL   (3,MKNPEP)
C
      IF(LKNEL(3).NE.0.AND.LKNEL(4).NE.0)THEN
C +++++
C VOLL-KUBISCHER STAB
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z1D16*(((-Z9 *R+Z9)*R+Z1 )*R-Z1)
          FWE( 2)=Z1D16*((( Z9 *R+Z9)*R-Z1 )*R-Z1)
          FWE( 3)=Z1D16*((( Z27*R-Z9)*R-Z27)*R+Z9)
          FWE( 4)=Z1D16*(((-Z27*R-Z9)*R+Z27)*R+Z9)
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL( 1, 1)=Z1D16*((-Z27*R+Z18)*R+Z1)
          FAL( 1, 2)=Z1D16*(( Z27*R+Z18)*R-Z1)
          FAL( 1, 3)=Z1D16*(( Z81*R-Z18)*R-Z27)
          FAL( 1, 4)=Z1D16*((-Z81*R-Z18)*R+Z27)
        END IF
      ELSE IF(LKNEL(3).NE.0)THEN
C +++++
C VOLL-QUADRATISCHER STAB
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z05*R*(R-Z1)
          FWE( 2)=Z05*R*(R+Z1)
          FWE( 3)=Z1-R*R
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL( 1, 1)=R-Z05
          FAL( 1, 2)=R+Z05
          FAL( 1, 3)=-Z2*R
        END IF
      ELSE
C +++++
C VOLL-LINEARER STAB
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z05*(Z1-R)
          FWE( 2)=Z05*(Z1+R)
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL( 1, 1)=-Z05
          FAL( 1, 2)= Z05
        END IF
      END IF
      R E T U R N
      END
