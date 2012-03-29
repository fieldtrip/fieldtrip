      SUBROUTINE FOFU02(IBEWER,IBEABL,R,S,LKNEL ,FWE   ,FAL   )         NEU
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
C FOFU02:
C =======
C
C BERECHNUNG DER WERTE UND DER LOKALEN ABLEITUNGEN DER FORMFUNKTIONEN
C DES ISOPARAMETRISCHEN DREIECKS (ELEMENTTYP XXX02)
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C IBEWER    I   I                 1: WERTE DER FORMF. BERECHNEN
C IBEABL    I   I                 1: WERTE DER ABLEITUNGEN DER FORMF. BER.
C R         D   I                 1. LOKALE KOORDINATE DES AUSWERTEPUNKTES
C S         D   I                 2. LOKALE KOORDINATE DES AUSWERTEPUNKTES
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
      ISUMQ=0
      ISUMK=0
      DO 10 IKNEL=4,6
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMQ=ISUMQ+1
      ELSE
        F(IKNEL)=Z0
      END IF
   10 CONTINUE
      DO 20 IKNEL=7,9
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMK=ISUMK+1
      ELSE
        F(IKNEL)=Z0
      END IF
   20 CONTINUE
      U=Z1-R-S
      RMI=Z1-R
      SMI=Z1-S
      UMI=Z1-U
      R13M=Z1D3-R
      S13M=Z1D3-S
      U13M=Z1D3-U
      R23M=Z2D3-R
      S23M=Z2D3-S
      U23M=Z2D3-U
      RS  =R*S
      RQ  =R*R
      SQ  =S*S
      RSU =Z9D2*R*S*U
      DXDR=S-SQ-Z2*RS
      DXDS=R-RQ-Z2*RS
      IF(ISUMQ.EQ.3.AND.ISUMK.EQ.3)THEN
C +++++
C VOLL-KUBISCHES DREIECK
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(1)=Z9D2*R*R13M*R23M
          FWE(2)=Z9D2*S*S13M*S23M
          FWE(3)=Z9D2*U*U13M*U23M
          FWE(4)=Z9D2*R*S*(Z2*R-S)
          FWE(5)=Z9D2*S*U*(Z2*S-U)
          FWE(6)=Z9D2*R*U*(Z2*U-R)
          FWE(7)=Z9D2*R*S*(Z2*S-R)
          FWE(8)=Z9D2*S*U*(Z2*U-S)
          FWE(9)=Z9D2*R*U*(Z2*R-U)
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,3)=Z9D2*(-Z11D9-Z3*RQ+Z4 *R-Z6 *RS+Z4 *S-Z3*SQ)
          FAL(2,3)=FAL(1,3)
          FAL(1,1)=Z9D2*( Z2D9 +Z3*RQ-Z2 *R                   )
          FAL(2,1)=Z0
          FAL(1,2)=Z0
          FAL(2,2)=Z9D2*( Z2D9                    -Z2 *S+Z3*SQ)
          FAL(1,6)=Z9D2*( Z2   +Z9*RQ-Z10*R+Z10*RS-Z4 *S+Z2*SQ)
          FAL(2,6)=Z9D2*(       Z5*RQ-Z4 *R+Z4 *RS            )
          FAL(1,4)=Z9D2*(                   Z4 *RS      -   SQ)
          FAL(2,4)=Z9D2*(       Z2*RQ      -Z2 *RS            )
          FAL(1,5)=Z9D2*(                  -Z2 *RS+Z2 *S-Z4*SQ)
          FAL(2,5)=Z9D2*(-Z1   -   RQ+Z2 *R-Z8 *RS+Z8 *S-Z9*SQ)
          FAL(1,9)=Z9D2*(-Z1   -Z9*RQ+Z8 *R-Z8 *RS+Z2 *S-   SQ)
          FAL(2,9)=Z9D2*(      -Z4*RQ+Z2 *R-Z2 *RS            )
          FAL(1,7)=Z9D2*(                  -Z2 *RS      +Z2*SQ)
          FAL(2,7)=Z9D2*(      -   RQ      +Z4 *RS            )
          FAL(1,8)=Z9D2*(                   Z4 *RS-Z4 *S+Z5*SQ)
          FAL(2,8)=Z9D2*( Z2   +Z2*RQ-Z4 *R+Z10*RS-Z10*S+Z9*SQ)
        END IF
      ELSE IF(ISUMQ.EQ.3.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-QUADRATISCHES DREIECK
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(1)=R*(Z2*R-Z1)
          FWE(2)=S*(Z2*S-Z1)
          FWE(3)=U*(Z2*U-Z1)
          FWE(4)=Z4*R*S
          FWE(5)=Z4*S*U
          FWE(6)=Z4*R*U
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,1)= Z4*R-Z1
          FAL(2,1)= Z0
          FAL(1,2)= Z0
          FAL(2,2)= Z4*S-Z1
          FAL(1,3)=-Z4*U+Z1
          FAL(2,3)=-Z4*U+Z1
          FAL(1,4)= Z4*S
          FAL(2,4)= Z4*R
          FAL(1,5)=-Z4*S
          FAL(2,5)=Z4*(U-S)
          FAL(1,6)=Z4*(U-R)
          FAL(2,6)=-Z4*R
        END IF
      ELSE IF(ISUMQ.EQ.0.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-LINEARES DREIECK
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(1)=R
          FWE(2)=S
          FWE(3)=U
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,1)=Z1
          FAL(2,1)=Z0
          FAL(1,2)=Z0
          FAL(2,2)=Z1
          FAL(1,3)=-Z1
          FAL(2,3)=-Z1
        END IF
      ELSE
C +++++
C DREIECK MIT VARIABLER KNOTENZAHL
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(7)=F(7)*Z9D2*R*S*(Z2*S-R)
          FWE(8)=F(8)*Z9D2*S*U*(Z2*U-S)
          FWE(9)=F(9)*Z9D2*R*U*(Z2*R-U)
          FWE(4)=F(4)*(Z4*R*S + F(7)*(Z05*R*S-FWE(7)-RSU) )
          FWE(5)=F(5)*(Z4*S*U + F(8)*(Z05*S*U-FWE(8)-RSU) )
          FWE(6)=F(6)*(Z4*R*U + F(9)*(Z05*R*U-FWE(9)-RSU) )
          FWE(1)=R - F(4)*Z2*R*S - F(6)*Z2*R*U
     *             - F(7)*(Z2D9*FWE(4)-Z1D9*FWE(7)-Z4D9*RSU)
     *             - F(9)*(Z2D9*FWE(9)-Z1D9*FWE(6)-Z4D9*RSU)
          FWE(2)=S - F(5)*Z2*S*U - F(4)*Z2*R*S
     *             - F(7)*(Z2D9*FWE(7)-Z1D9*FWE(4)-Z4D9*RSU)
     *             - F(8)*(Z2D9*FWE(5)-Z1D9*FWE(8)-Z4D9*RSU)
          FWE(3)=U - F(6)*Z2*R*U - F(5)*Z2*S*U
     *             - F(8)*(Z2D9*FWE(8)-Z1D9*FWE(5)-Z4D9*RSU)
     *             - F(9)*(Z2D9*FWE(6)-Z1D9*FWE(9)-Z4D9*RSU)
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,9)= F(9)*Z9D2*(-Z1   -Z9*RQ+Z8 *R-Z8 *RS+Z2 *S-   SQ)
          FAL(2,9)= F(9)*Z9D2*(      -Z4*RQ+Z2 *R-Z2 *RS            )
          FAL(1,7)= F(7)*Z9D2*(                  -Z2 *RS      +Z2*SQ)
          FAL(2,7)= F(7)*Z9D2*(      -   RQ      +Z4 *RS            )
          FAL(1,8)= F(8)*Z9D2*(                   Z4 *RS-Z4 *S+Z5*SQ)
          FAL(2,8)= F(8)*Z9D2*( Z2   +Z2*RQ-Z4 *R+Z10*RS-Z10*S+Z9*SQ)
          FAL(1,4)= F(4)*( Z4*S    +F(7)*( Z05*S    -FAL(1,7)-DXDR) )
          FAL(2,4)= F(4)*( Z4*R    +F(7)*( Z05*R    -FAL(2,7)-DXDS) )
          FAL(1,5)= F(5)*(-Z4*S    +F(8)*(-Z05*S    -FAL(1,8)-DXDR) )
          FAL(2,5)= F(5)*( Z4*(U-S)+F(8)*( Z05*(U-S)-FAL(2,8)-DXDS) )
          FAL(1,6)= F(6)*( Z4*(U-R)+F(9)*( Z05*(U-R)-FAL(1,9)-DXDR) )
          FAL(2,6)= F(6)*(-Z4*R    +F(9)*(-Z05*R    -FAL(2,9)-DXDS) )
          FAL(1,1)= Z1-F(4)*Z2*S    -F(6)*Z2*(U-R)
     *             -F(7)*(Z2D9*FAL(1,4)-Z1D9*FAL(1,7)-Z4D9*DXDR)
     *             -F(9)*(Z2D9*FAL(1,9)-Z1D9*FAL(1,6)-Z4D9*DXDR)
          FAL(2,1)=   -F(4)*Z2*R    +F(6)*Z2*R
     *             -F(7)*(Z2D9*FAL(2,4)-Z1D9*FAL(2,7)-Z4D9*DXDS)
     *             -F(9)*(Z2D9*FAL(2,9)-Z1D9*FAL(2,6)-Z4D9*DXDS)
          FAL(1,2)=    F(5)*Z2*S    -F(4)*Z2*S
     *             -F(7)*(Z2D9*FAL(1,7)-Z1D9*FAL(1,4)-Z4D9*DXDR)
     *             -F(8)*(Z2D9*FAL(1,5)-Z1D9*FAL(1,8)-Z4D9*DXDR)
          FAL(2,2)= Z1-F(5)*Z2*(U-S)-F(4)*Z2*R
     *             -F(7)*(Z2D9*FAL(2,7)-Z1D9*FAL(2,4)-Z4D9*DXDS)
     *             -F(8)*(Z2D9*FAL(2,5)-Z1D9*FAL(2,8)-Z4D9*DXDS)
          FAL(1,3)=-Z1-F(6)*Z2*(U-R)+F(5)*Z2*S
     *             -F(8)*(Z2D9*FAL(1,8)-Z1D9*FAL(1,5)-Z4D9*DXDR)
     *             -F(9)*(Z2D9*FAL(1,6)-Z1D9*FAL(1,9)-Z4D9*DXDR)
          FAL(2,3)=-Z1+F(6)*Z2*R    -F(5)*Z2*(U-S)
     *             -F(8)*(Z2D9*FAL(2,8)-Z1D9*FAL(2,5)-Z4D9*DXDS)
     *             -F(9)*(Z2D9*FAL(2,6)-Z1D9*FAL(2,9)-Z4D9*DXDS)
        END IF
      END IF
      R E T U R N
      END
