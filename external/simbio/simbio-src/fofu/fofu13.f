      SUBROUTINE FOFU13(IBEWER,IBEABL,R,S,T,LKNEL ,FWE   ,FAL   )       NEU
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
C FOFU13:
C =======
C
C BERECHNUNG DER WERTE UND DER LOKALEN ABLEITUNGEN DER FORMFUNKTIONEN
C DES ISOPARAMETRISCHEN KEILS (ELEMENTTYP XXX13)
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C IBEWER    I   I                 1: WERTE DER FORMF. BERECHNEN
C IBEABL    I   I                 1: WERTE DER ABLEITUNGEN DER FORMF. BER.
C R         D   I                 1. LOKALE KOORDINATE DES AUSWERTEPUNKTES
C S         D   I                 2. LOKALE KOORDINATE DES AUSWERTEPUNKTES
C T         D   I                 3. LOKALE KOORDINATE DES AUSWERTEPUNKTES
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
      DO 10 IKNEL=7,15
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMQ=ISUMQ+1
      ELSE
        F(IKNEL)=Z0
      END IF
   10 CONTINUE
      DO 20 IKNEL=16,24
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMK=ISUMK+1
      ELSE
        F(IKNEL)=Z0
      END IF
   20 CONTINUE
      U=Z1-R-S
      TMI=Z05*(Z1-T)
      TPL=Z05*(Z1+T)
      TMQ=Z05*(Z1-T*T)
      RMI=Z1-R
      SMI=Z1-S
      UMI=Z1-U
      RMQ=Z1-R*R
      SMQ=Z1-S*S
      UMQ=Z1-U*U
      R13M=Z1D3-R
      R23M=Z2D3-R
      S13M=Z1D3-S
      S23M=Z2D3-S
      U13M=Z1D3-U
      U23M=Z2D3-U
      T13M=Z1D3-T
      T13P=Z1D3+T
      T23M=Z2D3-T
      T23P=Z2D3+T
      RSU=R*S*U*Z9D2
      IF(ISUMQ.EQ.9.AND.ISUMK.EQ.9)THEN
C +++++
C VOLL-KUBISCHER KEIL
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z9D4*R*TMI*(Z2*R23M*R13M-TMQ)
          FWE( 2)=Z9D4*S*TMI*(Z2*S23M*S13M-TMQ)
          FWE( 3)=Z9D4*U*TMI*(Z2*U23M*U13M-TMQ)
          FWE( 4)=Z9D4*R*TPL*(Z2*R23M*R13M-TMQ)
          FWE( 5)=Z9D4*S*TPL*(Z2*S23M*S13M-TMQ)
          FWE( 6)=Z9D4*U*TPL*(Z2*U23M*U13M-TMQ)
          FWE( 7)=Z9D2*R*S*(Z2*R-S)*TMI
          FWE( 8)=Z9D2*S*U*(Z2*S-U)*TMI
          FWE( 9)=Z9D2*U*R*(Z2*U-R)*TMI
          FWE(10)=Z9D2*R*S*(Z2*R-S)*TPL
          FWE(11)=Z9D2*S*U*(Z2*S-U)*TPL
          FWE(12)=Z9D2*U*R*(Z2*U-R)*TPL
          FWE(13)=Z27D4*R*TMI*TPL*T13M
          FWE(14)=Z27D4*S*TMI*TPL*T13M
          FWE(15)=Z27D4*U*TMI*TPL*T13M
          FWE(16)=Z9D2*R*S*(Z2*S-R)*TMI
          FWE(17)=Z9D2*S*U*(Z2*U-S)*TMI
          FWE(18)=Z9D2*U*R*(Z2*R-U)*TMI
          FWE(19)=Z9D2*R*S*(Z2*S-R)*TPL
          FWE(20)=Z9D2*S*U*(Z2*U-S)*TPL
          FWE(21)=Z9D2*U*R*(Z2*R-U)*TPL
          FWE(22)=Z27D4*R*TMI*TPL*T13P
          FWE(23)=Z27D4*S*TMI*TPL*T13P
          FWE(24)=Z27D4*U*TMI*TPL*T13P
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1, 1)=Z9D8*TMI*(T*T-Z1D9-Z12*R*R23M)
          FAL(2, 1)=Z0
          FAL(3, 1)=Z9D16*R*(Z1D9+Z4*R*RMI+Z3*T*T23M)
          FAL(1, 2)=Z0
          FAL(2, 2)=Z9D8*TMI*(T*T-Z1D9-Z12*S*S23M)
          FAL(3, 2)=Z9D16*S*(Z1D9+Z4*S*SMI+Z3*T*T23M)
          FAL(1, 3)=-Z9D8*TMI*(T*T-Z1D9-Z12*U*U23M)
          FAL(2, 3)=FAL(1, 3)
          FAL(3, 3)=Z9D16*U*(Z1D9+Z4*U*UMI+Z3*T*T23M)
          FAL(1, 4)=Z9D8*TPL*(T*T-Z1D9-Z12*R*R23M)
          FAL(2, 4)=Z0
          FAL(3, 4)=Z9D16*R*(-Z1D9-Z4*R*RMI+Z3*T*T23P)
          FAL(1, 5)=Z0
          FAL(2, 5)=Z9D8*TPL*(T*T-Z1D9-Z12*S*S23M)
          FAL(3, 5)=Z9D16*S*(-Z1D9-Z4*S*SMI+Z3*T*T23P)
          FAL(1, 6)=-Z9D8*TPL*(T*T-Z1D9-Z12*U*U23M)
          FAL(2, 6)=FAL(1, 6)
          FAL(3, 6)=Z9D16*U*(-Z1D9-Z4*U*UMI+Z3*T*T23P)
          FAL(1, 7)=Z9D2*TMI*S*(Z4*R-S)
          FAL(2, 7)=Z9*TMI*R*(R-S)
          FAL(3, 7)=Z9D4*R*S*(S-Z2*R)
          FAL(1, 8)=Z9*S*TMI*(U-S)
          FAL(2, 8)=Z9D2*TMI*(Z8*S*RMI-RMI*RMI-Z9*S*S)
          FAL(3, 8)=Z9D4*S*U*(U-Z2*S)
          FAL(1, 9)=Z9D2*TMI*(Z2*SMI*SMI-Z10*R*SMI+Z9*R*R)
          FAL(2, 9)=Z9D2*TMI*R*(Z5*R-Z4*SMI)
          FAL(3, 9)=Z9D4*U*R*(R-Z2*U)
          FAL(1,10)=Z9D2*TPL*S*(Z4*R-S)
          FAL(2,10)=Z9*TPL*R*(R-S)
          FAL(3,10)=-FAL(3, 7)
          FAL(1,11)=Z9*S*TPL*(U-S)
          FAL(2,11)=Z9D2*TPL*(Z8*S*RMI-RMI*RMI-Z9*S*S)
          FAL(3,11)=-FAL(3, 8)
          FAL(1,12)=Z9D2*TPL*(Z2*SMI*SMI-Z10*R*SMI+Z9*R*R)
          FAL(2,12)=Z9D2*TPL*R*(Z5*R-Z4*SMI)
          FAL(3,12)=-FAL(3, 9)
          FAL(1,13)=Z27D4*TMI*TPL*T13M
          FAL(2,13)=Z0
          FAL(3,13)=Z27D16*R*(-Z1-Z2D3*T+Z3*T*T)
          FAL(1,14)=Z0
          FAL(2,14)=FAL(1,13)
          FAL(3,14)=Z27D16*S*(-Z1-Z2D3*T+Z3*T*T)
          FAL(1,15)=-FAL(1,13)
          FAL(2,15)=-FAL(2,14)
          FAL(3,15)=Z27D16*U*(-Z1-Z2D3*T+Z3*T*T)
          FAL(1,16)=Z9*TMI*S*(S-R)
          FAL(2,16)=Z9D2*TMI*R*(Z4*S-R)
          FAL(3,16)=Z9D4*R*S*(R-Z2*S)
          FAL(1,17)=Z9D2*TMI*S*(Z5*S-Z4*RMI)
          FAL(2,17)=Z9D2*TMI*(Z2*RMI*RMI-Z10*S*RMI+Z9*S*S)
          FAL(3,17)=Z9D4*S*U*(S-Z2*U)
          FAL(1,18)=Z9D2*TMI*(Z8*R*SMI-SMI*SMI-Z9*R*R)
          FAL(2,18)=Z9*R*TMI*(U-R)
          FAL(3,18)=Z9D4*R*U*(U-Z2*R)
          FAL(1,19)=Z9*TPL*S*(S-R)
          FAL(2,19)=Z9D2*TPL*(Z4*R*S-R*R)
          FAL(3,19)=-FAL(3,16)
          FAL(1,20)=Z9D2*TPL*S*(Z5*S-Z4*RMI)
          FAL(2,20)=Z9D2*TPL*(Z2*RMI*RMI-Z10*S*RMI+Z9*S*S)
          FAL(3,20)=-FAL(3,17)
          FAL(1,21)=Z9D2*TPL*(Z8*R*SMI-SMI*SMI-Z9*R*R)
          FAL(2,21)=Z9*R*TPL*(U-R)
          FAL(3,21)=-FAL(3,18)
          FAL(1,22)=Z27D16*TMI*TPL*T13P
          FAL(2,22)=Z0
          FAL(3,22)=Z27D16*R*(Z1-Z2D3*T-Z3*T*T)
          FAL(1,23)=Z0
          FAL(2,23)=FAL(1,22)
          FAL(3,23)=Z27D16*S*(Z1-Z2D3*T-Z3*T*T)
          FAL(1,24)=-FAL(1,22)
          FAL(2,24)=-FAL(2,23)
          FAL(3,24)=Z27D16*U*(Z1-Z2D3*T-Z3*T*T)
        END IF
      ELSE IF(ISUMQ.EQ.9.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-QUADRATISCHER KEIL
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=(Z2*R-Z1)*R*TMI-R*TMQ
          FWE( 2)=(Z2*S-Z1)*S*TMI-S*TMQ
          FWE( 3)=(Z2*U-Z1)*U*TMI-U*TMQ
          FWE( 4)=(Z2*R-Z1)*R*TPL-R*TMQ
          FWE( 5)=(Z2*S-Z1)*S*TPL-S*TMQ
          FWE( 6)=(Z2*U-Z1)*U*TPL-U*TMQ
          FWE( 7)=Z4*R*S*TMI
          FWE( 8)=Z4*S*U*TMI
          FWE( 9)=Z4*R*U*TMI
          FWE(10)=Z4*R*S*TPL
          FWE(11)=Z4*S*U*TPL
          FWE(12)=Z4*R*U*TPL
          FWE(13)=Z2*R*  TMQ
          FWE(14)=Z2*S*  TMQ
          FWE(15)=Z2*U*  TMQ
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1, 1)= TMI*(Z4*R-Z1 )-TMQ
          FAL(2, 1)= Z0
          FAL(3, 1)=-R  *(   R-Z05)+R  *T
          FAL(1, 2)= Z0
          FAL(2, 2)= TMI*(Z4*S-Z1 )-TMQ
          FAL(3, 2)=-S  *(   S-Z05)+S  *T
          FAL(1, 3)=-TMI*(Z4*U-Z1 )+TMQ
          FAL(2, 3)=-TMI*(Z4*U-Z1 )+TMQ
          FAL(3, 3)=-U  *(   U-Z05)+U  *T
          FAL(1, 4)= TPL*(Z4*R-Z1 )-TMQ
          FAL(2, 4)= Z0
          FAL(3, 4)= R  *(   R-Z05)+R  *T
          FAL(1, 5)= Z0
          FAL(2, 5)= TPL*(Z4*S-Z1 )-TMQ
          FAL(3, 5)= S  *(   S-Z05)+S  *T
          FAL(1, 6)=-TPL*(Z4*U-Z1 )+TMQ
          FAL(2, 6)=-TPL*(Z4*U-Z1 )+TMQ
          FAL(3, 6)= U  *(   U-Z05)+U  *T
          FAL(1, 7)= Z4 *    S     *TMI
          FAL(2, 7)= Z4 *    R     *TMI
          FAL(3, 7)=-Z2 *    R     *S
          FAL(1, 8)=-Z4 *    S     *TMI
          FAL(2, 8)= Z4 *(   U-S )*TMI
          FAL(3, 8)=-Z2 *    S     *U
          FAL(1, 9)= Z4 *(   U-R )*TMI
          FAL(2, 9)=-Z4 *    R     *TMI
          FAL(3, 9)=-Z2 *    R     *U
          FAL(1,10)= Z4 *    S     *TPL
          FAL(2,10)= Z4 *    R     *TPL
          FAL(3,10)= Z2 *    R     *S
          FAL(1,11)=-Z4 *    S     *TPL
          FAL(2,11)= Z4 *(   U-S )*TPL
          FAL(3,11)= Z2 *    S     *U
          FAL(1,12)= Z4 *(   U-R )*TPL
          FAL(2,12)=-Z4 *    R     *TPL
          FAL(3,12)= Z2 *    R     *U
          FAL(1,13)= Z2            *TMQ
          FAL(2,13)= Z0
          FAL(3,13)=-Z2 *    R     *T
          FAL(1,14)= Z0
          FAL(2,14)= Z2            *TMQ
          FAL(3,14)=-Z2 *    S     *T
          FAL(1,15)=-Z2            *TMQ
          FAL(2,15)=-Z2            *TMQ
          FAL(3,15)=-Z2 *    U     *T
        END IF
      ELSE IF(ISUMQ.EQ.0.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-LINEARER KEIL
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=R*TMI
          FWE( 2)=S*TMI
          FWE( 3)=U*TMI
          FWE( 4)=R*TPL
          FWE( 5)=S*TPL
          FWE( 6)=U*TPL
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,1)= TMI
          FAL(2,1)= Z0
          FAL(3,1)=-Z05*R
          FAL(1,2)= Z0
          FAL(2,2)= TMI
          FAL(3,2)=-Z05*S
          FAL(1,3)=-TMI
          FAL(2,3)=-TMI
          FAL(3,3)=-Z05*U
          FAL(1,4)= TPL
          FAL(2,4)= Z0
          FAL(3,4)= Z05*R
          FAL(1,5)= Z0
          FAL(2,5)= TPL
          FAL(3,5)= Z05*S
          FAL(1,6)=-TPL
          FAL(2,6)=-TPL
          FAL(3,6)= Z05*U
        END IF
      ELSE
C +++++
C KEIL MIT VARIABLER KNOTENZAHL
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(16)=F(16)*Z9D2*R*S*(Z2*S-R)*TMI
          FWE(17)=F(17)*Z9D2*S*U*(Z2*U-S)*TMI
          FWE(18)=F(18)*Z9D2*U*R*(Z2*R-U)*TMI
          FWE(19)=F(19)*Z9D2*R*S*(Z2*S-R)*TPL
          FWE(20)=F(20)*Z9D2*S*U*(Z2*U-S)*TPL
          FWE(21)=F(21)*Z9D2*U*R*(Z2*R-U)*TPL
          FWE(22)=F(22)*Z27D4*R*TMI*TPL*T13P
          FWE(23)=F(23)*Z27D4*S*TMI*TPL*T13P
          FWE(24)=F(24)*Z27D4*U*TMI*TPL*T13P
          FWE( 7)=F( 7)*(Z4*R*S*TMI+F(16)*(Z05*R*S*TMI - FWE(16)
     *                   - RSU*TMI))
          FWE( 8)=F( 8)*(Z4*S*U*TMI+F(17)*(Z05*S*U*TMI - FWE(17)
     *                   - RSU*TMI))
          FWE( 9)=F( 9)*(Z4*U*R*TMI+F(18)*(Z05*U*R*TMI - FWE(18)
     *                   - RSU*TMI))
          FWE(10)=F(10)*(Z4*R*S*TPL+F(19)*(Z05*R*S*TPL - FWE(19)
     *                   - RSU*TPL))
          FWE(11)=F(11)*(Z4*S*U*TPL+F(20)*(Z05*S*U*TPL - FWE(20)
     *                   - RSU*TPL))
          FWE(12)=F(12)*(Z4*U*R*TPL+F(21)*(Z05*U*R*TPL - FWE(21)
     *                   - RSU*TPL))
          FWE(13)=F(13)*(Z2*R*TMQ + F(22)*(Z1D4*R*TMQ - FWE(22)))
          FWE(14)=F(14)*(Z2*S*TMQ + F(23)*(Z1D4*S*TMQ - FWE(23)))
          FWE(15)=F(15)*(Z2*U*TMQ + F(24)*(Z1D4*U*TMQ - FWE(24)))
          FWE( 1)=R*TMI - F( 7)*Z2*R*S*TMI - F( 9)*Z2*R*U*TMI
     *                  - F(13)*R*TMQ
     *                 - F(16)*(Z2D9*FWE( 7)-Z1D9*FWE(16)-Z4D9*RSU*TMI)
     *                 - F(18)*(Z2D9*FWE(18)-Z1D9*FWE( 9)-Z4D9*RSU*TMI)
     *                 - F(22)*(Z2D9*FWE(13)-Z1D9*FWE(22))
          FWE( 2)=S*TMI - F( 8)*Z2*U*S*TMI - F( 7)*Z2*R*S*TMI
     *                  - F(14)*S*TMQ
     *                 - F(17)*(Z2D9*FWE( 8)-Z1D9*FWE(17)-Z4D9*RSU*TMI)
     *                 - F(16)*(Z2D9*FWE(16)-Z1D9*FWE( 7)-Z4D9*RSU*TMI)
     *                 - F(23)*(Z2D9*FWE(14)-Z1D9*FWE(23))
          FWE( 3)=U*TMI - F( 9)*Z2*R*U*TMI - F( 8)*Z2*U*S*TMI
     *                  - F(15)*U*TMQ
     *                 - F(18)*(Z2D9*FWE( 9)-Z1D9*FWE(18)-Z4D9*RSU*TMI)
     *                 - F(17)*(Z2D9*FWE(17)-Z1D9*FWE( 8)-Z4D9*RSU*TMI)
     *                 - F(24)*(Z2D9*FWE(15)-Z1D9*FWE(24))
          FWE( 4)=R*TPL - F(10)*Z2*R*S*TPL - F(12)*Z2*R*U*TPL
     *                  - F(13)*R*TMQ
     *                 - F(19)*(Z2D9*FWE(10)-Z1D9*FWE(19)-Z4D9*RSU*TPL)
     *                 - F(21)*(Z2D9*FWE(21)-Z1D9*FWE(12)-Z4D9*RSU*TPL)
     *                 - F(22)*(Z2D9*FWE(22)-Z1D9*FWE(13))
          FWE( 5)=S*TPL - F(11)*Z2*S*U*TPL - F(10)*Z2*R*S*TPL
     *                  - F(14)*S*TMQ
     *                 - F(20)*(Z2D9*FWE(11)-Z1D9*FWE(20)-Z4D9*RSU*TPL)
     *                 - F(19)*(Z2D9*FWE(19)-Z1D9*FWE(10)-Z4D9*RSU*TPL)
     *                 - F(23)*(Z2D9*FWE(23)-Z1D9*FWE(14))
          FWE( 6)=U*TPL - F(12)*Z2*R*U*TPL - F(11)*Z2*S*U*TPL
     *                  - F(15)*U*TMQ
     *                 - F(21)*(Z2D9*FWE(12)-Z1D9*FWE(21)-Z4D9*RSU*TPL)
     *                 - F(20)*(Z2D9*FWE(20)-Z1D9*FWE(11)-Z4D9*RSU*TPL)
     *                 - F(24)*(Z2D9*FWE(24)-Z1D9*FWE(15))
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,16)= F(16)*Z9*TMI*S*(S-R)
          FAL(2,16)= F(16)*Z9D2*TMI*R*(Z4*S-R)
          FAL(3,16)= F(16)*Z9D4*R*S*(R-Z2*S)
          FAL(1,17)= F(17)*Z9D2*TMI*S*(Z5*S-Z4*RMI)
          FAL(2,17)= F(17)*Z9D2*TMI*(Z2*RMI*RMI-Z10*S*RMI+Z9*S*S)
          FAL(3,17)= F(17)*Z9D4*S*U*(S-Z2*U)
          FAL(1,18)= F(18)*Z9D2*TMI*(Z8*R*SMI-SMI*SMI-Z9*R*R)
          FAL(2,18)= F(18)*Z9*R*TMI*(U-R)
          FAL(3,18)= F(18)*Z9D4*R*U*(U-Z2*R)
          FAL(1,19)= F(19)*Z9*TPL*S*(S-R)
          FAL(2,19)= F(19)*Z9D2*TPL*(Z4*R*S-R*R)
          FAL(3,19)= F(19)*Z9D4*R*S*(Z2*S-R)
          FAL(1,20)= F(20)*Z9D2*TPL*S*(Z5*S-Z4*RMI)
          FAL(2,20)= F(20)*Z9D2*TPL*(Z2*RMI*RMI-Z10*S*RMI+Z9*S*S)
          FAL(3,20)= F(20)*Z9D4*S*U*(Z2*U-S)
          FAL(1,21)= F(21)*Z9D2*TPL*(Z8*R*SMI-SMI*SMI-Z9*R*R)
          FAL(2,21)= F(21)*Z9*R*TPL*(U-R)
          FAL(3,21)= F(21)*Z9D4*R*U*(Z2*R-U)
          FAL(1,22)= F(22)*Z27D16*TMI*TPL*T13P
          FAL(2,22)=       Z0
          FAL(3,22)= F(22)*Z27D16*R*(Z1-Z2D3*T-Z3*T*T)
          FAL(1,23)=       Z0
          FAL(2,23)= F(23)*Z27D16*TMI*TPL*T13P
          FAL(3,23)= F(23)*Z27D16*S*(Z1-Z2D3*T-Z3*T*T)
          FAL(1,24)=-F(24)*Z27D16*TMI*TPL*T13P
          FAL(2,24)=        FAL(1,24)
          FAL(3,24)= F(24)*Z27D16*U*(Z1-Z2D3*T-Z3*T*T)
          FAL(1, 7)= F( 7)*(Z4*S*TMI+F(16)*(Z05*S*TMI-FAL(1,16)
     *                                     -Z9D2*S*(U-R)*TMI))
          FAL(2, 7)= F( 7)*(Z4*R*TMI+F(16)*(Z05*R*TMI-FAL(2,16)
     *                                     -Z9D2*R*(U-S)*TMI))
          FAL(3, 7)=-F( 7)*(Z2*R*S  +F(16)*(Z025*R*S +FAL(3,16)
     *                                     -Z05*RSU))
          FAL(1, 8)=-F( 8)*(Z4*S*TMI+F(17)*(Z05*S*TMI+FAL(1,17)
     *                                     +Z9D2*S*(U-R)*TMI))
          FAL(2, 8)= F( 8)*(Z4*(U-S)*TMI+F(17)*(Z05*(U-S)*TMI
     *                           -FAL(2,17)-Z9D2*R*(U-S)*TMI))
          FAL(3, 8)=-F( 8)*(Z2*S*U  +F(17)*(Z025*S*U +FAL(3,17)
     *                                     -Z05*RSU))
          FAL(1, 9)= F( 9)*(Z4*(U-R)*TMI+F(18)*(Z05*(U-R)*TMI
     *                           -FAL(1,18)-Z9D2*S*(U-R)*TMI))
          FAL(2, 9)=-F( 9)*(Z4*R*TMI+F(18)*(Z05*R*TMI+FAL(2,18)
     *                                     +Z9D2*R*(U-S)*TMI))
          FAL(3, 9)=-F( 9)*(Z2*U*R  +F(18)*(Z025*U*R +FAL(3,18)
     *                                     -Z05*RSU))
          FAL(1,10)= F(10)*(Z4*S*TPL+F(19)*(Z05*S*TPL-FAL(1,19)
     *                                     -Z9D2*S*(U-R)*TPL))
          FAL(2,10)= F(10)*(Z4*R*TPL+F(19)*(Z05*R*TPL-FAL(2,19)
     *                                     -Z9D2*R*(U-S)*TPL))
          FAL(3,10)= F(10)*(Z2*R*S  +F(19)*(Z025*R*S -FAL(3,19)
     *                                     -Z05*RSU))
          FAL(1,11)=-F(11)*(Z4*S*TPL+F(20)*(Z05*S*TPL+FAL(1,20)
     *                                     +Z9D2*S*(U-R)*TPL))
          FAL(2,11)= F(11)*(Z4*(U-S)*TPL+F(20)*(Z05*(U-S)*TPL
     *                           -FAL(2,20)-Z9D2*R*(U-S)*TPL))
          FAL(3,11)= F(11)*(Z2*S*U  +F(20)*(Z025*S*U -FAL(3,20)
     *                                     -Z05*RSU))
          FAL(1,12)= F(12)*(Z4*(U-R)*TPL+F(21)*(Z05*(U-R)*TPL
     *                           -FAL(1,21)-Z9D2*S*(U-R)*TPL))
          FAL(2,12)=-F(12)*(Z4*R*TPL+F(21)*(Z05*R*TPL+FAL(2,21)
     *                                     +Z9D2*R*(U-S)*TPL))
          FAL(3,12)= F(12)*(Z2*U*R  +F(21)*(Z025*U*R -FAL(3,21)
     *                                     -Z05*RSU))
          FAL(1,13)= F(13)*Z2      *TMQ
          FAL(1,13)= F(13)*(Z2*TMQ + F(22)*(Z1D4*TMQ - FAL(1,22)))
          FAL(2,13)=        Z0
          FAL(3,13)=-F(13)*(Z2*R*T + F(22)*(Z1D4*R*T + FAL(3,22)))
          FAL(1,14)=        Z0
          FAL(2,14)= F(14)*(Z2*TMQ + F(23)*(Z1D4*TMQ - FAL(2,23)))
          FAL(3,14)=-F(14)*(Z2*S*T + F(23)*(Z1D4*S*T + FAL(3,23)))
          FAL(1,15)=-F(15)*(Z2*TMQ + F(24)*(Z1D4*TMQ + FAL(1,24)))
          FAL(2,15)=        FAL(1,15)
          FAL(3,15)=-F(15)*(Z2*U*T + F(24)*(Z1D4*U*T + FAL(3,24)))
          FAL(1, 1)=TMI - F( 7)*Z2*S*TMI - F( 9)*Z2*(U-R)*TMI
     *                  - F(13)*TMQ
     *            - F(16)*(Z2D9*FAL(1, 7)-Z1D9*FAL(1,16)-Z2*S*(U-R)*TMI)
     *            - F(18)*(Z2D9*FAL(1,18)-Z1D9*FAL(1, 9)-Z2*S*(U-R)*TMI)
     *            - F(22)*(Z2D9*FAL(1,13)-Z1D9*FAL(1,22))
          FAL(2, 1)=    - F( 7)*Z2*R*TMI + F( 9)*Z2*R*TMI
     *            - F(16)*(Z2D9*FAL(2, 7)-Z1D9*FAL(2,16)-Z2*R*(U-S)*TMI)
     *            - F(18)*(Z2D9*FAL(2,18)-Z1D9*FAL(2, 9)-Z2*R*(U-S)*TMI)
     *            - F(22)*(Z2D9*FAL(2,13)-Z1D9*FAL(2,22))
          FAL(3, 1)=-Z05*R + F( 7)*R*S + F( 9)*R*U + F(13)*R*T
     *            - F(16)*(Z2D9*FAL(3, 7)-Z1D9*FAL(3,16)+Z2D9*RSU)
     *            - F(18)*(Z2D9*FAL(3,18)-Z1D9*FAL(3, 9)+Z2D9*RSU)
     *            - F(22)*(Z2D9*FAL(3,13)-Z1D9*FAL(3,22))
          FAL(1, 2)=         F( 8)*Z2*S*TMI - F( 7)*Z2*S*TMI
     *            - F(17)*(Z2D9*FAL(1, 8)-Z1D9*FAL(1,17)-Z2*S*(U-R)*TMI)
     *            - F(16)*(Z2D9*FAL(1,16)-Z1D9*FAL(1, 7)-Z2*S*(U-R)*TMI)
     *            - F(23)*(Z2D9*FAL(1,14)-Z1D9*FAL(1,23))
          FAL(2, 2)= TMI - F( 8)*Z2*(U-S)*TMI - F( 7)*Z2*R*TMI
     *                   - F(14)*TMQ
     *           - F(17)*(Z2D9*FAL(2, 8)-Z1D9*FAL(2,17)-Z2*R*(U-S)*TMI)
     *           - F(16)*(Z2D9*FAL(2,16)-Z1D9*FAL(2, 7)-Z2*R*(U-S)*TMI)
     *           - F(23)*(Z2D9*FAL(2,14)-Z1D9*FAL(2,23))
          FAL(3, 2)=-Z05*S + F( 8)*U*S + F( 7)*R*S + F(14)*S*T
     *           - F(17)*(Z2D9*FAL(3, 8)-Z1D9*FAL(3,17)+Z2D9*RSU)
     *           - F(16)*(Z2D9*FAL(3,16)-Z1D9*FAL(3, 7)+Z2D9*RSU)
     *           - F(23)*(Z2D9*FAL(3,14)-Z1D9*FAL(3,23))
          FAL(1, 3)=-TMI - F( 9)*Z2*(U-R)*TMI + F( 8)*Z2*S*TMI
     *                   + F(15)*TMQ
     *            - F(18)*(Z2D9*FAL(1, 9)-Z1D9*FAL(1,18)-Z2*S*(U-R)*TMI)
     *            - F(17)*(Z2D9*FAL(1,17)-Z1D9*FAL(1, 8)-Z2*S*(U-R)*TMI)
     *            - F(24)*(Z2D9*FAL(1,15)-Z1D9*FAL(1,24))
          FAL(2, 3)=-TMI + F( 9)*Z2*R*TMI - F( 8)*Z2*(U-S)*TMI
     *                   + F(15)*TMQ
     *            - F(18)*(Z2D9*FAL(2, 9)-Z1D9*FAL(2,18)-Z2*R*(U-S)*TMI)
     *            - F(17)*(Z2D9*FAL(2,17)-Z1D9*FAL(2, 8)-Z2*R*(U-S)*TMI)
     *            - F(24)*(Z2D9*FAL(2,15)-Z1D9*FAL(2,24))
          FAL(3, 3)=-Z05*U+ F( 9)*R*U + F( 8)*U*S + F(15)*U*T
     *            - F(18)*(Z2D9*FAL(3, 9)-Z1D9*FAL(3,18)+Z2D9*RSU)
     *            - F(17)*(Z2D9*FAL(3,17)-Z1D9*FAL(3, 8)+Z2D9*RSU)
     *            - F(24)*(Z2D9*FAL(3,15)-Z1D9*FAL(3,24))
          FAL(1, 4)= TPL - F(10)*Z2*S*TPL - F(12)*Z2*(U-R)*TPL
     *                   - F(13)*TMQ
     *            - F(19)*(Z2D9*FAL(1,10)-Z1D9*FAL(1,19)-Z2*S*(U-R)*TPL)
     *            - F(21)*(Z2D9*FAL(1,21)-Z1D9*FAL(1,12)-Z2*S*(U-R)*TPL)
     *            - F(22)*(Z2D9*FAL(1,22)-Z1D9*FAL(1,13))
          FAL(2, 4)=     - F(10)*Z2*R*TPL + F(12)*Z2*R*TPL
     *            - F(19)*(Z2D9*FAL(2,10)-Z1D9*FAL(2,19)-Z2*R*(U-S)*TPL)
     *            - F(21)*(Z2D9*FAL(2,21)-Z1D9*FAL(2,12)-Z2*R*(U-S)*TPL)
     *            - F(22)*(Z2D9*FAL(2,22)-Z1D9*FAL(2,13))
          FAL(3, 4)= Z05*R - F(10)*R*S      - F(12)*R*U + F(13)*R*T
     *            - F(19)*(Z2D9*FAL(3,10)-Z1D9*FAL(3,19)-Z2D9*RSU)
     *            - F(21)*(Z2D9*FAL(3,21)-Z1D9*FAL(3,12)-Z2D9*RSU)
     *            - F(22)*(Z2D9*FAL(3,22)-Z1D9*FAL(3,13))
          FAL(1, 5)=       F(11)*Z2*S*TPL - F(10)*Z2*S*TPL
     *            - F(20)*(Z2D9*FAL(1,11)-Z1D9*FAL(1,20)-Z2*S*(U-R)*TPL)
     *            - F(19)*(Z2D9*FAL(1,19)-Z1D9*FAL(1,10)-Z2*S*(U-R)*TPL)
     *            - F(23)*(Z2D9*FAL(1,23)-Z1D9*FAL(1,14))
          FAL(2, 5)= TPL - F(11)*Z2*(U-S)*TPL - F(10)*Z2*R*TPL
     *                   - F(14)*TMQ
     *            - F(20)*(Z2D9*FAL(2,11)-Z1D9*FAL(2,20)-Z2*R*(U-S)*TPL)
     *            - F(19)*(Z2D9*FAL(2,19)-Z1D9*FAL(2,10)-Z2*R*(U-S)*TPL)
     *            - F(23)*(Z2D9*FAL(2,23)-Z1D9*FAL(2,14))
          FAL(3, 5)= Z05*S - F(11)*S*U - F(10)*R*S + F(14)*S*T
     *            - F(20)*(Z2D9*FAL(3,11)-Z1D9*FAL(3,20)-Z2D9*RSU)
     *            - F(19)*(Z2D9*FAL(3,19)-Z1D9*FAL(3,10)-Z2D9*RSU)
     *            - F(23)*(Z2D9*FAL(3,23)-Z1D9*FAL(3,14))
          FAL(1, 6)=-TPL - F(12)*Z2*(U-R)*TPL + F(11)*Z2*S*TPL
     *                   + F(15)*TMQ
     *            - F(21)*(Z2D9*FAL(1,12)-Z1D9*FAL(1,21)-Z2*S*(U-R)*TPL)
     *            - F(20)*(Z2D9*FAL(1,20)-Z1D9*FAL(1,11)-Z2*S*(U-R)*TPL)
     *            - F(24)*(Z2D9*FAL(1,24)-Z1D9*FAL(1,15))
          FAL(2, 6)=-TPL + F(12)*Z2*R*TPL - F(11)*Z2*(U-S)*TPL
     *                   + F(15)*TMQ
     *            - F(21)*(Z2D9*FAL(2,12)-Z1D9*FAL(2,21)-Z2*R*(U-S)*TPL)
     *            - F(20)*(Z2D9*FAL(2,20)-Z1D9*FAL(2,11)-Z2*R*(U-S)*TPL)
     *            - F(24)*(Z2D9*FAL(2,24)-Z1D9*FAL(2,15))
          FAL(3, 6)= Z05*U - F(12)*R*U - F(11)*S*U + F(15)*U*T
     *            - F(21)*(Z2D9*FAL(3,12)-Z1D9*FAL(3,21)-Z2D9*RSU)
     *            - F(20)*(Z2D9*FAL(3,20)-Z1D9*FAL(3,11)-Z2D9*RSU)
     *            - F(24)*(Z2D9*FAL(3,24)-Z1D9*FAL(3,15))
        END IF
      END IF
      R E T U R N
      END
