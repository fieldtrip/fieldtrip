      SUBROUTINE FOFU23(IBEWER,IBEABL,R,S,T,LKNEL ,FWE   ,FAL   )       NEU
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
C FOFU23:
C =======
C
C BERECHNUNG DER WERTE UND DER LOKALEN ABLEITUNGEN DER FORMFUNKTIONEN
C DES ISOPARAMETRISCHEN QUADERS (ELEMENTTYP XXX23)
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
      DO 10 IKNEL=9,20
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMQ=ISUMQ+1
      ELSE
        F(IKNEL)=Z0
      END IF
   10 CONTINUE
      DO 20 IKNEL=21,32
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMK=ISUMK+1
      ELSE
        F(IKNEL)=Z0
      END IF
   20 CONTINUE
      RMI=Z1-R
      SMI=Z1-S
      TMI=Z1-T
      RPL=Z1+R
      SPL=Z1+S
      TPL=Z1+T
      RMQ=Z1-R*R
      SMQ=Z1-S*S
      TMQ=Z1-T*T
      RST=R*R+S*S+T*T-Z19D9
      R3P=Z1+Z3*R
      S3P=Z1+Z3*S
      T3P=Z1+Z3*T
      R3M=Z1-Z3*R
      S3M=Z1-Z3*S
      T3M=Z1-Z3*T
      IF(ISUMQ.EQ.12.AND.ISUMK.EQ.12)THEN
C +++++
C VOLL-KUBISCHER QUADER
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z9D64*RMI*SMI*TMI*RST
          FWE( 2)=Z9D64*RPL*SMI*TMI*RST
          FWE( 3)=Z9D64*RPL*SPL*TMI*RST
          FWE( 4)=Z9D64*RMI*SPL*TMI*RST
          FWE( 5)=Z9D64*RMI*SMI*TPL*RST
          FWE( 6)=Z9D64*RPL*SMI*TPL*RST
          FWE( 7)=Z9D64*RPL*SPL*TPL*RST
          FWE( 8)=Z9D64*RMI*SPL*TPL*RST
          FWE( 9)=Z9D64*RMQ*R3M*SMI*TMI
          FWE(10)=Z9D64*SMQ*S3M*RPL*TMI
          FWE(11)=Z9D64*RMQ*R3P*SPL*TMI
          FWE(12)=Z9D64*SMQ*S3P*RMI*TMI
          FWE(13)=Z9D64*RMQ*R3M*SMI*TPL
          FWE(14)=Z9D64*SMQ*S3M*RPL*TPL
          FWE(15)=Z9D64*RMQ*R3P*SPL*TPL
          FWE(16)=Z9D64*SMQ*S3P*RMI*TPL
          FWE(17)=Z9D64*TMQ*T3M*RMI*SMI
          FWE(18)=Z9D64*TMQ*T3M*RPL*SMI
          FWE(19)=Z9D64*TMQ*T3M*RPL*SPL
          FWE(20)=Z9D64*TMQ*T3M*RMI*SPL
          FWE(21)=Z9D64*RMQ*R3P*SMI*TMI
          FWE(22)=Z9D64*SMQ*S3P*RPL*TMI
          FWE(23)=Z9D64*RMQ*R3M*SPL*TMI
          FWE(24)=Z9D64*SMQ*S3M*RMI*TMI
          FWE(25)=Z9D64*RMQ*R3P*SMI*TPL
          FWE(26)=Z9D64*SMQ*S3P*RPL*TPL
          FWE(27)=Z9D64*RMQ*R3M*SPL*TPL
          FWE(28)=Z9D64*SMQ*S3M*RMI*TPL
          FWE(29)=Z9D64*TMQ*T3P*RMI*SMI
          FWE(30)=Z9D64*TMQ*T3P*RPL*SMI
          FWE(31)=Z9D64*TMQ*T3P*RPL*SPL
          FWE(32)=Z9D64*TMQ*T3P*RMI*SPL
        END IF
        IF(IBEABL.EQ.1)THEN
          CONTINUE
          FAL(1, 1)= Z9D64*SMI*TMI*(Z2*R*RMI-RST)
          FAL(2, 1)= Z9D64*RMI*TMI*(Z2*S*SMI-RST)
          FAL(3, 1)= Z9D64*RMI*SMI*(Z2*T*TMI-RST)
          FAL(1, 2)= Z9D64*SMI*TMI*(Z2*R*RPL+RST)
          FAL(2, 2)= Z9D64*RPL*TMI*(Z2*S*SMI-RST)
          FAL(3, 2)= Z9D64*RPL*SMI*(Z2*T*TMI-RST)
          FAL(1, 3)= Z9D64*SPL*TMI*(Z2*R*RPL+RST)
          FAL(2, 3)= Z9D64*RPL*TMI*(Z2*S*SPL+RST)
          FAL(3, 3)= Z9D64*RPL*SPL*(Z2*T*TMI-RST)
          FAL(1, 4)= Z9D64*SPL*TMI*(Z2*R*RMI-RST)
          FAL(2, 4)= Z9D64*RMI*TMI*(Z2*S*SPL+RST)
          FAL(3, 4)= Z9D64*RMI*SPL*(Z2*T*TMI-RST)
          FAL(1, 5)= Z9D64*SMI*TPL*(Z2*R*RMI-RST)
          FAL(2, 5)= Z9D64*RMI*TPL*(Z2*S*SMI-RST)
          FAL(3, 5)= Z9D64*RMI*SMI*(Z2*T*TPL+RST)
          FAL(1, 6)= Z9D64*SMI*TPL*(Z2*R*RPL+RST)
          FAL(2, 6)= Z9D64*RPL*TPL*(Z2*S*SMI-RST)
          FAL(3, 6)= Z9D64*RPL*SMI*(Z2*T*TPL+RST)
          FAL(1, 7)= Z9D64*SPL*TPL*(Z2*R*RPL+RST)
          FAL(2, 7)= Z9D64*RPL*TPL*(Z2*S*SPL+RST)
          FAL(3, 7)= Z9D64*RPL*SPL*(Z2*T*TPL+RST)
          FAL(1, 8)= Z9D64*SPL*TPL*(Z2*R*RMI-RST)
          FAL(2, 8)= Z9D64*RMI*TPL*(Z2*S*SPL+RST)
          FAL(3, 8)= Z9D64*RMI*SPL*(Z2*T*TPL+RST)
          FAL(1, 9)=-Z9D64*SMI*TMI*(R3P*R3M+Z2*RPL)
          FAL(2, 9)=-Z9D64*RMQ*R3M*TMI
          FAL(3, 9)=-Z9D64*RMQ*R3M*SMI
          FAL(1,10)= Z9D64*SMQ*S3M*TMI
          FAL(2,10)=-Z9D64*RPL*TMI*(S3P*S3M+Z2*SPL)
          FAL(3,10)=-Z9D64*SMQ*S3M*RPL
          FAL(1,11)= Z9D64*SPL*TMI*(R3P*R3M+Z2*RMI)
          FAL(2,11)= Z9D64*RMQ*R3P*TMI
          FAL(3,11)=-Z9D64*RMQ*R3P*SPL
          FAL(1,12)=-Z9D64*SMQ*S3P*TMI
          FAL(2,12)= Z9D64*RMI*TMI*(S3P*S3M+Z2*SMI)
          FAL(3,12)=-Z9D64*SMQ*S3P*RMI
          FAL(1,13)=-Z9D64*SMI*TPL*(R3P*R3M+Z2*RPL)
          FAL(2,13)=-Z9D64*RMQ*R3M*TPL
          FAL(3,13)= Z9D64*RMQ*R3M*SMI
          FAL(1,14)= Z9D64*SMQ*S3M*TPL
          FAL(2,14)=-Z9D64*RPL*TPL*(S3P*S3M+Z2*SPL)
          FAL(3,14)= Z9D64*SMQ*S3M*RPL
          FAL(1,15)= Z9D64*SPL*TPL*(R3P*R3M+Z2*RMI)
          FAL(2,15)= Z9D64*RMQ*R3P*TPL
          FAL(3,15)= Z9D64*RMQ*R3P*SPL
          FAL(1,16)=-Z9D64*SMQ*S3P*TPL
          FAL(2,16)= Z9D64*RMI*TPL*(S3P*S3M+Z2*SMI)
          FAL(3,16)= Z9D64*SMQ*S3P*RMI
          FAL(1,17)=-Z9D64*TMQ*T3M*SMI
          FAL(2,17)=-Z9D64*TMQ*T3M*RMI
          FAL(3,17)=-Z9D64*RMI*SMI*(T3P*T3M+Z2*TPL)
          FAL(1,18)= Z9D64*TMQ*T3M*SMI
          FAL(2,18)=-Z9D64*TMQ*T3M*RPL
          FAL(3,18)=-Z9D64*RPL*SMI*(T3P*T3M+Z2*TPL)
          FAL(1,19)= Z9D64*TMQ*T3M*SPL
          FAL(2,19)= Z9D64*TMQ*T3M*RPL
          FAL(3,19)=-Z9D64*RPL*SPL*(T3P*T3M+Z2*TPL)
          FAL(1,20)=-Z9D64*TMQ*T3M*SPL
          FAL(2,20)= Z9D64*TMQ*T3M*RMI
          FAL(3,20)=-Z9D64*RMI*SPL*(T3P*T3M+Z2*TPL)
          FAL(1,21)= Z9D64*SMI*TMI*(R3P*R3M+Z2*RMI)
          FAL(2,21)=-Z9D64*RMQ*R3P*TMI
          FAL(3,21)=-Z9D64*RMQ*R3P*SMI
          FAL(1,22)= Z9D64*SMQ*S3P*TMI
          FAL(2,22)= Z9D64*RPL*TMI*(S3P*S3M+Z2*SMI)
          FAL(3,22)=-Z9D64*SMQ*S3P*RPL
          FAL(1,23)=-Z9D64*SPL*TMI*(R3P*R3M+Z2*RPL)
          FAL(2,23)= Z9D64*RMQ*R3M*TMI
          FAL(3,23)=-Z9D64*RMQ*R3M*SPL
          FAL(1,24)=-Z9D64*SMQ*S3M*TMI
          FAL(2,24)=-Z9D64*RMI*TMI*(S3P*S3M+Z2*SPL)
          FAL(3,24)=-Z9D64*SMQ*S3M*RMI
          FAL(1,25)= Z9D64*SMI*TPL*(R3P*R3M+Z2*RMI)
          FAL(2,25)=-Z9D64*RMQ*R3P*TPL
          FAL(3,25)= Z9D64*RMQ*R3P*SMI
          FAL(1,26)= Z9D64*SMQ*S3P*TPL
          FAL(2,26)= Z9D64*RPL*TPL*(S3P*S3M+Z2*SMI)
          FAL(3,26)= Z9D64*SMQ*S3P*RPL
          FAL(1,27)=-Z9D64*SPL*TPL*(R3P*R3M+Z2*RPL)
          FAL(2,27)= Z9D64*RMQ*R3M*TPL
          FAL(3,27)= Z9D64*RMQ*R3M*SPL
          FAL(1,28)=-Z9D64*SMQ*S3M*TPL
          FAL(2,28)=-Z9D64*RMI*TPL*(S3P*S3M+Z2*SPL)
          FAL(3,28)= Z9D64*SMQ*S3M*RMI
          FAL(1,29)=-Z9D64*TMQ*T3P*SMI
          FAL(2,29)=-Z9D64*TMQ*T3P*RMI
          FAL(3,29)= Z9D64*RMI*SMI*(T3P*T3M+Z2*TMI)
          FAL(1,30)= Z9D64*TMQ*T3P*SMI
          FAL(2,30)=-Z9D64*TMQ*T3P*RPL
          FAL(3,30)= Z9D64*RPL*SMI*(T3P*T3M+Z2*TMI)
          FAL(1,31)= Z9D64*TMQ*T3P*SPL
          FAL(2,31)= Z9D64*TMQ*T3P*RPL
          FAL(3,31)= Z9D64*RPL*SPL*(T3P*T3M+Z2*TMI)
          FAL(1,32)=-Z9D64*TMQ*T3P*SPL
          FAL(2,32)= Z9D64*TMQ*T3P*RMI
          FAL(3,32)= Z9D64*RMI*SPL*(T3P*T3M+Z2*TMI)
        END IF
      ELSE IF(ISUMQ.EQ.12.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-QUADRATISCHER QUADER
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z1D8*RMI*SMI*TMI*(-R-S-T-Z2)
          FWE( 2)=Z1D8*RPL*SMI*TMI*( R-S-T-Z2)
          FWE( 3)=Z1D8*RPL*SPL*TMI*( R+S-T-Z2)
          FWE( 4)=Z1D8*RMI*SPL*TMI*(-R+S-T-Z2)
          FWE( 5)=Z1D8*RMI*SMI*TPL*(-R-S+T-Z2)
          FWE( 6)=Z1D8*RPL*SMI*TPL*( R-S+T-Z2)
          FWE( 7)=Z1D8*RPL*SPL*TPL*( R+S+T-Z2)
          FWE( 8)=Z1D8*RMI*SPL*TPL*(-R+S+T-Z2)
          FWE( 9)=Z025*RMQ*SMI*TMI
          FWE(10)=Z025*RPL*SMQ*TMI
          FWE(11)=Z025*RMQ*SPL*TMI
          FWE(12)=Z025*RMI*SMQ*TMI
          FWE(13)=Z025*RMQ*SMI*TPL
          FWE(14)=Z025*RPL*SMQ*TPL
          FWE(15)=Z025*RMQ*SPL*TPL
          FWE(16)=Z025*RMI*SMQ*TPL
          FWE(17)=Z025*RMI*SMI*TMQ
          FWE(18)=Z025*RPL*SMI*TMQ
          FWE(19)=Z025*RPL*SPL*TMQ
          FWE(20)=Z025*RMI*SPL*TMQ
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1, 1)= Z1D8     *SMI*TMI*( Z2*R+   S+   T+Z1)
          FAL(2, 1)= Z1D8*RMI     *TMI*(    R+Z2*S+   T+Z1)
          FAL(3, 1)= Z1D8*RMI*SMI     *(    R+   S+Z2*T+Z1)
          FAL(1, 2)= Z1D8     *SMI*TMI*( Z2*R-   S-   T-Z1)
          FAL(2, 2)= Z1D8*RPL     *TMI*(-   R+Z2*S+   T+Z1)
          FAL(3, 2)= Z1D8*RPL*SMI     *(-   R+   S+Z2*T+Z1)
          FAL(1, 3)= Z1D8     *SPL*TMI*( Z2*R+   S-   T-Z1)
          FAL(2, 3)= Z1D8*RPL     *TMI*(    R+Z2*S-   T-Z1)
          FAL(3, 3)= Z1D8*RPL*SPL     *(-   R-   S+Z2*T+Z1)
          FAL(1, 4)= Z1D8     *SPL*TMI*( Z2*R-   S+   T+Z1)
          FAL(2, 4)= Z1D8*RMI     *TMI*(-   R+Z2*S-   T-Z1)
          FAL(3, 4)= Z1D8*RMI*SPL     *(    R-   S+Z2*T+Z1)
          FAL(1, 5)= Z1D8     *SMI*TPL*( Z2*R+   S-   T+Z1)
          FAL(2, 5)= Z1D8*RMI     *TPL*(    R+Z2*S-   T+Z1)
          FAL(3, 5)= Z1D8*RMI*SMI     *(-   R-   S+Z2*T-Z1)
          FAL(1, 6)= Z1D8     *SMI*TPL*( Z2*R-   S+   T-Z1)
          FAL(2, 6)= Z1D8*RPL     *TPL*(-   R+Z2*S-   T+Z1)
          FAL(3, 6)= Z1D8*RPL*SMI     *(    R-   S+Z2*T-Z1)
          FAL(1, 7)= Z1D8     *SPL*TPL*( Z2*R+   S+   T-Z1)
          FAL(2, 7)= Z1D8*RPL     *TPL*(    R+Z2*S+   T-Z1)
          FAL(3, 7)= Z1D8*RPL*SPL     *(    R+   S+Z2*T-Z1)
          FAL(1, 8)= Z1D8     *SPL*TPL*( Z2*R-   S-   T+Z1)
          FAL(2, 8)= Z1D8*RMI     *TPL*(-   R+Z2*S+   T-Z1)
          FAL(3, 8)= Z1D8*RMI*SPL     *(-   R+   S+Z2*T-Z1)
          FAL(1, 9)=-Z05 *R  *SMI*TMI
          FAL(2, 9)=-Z025*RMQ    *TMI
          FAL(3, 9)=-Z025*RMQ*SMI
          FAL(1,10)= Z025    *SMQ*TMI
          FAL(2,10)=-Z05 *RPL*S  *TMI
          FAL(3,10)=-Z025*RPL*SMQ
          FAL(1,11)=-Z05 *R  *SPL*TMI
          FAL(2,11)= Z025*RMQ    *TMI
          FAL(3,11)=-Z025*RMQ*SPL
          FAL(1,12)=-Z025    *SMQ*TMI
          FAL(2,12)=-Z05 *RMI*S  *TMI
          FAL(3,12)=-Z025*RMI*SMQ
          FAL(1,13)=-Z05 *R  *SMI*TPL
          FAL(2,13)=-Z025*RMQ    *TPL
          FAL(3,13)= Z025*RMQ*SMI
          FAL(1,14)= Z025    *SMQ*TPL
          FAL(2,14)=-Z05 *RPL*S  *TPL
          FAL(3,14)= Z025*RPL*SMQ
          FAL(1,15)=-Z05 *R  *SPL*TPL
          FAL(2,15)= Z025*RMQ    *TPL
          FAL(3,15)= Z025*RMQ*SPL
          FAL(1,16)=-Z025    *SMQ*TPL
          FAL(2,16)=-Z05 *RMI*S  *TPL
          FAL(3,16)= Z025*RMI*SMQ
          FAL(1,17)=-Z025    *SMI*TMQ
          FAL(2,17)=-Z025*RMI    *TMQ
          FAL(3,17)=-Z05 *RMI*SMI*T
          FAL(1,18)= Z025    *SMI*TMQ
          FAL(2,18)=-Z025*RPL    *TMQ
          FAL(3,18)=-Z05 *RPL*SMI*T
          FAL(1,19)= Z025    *SPL*TMQ
          FAL(2,19)= Z025*RPL    *TMQ
          FAL(3,19)=-Z05 *RPL*SPL*T
          FAL(1,20)=-Z025    *SPL*TMQ
          FAL(2,20)= Z025*RMI    *TMQ
          FAL(3,20)=-Z05 *RMI*SPL*T
        END IF
      ELSE IF(ISUMQ.EQ.0.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-LINEARER QUADER
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z1D8*RMI*SMI*TMI
          FWE( 2)=Z1D8*RPL*SMI*TMI
          FWE( 3)=Z1D8*RPL*SPL*TMI
          FWE( 4)=Z1D8*RMI*SPL*TMI
          FWE( 5)=Z1D8*RMI*SMI*TPL
          FWE( 6)=Z1D8*RPL*SMI*TPL
          FWE( 7)=Z1D8*RPL*SPL*TPL
          FWE( 8)=Z1D8*RMI*SPL*TPL
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,1)=-Z1D8*     SMI*TMI
          FAL(2,1)=-Z1D8*RMI*     TMI
          FAL(3,1)=-Z1D8*RMI*SMI
          FAL(1,2)= Z1D8*     SMI*TMI
          FAL(2,2)=-Z1D8*RPL*     TMI
          FAL(3,2)=-Z1D8*RPL*SMI
          FAL(1,3)= Z1D8*     SPL*TMI
          FAL(2,3)= Z1D8*RPL*     TMI
          FAL(3,3)=-Z1D8*RPL*SPL
          FAL(1,4)=-Z1D8*     SPL*TMI
          FAL(2,4)= Z1D8*RMI*     TMI
          FAL(3,4)=-Z1D8*RMI*SPL
          FAL(1,5)=-Z1D8*     SMI*TPL
          FAL(2,5)=-Z1D8*RMI*     TPL
          FAL(3,5)= Z1D8*RMI*SMI
          FAL(1,6)= Z1D8*     SMI*TPL
          FAL(2,6)=-Z1D8*RPL*     TPL
          FAL(3,6)= Z1D8*RPL*SMI
          FAL(1,7)= Z1D8*     SPL*TPL
          FAL(2,7)= Z1D8*RPL*     TPL
          FAL(3,7)= Z1D8*RPL*SPL
          FAL(1,8)=-Z1D8*     SPL*TPL
          FAL(2,8)= Z1D8*RMI*     TPL
          FAL(3,8)= Z1D8*RMI*SPL
        END IF
      ELSE
C +++++
C QUADER MIT VARIABLER KNOTENZAHL
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(21)=F(21)*Z9D64*RMQ*R3P*SMI*TMI
          FWE(22)=F(22)*Z9D64*SMQ*S3P*RPL*TMI
          FWE(23)=F(23)*Z9D64*RMQ*R3M*SPL*TMI
          FWE(24)=F(24)*Z9D64*SMQ*S3M*RMI*TMI
          FWE(25)=F(25)*Z9D64*RMQ*R3P*SMI*TPL
          FWE(26)=F(26)*Z9D64*SMQ*S3P*RPL*TPL
          FWE(27)=F(27)*Z9D64*RMQ*R3M*SPL*TPL
          FWE(28)=F(28)*Z9D64*SMQ*S3M*RMI*TPL
          FWE(29)=F(29)*Z9D64*TMQ*T3P*RMI*SMI
          FWE(30)=F(30)*Z9D64*TMQ*T3P*RPL*SMI
          FWE(31)=F(31)*Z9D64*TMQ*T3P*RPL*SPL
          FWE(32)=F(32)*Z9D64*TMQ*T3P*RMI*SPL
          FWE( 9)=F( 9)*(Z025*SMI*TMI*(RMQ+F(21)*Z1D8*RMQ)-FWE(21))
          FWE(10)=F(10)*(Z025*RPL*TMI*(SMQ+F(22)*Z1D8*SMQ)-FWE(22))
          FWE(11)=F(11)*(Z025*SPL*TMI*(RMQ+F(23)*Z1D8*RMQ)-FWE(23))
          FWE(12)=F(12)*(Z025*RMI*TMI*(SMQ+F(24)*Z1D8*SMQ)-FWE(24))
          FWE(13)=F(13)*(Z025*SMI*TPL*(RMQ+F(25)*Z1D8*RMQ)-FWE(25))
          FWE(14)=F(14)*(Z025*RPL*TPL*(SMQ+F(26)*Z1D8*SMQ)-FWE(26))
          FWE(15)=F(15)*(Z025*SPL*TPL*(RMQ+F(27)*Z1D8*RMQ)-FWE(27))
          FWE(16)=F(16)*(Z025*RMI*TPL*(SMQ+F(28)*Z1D8*SMQ)-FWE(28))
          FWE(17)=F(17)*(Z025*RMI*SMI*(TMQ+F(29)*Z1D8*TMQ)-FWE(29))
          FWE(18)=F(18)*(Z025*RPL*SMI*(TMQ+F(30)*Z1D8*TMQ)-FWE(30))
          FWE(19)=F(19)*(Z025*RPL*SPL*(TMQ+F(31)*Z1D8*TMQ)-FWE(31))
          FWE(20)=F(20)*(Z025*RMI*SPL*(TMQ+F(32)*Z1D8*TMQ)-FWE(32))
          FWE( 1)=Z1D8*(RMI*SMI*TMI
     *           -F( 9)*RMQ*SMI*TMI-F(12)*SMQ*RMI*TMI-F(17)*TMQ*RMI*SMI)
     *           -F(21)*(Z2D9*FWE( 9)-Z1D9*FWE(21))
     *           -F(24)*(Z2D9*FWE(24)-Z1D9*FWE(12))
     *           -F(29)*(Z2D9*FWE(17)-Z1D9*FWE(29))
          FWE( 2)=Z1D8*(RPL*SMI*TMI
     *           -F( 9)*RMQ*SMI*TMI-F(10)*SMQ*RPL*TMI-F(18)*TMQ*RPL*SMI)
     *           -F(21)*(Z2D9*FWE(21)-Z1D9*FWE( 9))
     *           -F(22)*(Z2D9*FWE(10)-Z1D9*FWE(22))
     *           -F(30)*(Z2D9*FWE(18)-Z1D9*FWE(30))
          FWE( 3)=Z1D8*(RPL*SPL*TMI
     *           -F(10)*SMQ*RPL*TMI-F(11)*RMQ*SPL*TMI-F(19)*TMQ*RPL*SPL)
     *           -F(22)*(Z2D9*FWE(22)-Z1D9*FWE(10))
     *           -F(23)*(Z2D9*FWE(11)-Z1D9*FWE(23))
     *           -F(31)*(Z2D9*FWE(19)-Z1D9*FWE(31))
          FWE( 4)=Z1D8*(RMI*SPL*TMI
     *           -F(11)*RMQ*SPL*TMI-F(12)*SMQ*RMI*TMI-F(20)*TMQ*RMI*SPL)
     *           -F(23)*(Z2D9*FWE(23)-Z1D9*FWE(11))
     *           -F(24)*(Z2D9*FWE(12)-Z1D9*FWE(24))
     *           -F(32)*(Z2D9*FWE(20)-Z1D9*FWE(32))
          FWE( 5)=Z1D8*(RMI*SMI*TPL
     *           -F(13)*RMQ*SMI*TPL-F(16)*SMQ*RMI*TPL-F(17)*TMQ*RMI*SMI)
     *           -F(25)*(Z2D9*FWE(13)-Z1D9*FWE(25))
     *           -F(28)*(Z2D9*FWE(28)-Z1D9*FWE(16))
     *           -F(29)*(Z2D9*FWE(29)-Z1D9*FWE(17))
          FWE( 6)=Z1D8*(RPL*SMI*TPL
     *           -F(13)*RMQ*SMI*TPL-F(14)*SMQ*RPL*TPL-F(18)*TMQ*RPL*SMI)
     *           -F(25)*(Z2D9*FWE(25)-Z1D9*FWE(13))
     *           -F(26)*(Z2D9*FWE(14)-Z1D9*FWE(26))
     *           -F(30)*(Z2D9*FWE(30)-Z1D9*FWE(18))
          FWE( 7)=Z1D8*(RPL*SPL*TPL
     *           -F(14)*SMQ*RPL*TPL-F(15)*RMQ*SPL*TPL-F(19)*TMQ*RPL*SPL)
     *           -F(26)*(Z2D9*FWE(26)-Z1D9*FWE(14))
     *           -F(27)*(Z2D9*FWE(15)-Z1D9*FWE(27))
     *           -F(31)*(Z2D9*FWE(31)-Z1D9*FWE(19))
          FWE( 8)=Z1D8*(RMI*SPL*TPL
     *           -F(15)*RMQ*SPL*TPL-F(16)*SMQ*RMI*TPL-F(20)*TMQ*RMI*SPL)
     *           -F(27)*(Z2D9*FWE(27)-Z1D9*FWE(15))
     *           -F(28)*(Z2D9*FWE(16)-Z1D9*FWE(28))
     *           -F(32)*(Z2D9*FWE(32)-Z1D9*FWE(20))
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,21)= F(21)*Z9D64*SMI*TMI*(R3P*R3M+Z2*RMI)
          FAL(2,21)=-F(21)*Z9D64*RMQ*R3P*TMI
          FAL(3,21)=-F(21)*Z9D64*RMQ*R3P*SMI
          FAL(1,22)= F(22)*Z9D64*SMQ*S3P*TMI
          FAL(2,22)= F(22)*Z9D64*RPL*TMI*(S3P*S3M+Z2*SMI)
          FAL(3,22)=-F(22)*Z9D64*SMQ*S3P*RPL
          FAL(1,23)=-F(23)*Z9D64*SPL*TMI*(R3P*R3M+Z2*RPL)
          FAL(2,23)= F(23)*Z9D64*RMQ*R3M*TMI
          FAL(3,23)=-F(23)*Z9D64*RMQ*R3M*SPL
          FAL(1,24)=-F(24)*Z9D64*SMQ*S3M*TMI
          FAL(2,24)=-F(24)*Z9D64*RMI*TMI*(S3P*S3M+Z2*SPL)
          FAL(3,24)=-F(24)*Z9D64*SMQ*S3M*RMI
          FAL(1,25)= F(25)*Z9D64*SMI*TPL*(R3P*R3M+Z2*RMI)
          FAL(2,25)=-F(25)*Z9D64*RMQ*R3P*TPL
          FAL(3,25)= F(25)*Z9D64*RMQ*R3P*SMI
          FAL(1,26)= F(26)*Z9D64*SMQ*S3P*TPL
          FAL(2,26)= F(26)*Z9D64*RPL*TPL*(S3P*S3M+Z2*SMI)
          FAL(3,26)= F(26)*Z9D64*SMQ*S3P*RPL
          FAL(1,27)=-F(27)*Z9D64*SPL*TPL*(R3P*R3M+Z2*RPL)
          FAL(2,27)= F(27)*Z9D64*RMQ*R3M*TPL
          FAL(3,27)= F(27)*Z9D64*RMQ*R3M*SPL
          FAL(1,28)=-F(28)*Z9D64*SMQ*S3M*TPL
          FAL(2,28)=-F(28)*Z9D64*RMI*TPL*(S3P*S3M+Z2*SPL)
          FAL(3,28)= F(28)*Z9D64*SMQ*S3M*RMI
          FAL(1,29)=-F(29)*Z9D64*TMQ*T3P*SMI
          FAL(2,29)=-F(29)*Z9D64*TMQ*T3P*RMI
          FAL(3,29)= F(29)*Z9D64*RMI*SMI*(T3P*T3M+Z2*TMI)
          FAL(1,30)= F(30)*Z9D64*TMQ*T3P*SMI
          FAL(2,30)=-F(30)*Z9D64*TMQ*T3P*RPL
          FAL(3,30)= F(30)*Z9D64*RPL*SMI*(T3P*T3M+Z2*TMI)
          FAL(1,31)= F(31)*Z9D64*TMQ*T3P*SPL
          FAL(2,31)= F(31)*Z9D64*TMQ*T3P*RPL
          FAL(3,31)= F(31)*Z9D64*RPL*SPL*(T3P*T3M+Z2*TMI)
          FAL(1,32)=-F(32)*Z9D64*TMQ*T3P*SPL
          FAL(2,32)= F(32)*Z9D64*TMQ*T3P*RMI
          FAL(3,32)= F(32)*Z9D64*RMI*SPL*(T3P*T3M+Z2*TMI)
          FAL(1, 9)=-F( 9)*(Z05 *SMI*TMI*(R+F(21)*Z1D8*R)+FAL(1,21))
          FAL(2, 9)=-F( 9)*(Z025*TMI*(RMQ+F(21)*Z1D8*RMQ)+FAL(2,21))
          FAL(3, 9)=-F( 9)*(Z025*SMI*(RMQ+F(21)*Z1D8*RMQ)+FAL(3,21))
          FAL(1,10)= F(10)*(Z025*TMI*(SMQ+F(22)*Z1D8*SMQ)-FAL(1,22))
          FAL(2,10)=-F(10)*(Z05 *RPL*TMI*(S+F(22)*Z1D8*S)+FAL(2,22))
          FAL(3,10)=-F(10)*(Z025*RPL*(SMQ+F(22)*Z1D8*SMQ)+FAL(3,22))
          FAL(1,11)=-F(11)*(Z05 *SPL*TMI*(R+F(23)*Z1D8*R)+FAL(1,23))
          FAL(2,11)= F(11)*(Z025*TMI*(RMQ+F(23)*Z1D8*RMQ)-FAL(2,23))
          FAL(3,11)=-F(11)*(Z025*SPL*(RMQ+F(23)*Z1D8*RMQ)+FAL(3,23))
          FAL(1,12)=-F(12)*(Z025*TMI*(SMQ+F(24)*Z1D8*SMQ)+FAL(1,24))
          FAL(2,12)=-F(12)*(Z05 *RMI*TMI*(S+F(24)*Z1D8*S)+FAL(2,24))
          FAL(3,12)=-F(12)*(Z025*RMI*(SMQ+F(24)*Z1D8*SMQ)+FAL(3,24))
          FAL(1,13)=-F(13)*(Z05 *SMI*TPL*(R+F(25)*Z1D8*R)+FAL(1,25))
          FAL(2,13)=-F(13)*(Z025*TPL*(RMQ+F(25)*Z1D8*RMQ)+FAL(2,25))
          FAL(3,13)= F(13)*(Z025*SMI*(RMQ+F(25)*Z1D8*RMQ)-FAL(3,25))
          FAL(1,14)= F(14)*(Z025*TPL*(SMQ+F(26)*Z1D8*SMQ)-FAL(1,26))
          FAL(2,14)=-F(14)*(Z05 *RPL*TPL*(S+F(26)*Z1D8*S)+FAL(2,26))
          FAL(3,14)= F(14)*(Z025*RPL*(SMQ+F(26)*Z1D8*SMQ)-FAL(3,26))
          FAL(1,15)=-F(15)*(Z05 *SPL*TPL*(R+F(27)*Z1D8*R)+FAL(1,27))
          FAL(2,15)= F(15)*(Z025*TPL*(RMQ+F(27)*Z1D8*RMQ)-FAL(2,27))
          FAL(3,15)= F(15)*(Z025*SPL*(RMQ+F(27)*Z1D8*RMQ)-FAL(3,27))
          FAL(1,16)=-F(16)*(Z025*TPL*(SMQ+F(28)*Z1D8*SMQ)+FAL(1,28))
          FAL(2,16)=-F(16)*(Z05 *RMI*TPL*(S+F(28)*Z1D8*S)+FAL(2,28))
          FAL(3,16)= F(16)*(Z025*RMI*(SMQ+F(28)*Z1D8*SMQ)-FAL(3,28))
          FAL(1,17)=-F(17)*(Z025*SMI*(TMQ+F(29)*Z1D8*TMQ)+FAL(1,29))
          FAL(2,17)=-F(17)*(Z025*RMI*(TMQ+F(29)*Z1D8*TMQ)+FAL(2,29))
          FAL(3,17)=-F(17)*(Z05 *RMI*SMI*(T+F(29)*Z1D8*T)+FAL(3,29))
          FAL(1,18)= F(18)*(Z025*SMI*(TMQ+F(30)*Z1D8*TMQ)-FAL(1,30))
          FAL(2,18)=-F(18)*(Z025*RPL*(TMQ+F(30)*Z1D8*TMQ)+FAL(2,30))
          FAL(3,18)=-F(18)*(Z05 *RPL*SMI*(T+F(30)*Z1D8*T)+FAL(3,30))
          FAL(1,19)= F(19)*(Z025*SPL*(TMQ+F(31)*Z1D8*TMQ)-FAL(1,31))
          FAL(2,19)= F(19)*(Z025*RPL*(TMQ+F(31)*Z1D8*TMQ)-FAL(2,31))
          FAL(3,19)=-F(19)*(Z05 *RPL*SPL*(T+F(31)*Z1D8*T)+FAL(3,31))
          FAL(1,20)=-F(20)*(Z025*SPL*(TMQ+F(32)*Z1D8*TMQ)+FAL(1,32))
          FAL(2,20)= F(20)*(Z025*RMI*(TMQ+F(32)*Z1D8*TMQ)-FAL(2,32))
          FAL(3,20)=-F(20)*(Z05 *RMI*SPL*(T+F(32)*Z1D8*T)+FAL(3,32))
          FAL(1, 1)=Z1D8*(-SMI*TMI+F( 9)*Z2*R*SMI*TMI
     *                            +F(12)*SMQ*TMI
     *                            +F(17)*TMQ*SMI)
     *             -F(21)*(Z2D9*FAL(1, 9)-Z1D9*FAL(1,21))
     *             -F(24)*(Z2D9*FAL(1,24)-Z1D9*FAL(1,12))
     *             -F(29)*(Z2D9*FAL(1,17)-Z1D9*FAL(1,29))
          FAL(2, 1)=Z1D8*(-RMI*TMI+F( 9)*RMQ*TMI
     *                            +F(12)*Z2*S*RMI*TMI
     *                            +F(17)*TMQ*RMI)
     *             -F(21)*(Z2D9*FAL(2, 9)-Z1D9*FAL(2,21))
     *             -F(24)*(Z2D9*FAL(2,24)-Z1D9*FAL(2,12))
     *             -F(29)*(Z2D9*FAL(2,17)-Z1D9*FAL(2,29))
          FAL(3, 1)=Z1D8*(-RMI*SMI+F( 9)*RMQ*SMI
     *                            +F(12)*SMQ*RMI
     *                            +F(17)*Z2*T*RMI*SMI)
     *             -F(21)*(Z2D9*FAL(3, 9)-Z1D9*FAL(3,21))
     *             -F(24)*(Z2D9*FAL(3,24)-Z1D9*FAL(3,12))
     *             -F(29)*(Z2D9*FAL(3,17)-Z1D9*FAL(3,29))
          FAL(1, 2)=Z1D8*( SMI*TMI+F( 9)*Z2*R*SMI*TMI
     *                            -F(10)*SMQ*TMI
     *                            -F(18)*TMQ*SMI)
     *             -F(21)*(Z2D9*FAL(1,21)-Z1D9*FAL(1, 9))
     *             -F(22)*(Z2D9*FAL(1,10)-Z1D9*FAL(1,22))
     *             -F(30)*(Z2D9*FAL(1,18)-Z1D9*FAL(1,30))
          FAL(2, 2)=Z1D8*(-RPL*TMI+F( 9)*RMQ*TMI
     *                            +F(10)*Z2*S*RPL*TMI
     *                            +F(18)*TMQ*RPL)
     *             -F(21)*(Z2D9*FAL(2,21)-Z1D9*FAL(2, 9))
     *             -F(22)*(Z2D9*FAL(2,10)-Z1D9*FAL(2,22))
     *             -F(30)*(Z2D9*FAL(2,18)-Z1D9*FAL(2,30))
          FAL(3, 2)=Z1D8*(-RPL*SMI+F( 9)*RMQ*SMI
     *                            +F(10)*SMQ*RPL
     *                            +F(18)*Z2*T*RPL*SMI)
     *             -F(21)*(Z2D9*FAL(3,21)-Z1D9*FAL(3, 9))
     *             -F(22)*(Z2D9*FAL(3,10)-Z1D9*FAL(3,22))
     *             -F(30)*(Z2D9*FAL(3,18)-Z1D9*FAL(3,30))
          FAL(1, 3)=Z1D8*( SPL*TMI-F(10)*SMQ*TMI
     *                            +F(11)*Z2*R*SPL*TMI
     *                            -F(19)*TMQ*SPL)
     *             -F(22)*(Z2D9*FAL(1,22)-Z1D9*FAL(1,10))
     *             -F(23)*(Z2D9*FAL(1,11)-Z1D9*FAL(1,23))
     *             -F(31)*(Z2D9*FAL(1,19)-Z1D9*FAL(1,31))
          FAL(2, 3)=Z1D8*( RPL*TMI+F(10)*Z2*S*RPL*TMI
     *                            -F(11)*RMQ*TMI
     *                            -F(19)*TMQ*RPL)
     *             -F(22)*(Z2D9*FAL(2,22)-Z1D9*FAL(2,10))
     *             -F(23)*(Z2D9*FAL(2,11)-Z1D9*FAL(2,23))
     *             -F(31)*(Z2D9*FAL(2,19)-Z1D9*FAL(2,31))
          FAL(3, 3)=Z1D8*(-RPL*SPL+F(10)*SMQ*RPL
     *                            +F(11)*RMQ*SPL
     *                            +F(19)*Z2*T*RPL*SPL)
     *             -F(22)*(Z2D9*FAL(3,22)-Z1D9*FAL(3,10))
     *             -F(23)*(Z2D9*FAL(3,11)-Z1D9*FAL(3,23))
     *             -F(31)*(Z2D9*FAL(3,19)-Z1D9*FAL(3,31))
          FAL(1, 4)=Z1D8*(-SPL*TMI+F(11)*Z2*R*SPL*TMI
     *                            +F(12)*SMQ*TMI
     *                            +F(20)*TMQ*SPL)
     *             -F(23)*(Z2D9*FAL(1,23)-Z1D9*FAL(1,11))
     *             -F(24)*(Z2D9*FAL(1,12)-Z1D9*FAL(1,24))
     *             -F(32)*(Z2D9*FAL(1,20)-Z1D9*FAL(1,32))
          FAL(2, 4)=Z1D8*( RMI*TMI-F(11)*RMQ*TMI
     *                            +F(12)*Z2*S*RMI*TMI
     *                            -F(20)*TMQ*RMI)
     *             -F(23)*(Z2D9*FAL(2,23)-Z1D9*FAL(2,11))
     *             -F(24)*(Z2D9*FAL(2,12)-Z1D9*FAL(2,24))
     *             -F(32)*(Z2D9*FAL(2,20)-Z1D9*FAL(2,32))
          FAL(3, 4)=Z1D8*(-RMI*SPL+F(11)*RMQ*SPL
     *                            +F(12)*SMQ*RMI
     *                            +F(20)*Z2*T*RMI*SPL)
     *             -F(23)*(Z2D9*FAL(3,23)-Z1D9*FAL(3,11))
     *             -F(24)*(Z2D9*FAL(3,12)-Z1D9*FAL(3,24))
     *             -F(32)*(Z2D9*FAL(3,20)-Z1D9*FAL(3,32))
          FAL(1, 5)=Z1D8*(-SMI*TPL+F(13)*Z2*R*SMI*TPL
     *                            +F(16)*SMQ*TPL
     *                            +F(17)*TMQ*SMI)
     *             -F(25)*(Z2D9*FAL(1,13)-Z1D9*FAL(1,25))
     *             -F(28)*(Z2D9*FAL(1,28)-Z1D9*FAL(1,16))
     *             -F(29)*(Z2D9*FAL(1,29)-Z1D9*FAL(1,17))
          FAL(2, 5)=Z1D8*(-RMI*TPL+F(13)*RMQ*TPL
     *                            +F(16)*Z2*S*RMI*TPL
     *                            +F(17)*TMQ*RMI)
     *             -F(25)*(Z2D9*FAL(2,13)-Z1D9*FAL(2,25))
     *             -F(28)*(Z2D9*FAL(2,28)-Z1D9*FAL(2,16))
     *             -F(29)*(Z2D9*FAL(2,29)-Z1D9*FAL(2,17))
          FAL(3, 5)=Z1D8*( RMI*SMI-F(13)*RMQ*SMI
     *                            -F(16)*SMQ*RMI
     *                            +F(17)*Z2*T*RMI*SMI)
     *             -F(25)*(Z2D9*FAL(3,13)-Z1D9*FAL(3,25))
     *             -F(28)*(Z2D9*FAL(3,28)-Z1D9*FAL(3,16))
     *             -F(29)*(Z2D9*FAL(3,29)-Z1D9*FAL(3,17))
          FAL(1, 6)=Z1D8*( SMI*TPL+F(13)*Z2*R*SMI*TPL
     *                            -F(14)*SMQ*TPL
     *                            -F(18)*TMQ*SMI)
     *             -F(25)*(Z2D9*FAL(1,25)-Z1D9*FAL(1,13))
     *             -F(26)*(Z2D9*FAL(1,14)-Z1D9*FAL(1,26))
     *             -F(30)*(Z2D9*FAL(1,30)-Z1D9*FAL(1,18))
          FAL(2, 6)=Z1D8*(-RPL*TPL+F(13)*RMQ*TPL
     *                            +F(14)*Z2*S*RPL*TPL
     *                            +F(18)*TMQ*RPL)
     *             -F(25)*(Z2D9*FAL(2,25)-Z1D9*FAL(2,13))
     *             -F(26)*(Z2D9*FAL(2,14)-Z1D9*FAL(2,26))
     *             -F(30)*(Z2D9*FAL(2,30)-Z1D9*FAL(2,18))
          FAL(3, 6)=Z1D8*( RPL*SMI-F(13)*RMQ*SMI
     *                            -F(14)*SMQ*RPL
     *                            +F(18)*Z2*T*RPL*SMI)
     *             -F(25)*(Z2D9*FAL(3,25)-Z1D9*FAL(3,13))
     *             -F(26)*(Z2D9*FAL(3,14)-Z1D9*FAL(3,26))
     *             -F(30)*(Z2D9*FAL(3,30)-Z1D9*FAL(3,18))
          FAL(1, 7)=Z1D8*( SPL*TPL-F(14)*SMQ*TPL
     *                            +F(15)*Z2*R*SPL*TPL
     *                            -F(19)*TMQ*SPL)
     *             -F(26)*(Z2D9*FAL(1,26)-Z1D9*FAL(1,14))
     *             -F(27)*(Z2D9*FAL(1,15)-Z1D9*FAL(1,27))
     *             -F(31)*(Z2D9*FAL(1,31)-Z1D9*FAL(1,19))
          FAL(2, 7)=Z1D8*( RPL*TPL+F(14)*Z2*S*RPL*TPL
     *                            -F(15)*RMQ*TPL
     *                            -F(19)*TMQ*RPL)
     *             -F(26)*(Z2D9*FAL(2,26)-Z1D9*FAL(2,14))
     *             -F(27)*(Z2D9*FAL(2,15)-Z1D9*FAL(2,27))
     *             -F(31)*(Z2D9*FAL(2,31)-Z1D9*FAL(2,19))
          FAL(3, 7)=Z1D8*( RPL*SPL-F(14)*SMQ*RPL
     *                            -F(15)*RMQ*SPL
     *                            +F(19)*Z2*T*RPL*SPL)
     *             -F(26)*(Z2D9*FAL(3,26)-Z1D9*FAL(3,14))
     *             -F(27)*(Z2D9*FAL(3,15)-Z1D9*FAL(3,27))
     *             -F(31)*(Z2D9*FAL(3,31)-Z1D9*FAL(3,19))
          FAL(1, 8)=Z1D8*(-SPL*TPL+F(15)*Z2*R*SPL*TPL
     *                            +F(16)*SMQ*TPL
     *                            +F(20)*TMQ*SPL)
     *             -F(27)*(Z2D9*FAL(1,27)-Z1D9*FAL(1,15))
     *             -F(28)*(Z2D9*FAL(1,16)-Z1D9*FAL(1,28))
     *             -F(32)*(Z2D9*FAL(1,32)-Z1D9*FAL(1,20))
          FAL(2, 8)=Z1D8*( RMI*TPL-F(15)*RMQ*TPL
     *                            +F(16)*Z2*S*RMI*TPL
     *                            -F(20)*TMQ*RMI)
     *             -F(27)*(Z2D9*FAL(2,27)-Z1D9*FAL(2,15))
     *             -F(28)*(Z2D9*FAL(2,16)-Z1D9*FAL(2,28))
     *             -F(32)*(Z2D9*FAL(2,32)-Z1D9*FAL(2,20))
          FAL(3, 8)=Z1D8*( RMI*SPL-F(15)*RMQ*SPL
     *                            -F(16)*SMQ*RMI
     *                            +F(20)*Z2*T*RMI*SPL)
     *             -F(27)*(Z2D9*FAL(3,27)-Z1D9*FAL(3,15))
     *             -F(28)*(Z2D9*FAL(3,16)-Z1D9*FAL(3,28))
     *             -F(32)*(Z2D9*FAL(3,32)-Z1D9*FAL(3,20))
        END IF
      END IF
      R E T U R N
      END
