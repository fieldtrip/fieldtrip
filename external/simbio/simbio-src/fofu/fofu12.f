      SUBROUTINE FOFU12(IBEWER,IBEABL,R,S,LKNEL ,FWE   ,FAL   )         NEU
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
C FOFU12:
C =======
C
C BERECHNUNG DER WERTE UND DER LOKALEN ABLEITUNGEN DER FORMFUNKTIONEN
C DES ISOPARAMETRISCHEN VIERECKES (ELEMENTTYP XXX12)
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
      DO 10 IKNEL=5,8
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMQ=ISUMQ+1
      ELSE
        F(IKNEL)=Z0
      END IF
   10 CONTINUE
      DO 20 IKNEL=9,12
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMK=ISUMK+1
      ELSE
        F(IKNEL)=Z0
      END IF
   20 CONTINUE
      RMI=Z1-R
      RPL=Z1+R
      RMQ=Z1-R*R
      SMI=Z1-S
      SPL=Z1+S
      SMQ=Z1-S*S
      S3M=Z1-Z3*S
      S3P=Z1+Z3*S
      R3M=Z1-Z3*R
      R3P=Z1+Z3*R
      RSQ=(Z9*R*R + Z9*S*S - Z10)
      IF(ISUMQ.EQ.4.AND.ISUMK.EQ.4)THEN
C +++++
C VOLL-KUBISCHES VIERECK
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z1D32*RMI*SMI*RSQ
          FWE( 2)=Z1D32*RPL*SMI*RSQ
          FWE( 3)=Z1D32*RPL*SPL*RSQ
          FWE( 4)=Z1D32*RMI*SPL*RSQ
          FWE( 5)=Z9D32*RMQ*R3M*SMI
          FWE( 6)=Z9D32            *SMQ*S3M*RPL
          FWE( 7)=Z9D32*RMQ*R3P*SPL
          FWE( 8)=Z9D32            *SMQ*S3P*RMI
          FWE( 9)=Z9D32*RMQ*R3P*SMI
          FWE(10)=Z9D32            *SMQ*S3P*RPL
          FWE(11)=Z9D32*RMQ*R3M*SPL
          FWE(12)=Z9D32            *SMQ*S3M*RMI
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1, 1)= Z9D32*SMI*(Z3*RMQ+Z2*R+SMQ-Z26D9)
          FAL(2, 1)= Z9D32*RMI*(Z3*SMQ+Z2*S+RMQ-Z26D9)
          FAL(1, 2)=-Z9D32*SMI*(Z3*RMQ-Z2*R+SMQ-Z26D9)
          FAL(2 ,2)= Z9D32*RPL*(Z3*SMQ+Z2*S+RMQ-Z26D9)
          FAL(1, 3)=-Z9D32*SPL*(Z3*RMQ-Z2*R+SMQ-Z26D9)
          FAL(2, 3)=-Z9D32*RPL*(Z3*SMQ-Z2*S+RMQ-Z26D9)
          FAL(1, 4)= Z9D32*SPL*(Z3*RMQ+Z2*R+SMQ-Z26D9)
          FAL(2, 4)=-Z9D32*RMI*(Z3*SMQ-Z2*S+RMQ-Z26D9)
          FAL(1, 5)=-Z9D32*SMI*(R3P*R3M+Z2*RPL)
          FAL(2, 5)=-Z9D32*RMQ*R3M
          FAL(1, 6)= Z9D32*SMQ*S3M
          FAL(2, 6)=-Z9D32*RPL*(S3P*S3M+Z2*SPL)
          FAL(1, 7)= Z9D32*SPL*(R3P*R3M+Z2*RMI)
          FAL(2, 7)= Z9D32*RMQ*R3P
          FAL(1, 8)=-Z9D32*SMQ*S3P
          FAL(2, 8)= Z9D32*RMI*(S3P*S3M+Z2*SMI)
          FAL(1, 9)= Z9D32*SMI*(R3P*R3M+Z2*RMI)
          FAL(2, 9)=-Z9D32*RMQ*R3P
          FAL(1,10)= Z9D32*SMQ*S3P
          FAL(2,10)= Z9D32*RPL*(S3P*S3M+Z2*SMI)
          FAL(1,11)=-Z9D32*SPL*(R3P*R3M+Z2*RPL)
          FAL(2,11)= Z9D32*RMQ*R3M
          FAL(1,12)=-Z9D32*SMQ*S3M
          FAL(2,12)=-Z9D32*RMI*(S3P*S3M+Z2*SPL)
        END IF
      ELSE IF(ISUMQ.EQ.4.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-QUADRATISCHES VIERECK
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(1)=Z025*RMI*SMI*(-Z1-R-S)
          FWE(2)=Z025*RPL*SMI*(-Z1+R-S)
          FWE(3)=Z025*RPL*SPL*(-Z1+R+S)
          FWE(4)=Z025*RMI*SPL*(-Z1-R+S)
          FWE(5)=Z05 *RMQ*SMI
          FWE(6)=Z05 *RPL*SMQ
          FWE(7)=Z05 *RMQ*SPL
          FWE(8)=Z05 *RMI*SMQ
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,1)=Z025*SMI*( Z2*R+   S)
          FAL(2,1)=Z025*RMI*(    R+Z2*S)
          FAL(1,2)=Z025*SMI*( Z2*R-   S)
          FAL(2,2)=Z025*RPL*(-   R+Z2*S)
          FAL(1,3)=Z025*SPL*( Z2*R+   S)
          FAL(2,3)=Z025*RPL*(    R+Z2*S)
          FAL(1,4)=Z025*SPL*( Z2*R-   S)
          FAL(2,4)=Z025*RMI*(-   R+Z2*S)
          FAL(1,5)=-    R*SMI
          FAL(2,5)=-Z05  *RMQ
          FAL(1,6)= Z05  *SMQ
          FAL(2,6)=-    S*RPL
          FAL(1,7)=-    R*SPL
          FAL(2,7)= Z05  *RMQ
          FAL(1,8)=-Z05  *SMQ
          FAL(2,8)=-    S*RMI
        END IF
      ELSE IF(ISUMQ.EQ.0.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-LINEARES VIERECK
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(1)= Z025*RMI*SMI
          FWE(2)= Z025*RPL*SMI
          FWE(3)= Z025*RPL*SPL
          FWE(4)= Z025*RMI*SPL
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,1)=-Z025*SMI
          FAL(2,1)=-Z025*RMI
          FAL(1,2)= Z025*SMI
          FAL(2,2)=-Z025*RPL
          FAL(1,3)= Z025*SPL
          FAL(2,3)= Z025*RPL
          FAL(1,4)=-Z025*SPL
          FAL(2,4)= Z025*RMI
        END IF
      ELSE
C +++++
C VIERECK MIT VARIABLER KNOTENZAHL
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 9)=F( 9)*Z9D32*RMQ*R3P*SMI
          FWE(10)=F(10)*Z9D32            *SMQ*S3P*RPL
          FWE(11)=F(11)*Z9D32*RMQ*R3M*SPL
          FWE(12)=F(12)*Z9D32            *SMQ*S3M*RMI

          FWE( 5)=F(5)*(Z05*SMI*(RMQ+F( 9)*Z0125*RMQ) - FWE( 9))
          FWE( 6)=F(6)*(Z05*RPL*(SMQ+F(10)*Z0125*SMQ) - FWE(10))
          FWE( 7)=F(7)*(Z05*SPL*(RMQ+F(11)*Z0125*RMQ) - FWE(11))
          FWE( 8)=F(8)*(Z05*RMI*(SMQ+F(12)*Z0125*SMQ) - FWE(12))

          FWE( 1)=Z025*RMI*SMI-F(5)*Z025*RMQ*SMI-F(8)*Z025*SMQ*RMI
     *           -F( 9)*(Z2D9*FWE( 5)-Z1D9*FWE( 9))
     *           -F(12)*(Z2D9*FWE(12)-Z1D9*FWE( 8))
          FWE( 2)=Z025*RPL*SMI-F(6)*Z025*SMQ*RPL-F(5)*Z025*RMQ*SMI
     *           -F(10)*(Z2D9*FWE( 6)-Z1D9*FWE(10))
     *           -F( 9)*(Z2D9*FWE( 9)-Z1D9*FWE( 5))
          FWE( 3)=Z025*RPL*SPL-F(7)*Z025*RMQ*SPL-F(6)*Z025*SMQ*RPL
     *           -F(11)*(Z2D9*FWE( 7)-Z1D9*FWE(11))
     *           -F(10)*(Z2D9*FWE(10)-Z1D9*FWE( 6))
          FWE( 4)=Z025*RMI*SPL-F(8)*Z025*SMQ*RMI-F(7)*Z025*RMQ*SPL
     *           -F(12)*(Z2D9*FWE( 8)-Z1D9*FWE(12))
     *           -F(11)*(Z2D9*FWE(11)-Z1D9*FWE( 7))
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1, 9)= F( 9)*Z9D32*SMI*(R3P*R3M+Z2*RMI)
          FAL(2, 9)=-F( 9)*Z9D32*RMQ*R3P
          FAL(1,10)= F(10)*Z9D32*SMQ*S3P
          FAL(2,10)= F(10)*Z9D32*RPL*(S3P*S3M+Z2*SMI)
          FAL(1,11)=-F(11)*Z9D32*SPL*(R3P*R3M+Z2*RPL)
          FAL(2,11)= F(11)*Z9D32*RMQ*R3M
          FAL(1,12)=-F(12)*Z9D32*SMQ*S3M
          FAL(2,12)=-F(12)*Z9D32*RMI*(S3P*S3M+Z2*SPL)
          FAL(1, 5)=-F( 5)*(SMI*(R  +F( 9)*Z0125*R  ) + FAL(1, 9))
          FAL(2, 5)=-F( 5)*(Z05*(RMQ+F( 9)*Z0125*RMQ) + FAL(2, 9))
          FAL(1, 6)= F( 6)*(Z05*(SMQ+F(10)*Z0125*SMQ) - FAL(1,10))
          FAL(2, 6)=-F( 6)*(RPL*(S  +F(10)*Z0125*S  ) + FAL(2,10))
          FAL(1, 7)=-F( 7)*(SPL*(R  +F(11)*Z0125*R  ) + FAL(1,11))
          FAL(2, 7)= F( 7)*(Z05*(RMQ+F(11)*Z0125*RMQ) - FAL(2,11))
          FAL(1, 8)=-F( 8)*(Z05*(SMQ+F(12)*Z0125*SMQ) + FAL(1,12))
          FAL(2, 8)=-F( 8)*(RMI*(S  +F(12)*Z0125*S  ) + FAL(2,12))
          FAL(1, 1)=-Z025*SMI + F(5)*Z05*R*SMI + F(8)*Z025*SMQ
     *              -F( 9)*(Z2D9*FAL(1, 5)-Z1D9*FAL(1,9))
     *              -F(12)*(Z2D9*FAL(1,12)-Z1D9*FAL(1, 8))
          FAL(2, 1)=-Z025*RMI + F(5)*Z025*RMQ + F(8)*Z05*S*RMI
     *              -F( 9)*(Z2D9*FAL(2, 5)-Z1D9*FAL(2, 9))
     *              -F(12)*(Z2D9*FAL(2,12)-Z1D9*FAL(2, 8))
          FAL(1, 2)= Z025*SMI + F(5)*Z05*R*SMI - F(6)*Z025*SMQ
     *              -F( 9)*(Z2D9*FAL(1, 9) - Z1D9*FAL(1, 5))
     *              -F(10)*(Z2D9*FAL(1, 6) - Z1D9*FAL(1,10))
          FAL(2, 2)=-Z025*RPL + F(5)*Z025*RMQ + F(6)*Z05*S*RPL
     *              -F( 9)*(Z2D9*FAL(2, 9) - Z1D9*FAL(2, 5))
     *              -F(10)*(Z2D9*FAL(2, 6) - Z1D9*FAL(2,10))
          FAL(1, 3)= Z025*SPL - F(6)*Z025*SMQ + F(7)*Z05*R*SPL
     *              -F(10)*(Z2D9*FAL(1,10) - Z1D9*FAL(1, 6))
     *              -F(11)*(Z2D9*FAL(1, 7) - Z1D9*FAL(1,11))
          FAL(2, 3)= Z025*RPL + F(6)*Z05*S*RPL - F(7)*Z025*RMQ
     *              -F(10)*(Z2D9*FAL(2,10) - Z1D9*FAL(2, 6))
     *              -F(11)*(Z2D9*FAL(2, 7) - Z1D9*FAL(2,11))
          FAL(1, 4)=-Z025*SPL + F(7)*Z05*R*SPL + F(8)*Z025*SMQ
     *              -F(11)*(Z2D9*FAL(1,11) - Z1D9*FAL(1, 7))
     *              -F(12)*(Z2D9*FAL(1, 8) - Z1D9*FAL(1,12))
          FAL(2, 4)= Z025*RMI - F(7)*Z025*RMQ + F(8)*Z05*S*RMI
     *              -F(11)*(Z2D9*FAL(2,11) - Z1D9*FAL(2, 7))
     *              -F(12)*(Z2D9*FAL(2, 8) - Z1D9*FAL(2,12))
        END IF
      END IF
      R E T U R N
      END
