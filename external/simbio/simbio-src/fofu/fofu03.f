      SUBROUTINE FOFU03(IBEWER,IBEABL,R,S,T,LKNEL ,FWE   ,FAL   )       NEU
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
C FOFU03:
C =======
C
C BERECHNUNG DER WERTE UND DER LOKALEN ABLEITUNGEN DER FORMFUNKTIONEN
C DES ISOPARAMETRISCHEN TETRAEDERS (ELEMENTTYP XXX03)
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
      DO 10 IKNEL=5,10
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMQ=ISUMQ+1
      ELSE
        F(IKNEL)=Z0
      END IF
   10 CONTINUE
      DO 20 IKNEL=11,16
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMK=ISUMK+1
      ELSE
        F(IKNEL)=Z0
      END IF
   20 CONTINUE
      U=Z1-R-S-T
      RQ=R*R
      SQ=S*S
      TQ=T*T
      UQ=U*U
      RST = Z2*R*S*T
      RSU = Z2*R*S*U
      STU = Z2*S*T*U
      RTU = Z2*R*T*U
      IF(ISUMQ.EQ.6.AND.ISUMK.EQ.6)THEN
C +++++
C VOLL-KUBISCHES TETRAEDER
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=Z9D2*R*(Z1D3-R)*(Z2D3-R)
          FWE( 2)=Z9D2*S*(Z1D3-S)*(Z2D3-S)
          FWE( 3)=Z9D2*T*(Z1D3-T)*(Z2D3-T)
          FWE( 4)=Z9D2*U*(Z1D3-U)*(Z2D3-U)
          FWE( 5)=Z9D2*R*S*(Z2*R-S)
          FWE( 6)=Z9D2*S*T*(Z2*S-T)
          FWE( 7)=Z9D2*R*T*(Z2*T-R)
          FWE( 8)=Z9D2*R*U*(Z2*R-U)
          FWE( 9)=Z9D2*S*U*(Z2*S-U)
          FWE(10)=Z9D2*T*U*(Z2*T-U)
          FWE(11)=Z9D2*R*S*(Z2*S-R)
          FWE(12)=Z9D2*S*T*(Z2*T-S)
          FWE(13)=Z9D2*R*T*(Z2*R-T)
          FWE(14)=Z9D2*R*U*(Z2*U-R)
          FWE(15)=Z9D2*S*U*(Z2*U-S)
          FWE(16)=Z9D2*T*U*(Z2*U-T)
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1, 1)=Z9D2*(Z2D9-Z2*R+Z3*RQ)
          FAL(2, 1)=Z0
          FAL(3, 1)=Z0
          FAL(1, 2)=Z0
          FAL(2, 2)=Z9D2*(Z2D9-Z2*S+Z3*SQ)
          FAL(3, 2)=Z0
          FAL(1, 3)=Z0
          FAL(2, 3)=Z0
          FAL(3, 3)=Z9D2*(Z2D9-Z2*T+Z3*TQ)
          FAL(1, 4)=-Z9D2*(Z2D9-Z2*U+Z3*UQ)
          FAL(2, 4)= FAL(1, 4)
          FAL(3, 4)= FAL(1, 4)
          FAL(1, 5)=Z9D2*(Z4*R*S-SQ)
          FAL(2, 5)=Z9  *(RQ-R*S)
          FAL(3, 5)=Z0
          FAL(1, 6)=Z0
          FAL(2, 6)=Z9D2*(Z4*S*T-TQ)
          FAL(3, 6)=Z9  *(SQ-S*T)
          FAL(1, 7)=Z9  *(TQ-R*T)
          FAL(2, 7)=Z0
          FAL(3, 7)=Z9D2*(Z4*T*R-RQ)
          FAL(1, 8)=Z9D2*((Z5*R-U)*(U-R)+Z3*RQ)
          FAL(2, 8)=-Z9 *(RQ-R*U)
          FAL(3, 8)= FAL(2, 8)
          FAL(1, 9)=-Z9 *(SQ-S*U)
          FAL(2, 9)=Z9D2*((Z5*S-U)*(U-S)+Z3*SQ)
          FAL(3, 9)= FAL(1, 9)
          FAL(1,10)=-Z9 *(TQ-T*U)
          FAL(2,10)= FAL(1,10)
          FAL(3,10)=Z9D2*((Z5*T-U)*(U-T)+Z3*TQ)
          FAL(1,11)=Z9  *(SQ-R*S)
          FAL(2,11)=Z9D2*(Z4*R*S-RQ)
          FAL(3,11)=Z0
          FAL(1,12)=Z0
          FAL(2,12)=Z9  *(TQ-S*T)
          FAL(3,12)=Z9D2*(Z4*S*T-SQ)
          FAL(1,13)=Z9D2*(Z4*R*T-TQ)
          FAL(2,13)=Z0
          FAL(3,13)=Z9  *(RQ-R*T)
          FAL(1,14)=Z9*((U-R)*(U-Z2*R)-Z3D2*RQ)
          FAL(2,14)=-Z9D2*(Z4*U*R-RQ)
          FAL(3,14)= FAL(2,14)
          FAL(1,15)=-Z9D2*(Z4*U*S-SQ)
          FAL(2,15)=Z9*((U-S)*(U-Z2*S)-Z3D2*SQ)
          FAL(3,15)= FAL(1,15)
          FAL(1,16)=-Z9D2*(Z4*U*T-TQ)
          FAL(2,16)= FAL(1,16)
          FAL(3,16)=Z9*((U-T)*(U-Z2*T)-Z3D2*TQ)
        END IF
      ELSE IF(ISUMQ.EQ.6.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-QUADRATISCHES TETRAEDER
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=(Z2*R-Z1)*R
          FWE( 2)=(Z2*S-Z1)*S
          FWE( 3)=(Z2*T-Z1)*T
          FWE( 4)=(Z2*U-Z1)*U
          FWE( 5)=Z4*R*S
          FWE( 6)=Z4*S*T
          FWE( 7)=Z4*R*T
          FWE( 8)=Z4*R*U
          FWE( 9)=Z4*S*U
          FWE(10)=Z4*T*U
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL( 1, 1)= Z4*R-Z1
          FAL( 2, 1)= Z0
          FAL( 3, 1)= Z0
          FAL( 1, 2)= Z0
          FAL( 2, 2)= Z4*S-Z1
          FAL( 3, 2)= Z0
          FAL( 1, 3)= Z0
          FAL( 2, 3)= Z0
          FAL( 3, 3)= Z4*T-Z1
          FAL( 1, 4)= Z1-Z4*U
          FAL( 2, 4)= Z1-Z4*U
          FAL( 3, 4)= Z1-Z4*U
          FAL( 1, 5)= Z4*S
          FAL( 2, 5)= Z4*R
          FAL( 3, 5)= Z0
          FAL( 1, 6)= Z0
          FAL( 2, 6)= Z4*T
          FAL( 3, 6)= Z4*S
          FAL( 1, 7)= Z4*T
          FAL( 2, 7)= Z0
          FAL( 3, 7)= Z4*R
          FAL( 1, 8)= Z4*(U-R)
          FAL( 2, 8)=-Z4*R
          FAL( 3, 8)=-Z4*R
          FAL( 1, 9)=-Z4*S
          FAL( 2, 9)= Z4*(U-S)
          FAL( 3, 9)=-Z4*S
          FAL( 1,10)=-Z4*T
          FAL( 2,10)=-Z4*T
          FAL( 3,10)= Z4*(U-T)
        END IF
      ELSE IF(ISUMQ.EQ.0.AND.ISUMK.EQ.0)THEN
C +++++
C VOLL-LINEARES TETRAEDER
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE( 1)=R
          FWE( 2)=S
          FWE( 3)=T
          FWE( 4)=U
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL( 1, 1)= Z1
          FAL( 2, 1)= Z0
          FAL( 3, 1)= Z0
          FAL( 1, 2)= Z0
          FAL( 2, 2)= Z1
          FAL( 3, 2)= Z0
          FAL( 1, 3)= Z0
          FAL( 2, 3)= Z0
          FAL( 3, 3)= Z1
          FAL( 1, 4)=-Z1
          FAL( 2, 4)=-Z1
          FAL( 3, 4)=-Z1
        END IF
      ELSE
C +++++
C TETRAEDER MIT VARIABLER KNOTENZAHL
C +++++
        IF(IBEWER.EQ.1)THEN
          FWE(11)=F(11)*Z9D2*R*S*(Z2*S-R)
          FWE(12)=F(12)*Z9D2*S*T*(Z2*T-S)
          FWE(13)=F(13)*Z9D2*R*T*(Z2*R-T)
          FWE(14)=F(14)*Z9D2*R*U*(Z2*U-R)
          FWE(15)=F(15)*Z9D2*S*U*(Z2*U-S)
          FWE(16)=F(16)*Z9D2*T*U*(Z2*U-T)
          FWE( 5)=F( 5)*(Z4*R*S-F(11)*(FWE(11)+Z9D2*R*S*(Z8D9-R-S)) )
          FWE( 6)=F( 6)*(Z4*S*T-F(12)*(FWE(12)+Z9D2*S*T*(Z8D9-S-T)) )
          FWE( 7)=F( 7)*(Z4*R*T-F(13)*(FWE(13)+Z9D2*R*T*(Z8D9-R-T)) )
          FWE( 8)=F( 8)*(Z4*R*U-F(14)*(FWE(14)+Z9D2*R*U*(Z8D9-R-U)) )
          FWE( 9)=F( 9)*(Z4*S*U-F(15)*(FWE(15)+Z9D2*S*U*(Z8D9-S-U)) )
          FWE(10)=F(10)*(Z4*T*U-F(16)*(FWE(16)+Z9D2*T*U*(Z8D9-T-U)) )
          FWE( 1)=R-F( 5)*Z2*R*S-F( 7)*Z2*R*T-F( 8)*Z2*R*U
     *             -F(11)*(Z2D9*FWE( 5)-Z1D9*FWE(11)-RST-RSU)
     *             -F(13)*(Z2D9*FWE(13)-Z1D9*FWE( 7)-RST-RTU)
     *             -F(14)*(Z2D9*FWE( 8)-Z1D9*FWE(14)-RTU-RSU)
          FWE( 2)=S-F( 5)*Z2*R*S-F( 6)*Z2*S*T-F( 9)*Z2*S*U
     *             -F(11)*(Z2D9*FWE(11)-Z1D9*FWE( 5)-RST-RSU)
     *             -F(12)*(Z2D9*FWE( 6)-Z1D9*FWE(12)-RST-STU)
     *             -F(15)*(Z2D9*FWE( 9)-Z1D9*FWE(15)-RSU-STU)
          FWE( 3)=T-F( 6)*Z2*S*T-F( 7)*Z2*R*T-F(10)*Z2*U*T
     *             -F(12)*(Z2D9*FWE(12)-Z1D9*FWE( 6)-RST-STU)
     *             -F(13)*(Z2D9*FWE( 7)-Z1D9*FWE(13)-RST-RTU)
     *             -F(16)*(Z2D9*FWE(10)-Z1D9*FWE(16)-RTU-STU)
          FWE( 4)=U-F( 8)*Z2*R*U-F( 9)*Z2*S*U-F(10)*Z2*T*U
     *             -F(14)*(Z2D9*FWE(14)-Z1D9*FWE( 8)-RTU-RSU)
     *             -F(15)*(Z2D9*FWE(15)-Z1D9*FWE( 9)-RSU-STU)
     *             -F(16)*(Z2D9*FWE(16)-Z1D9*FWE(10)-RTU-STU)
        END IF
        IF(IBEABL.EQ.1)THEN
          FAL(1,11)= F(11)*Z9  *(SQ-R*S)
          FAL(2,11)= F(11)*Z9D2*(Z4*R*S-RQ)
          FAL(3,11)=       Z0
          FAL(1,12)=       Z0
          FAL(2,12)= F(12)*Z9  *(TQ-S*T)
          FAL(3,12)= F(12)*Z9D2*(Z4*S*T-SQ)
          FAL(1,13)= F(13)*Z9D2*(Z4*R*T-TQ)
          FAL(2,13)=       Z0
          FAL(3,13)= F(13)*Z9  *(RQ-R*T)
          FAL(1,14)= F(14)*Z9*((U-R)*(U-Z2*R)-Z3D2*RQ)
          FAL(2,14)=-F(14)*Z9D2*(Z4*U*R-RQ)
          FAL(3,14)=        FAL(2,14)
          FAL(1,15)=-F(15)*Z9D2*(Z4*U*S-SQ)
          FAL(2,15)= F(15)*Z9*((U-S)*(U-Z2*S)-Z3D2*SQ)
          FAL(3,15)=        FAL(1,15)
          FAL(1,16)=-F(16)*Z9D2*(Z4*U*T-TQ)
          FAL(2,16)=        FAL(1,16)
          FAL(3,16)= F(16)*Z9*((U-T)*(U-Z2*T)-Z3D2*TQ)
          FAL(1, 5)= F( 5)*(Z4*S-F(11)*(FAL(1,11)
     *                     +Z9D2*S*(Z8D9-Z2*R-S)) )
          FAL(2, 5)= F( 5)*(Z4*R-F(11)*(FAL(2,11)
     *                     +Z9D2*R*(Z8D9-Z2*S-R)) )
          FAL(3, 5)=       Z0
          FAL(1, 6)=       Z0
          FAL(2, 6)= F( 6)*(Z4*T-F(12)*(FAL(2,12)
     *                     +Z9D2*T*(Z8D9-Z2*S-T)) )
          FAL(3, 6)= F( 6)*(Z4*S-F(12)*(FAL(3,12)
     *                     +Z9D2*S*(Z8D9-Z2*T-S)) )
          FAL(1, 7)= F( 7)*(Z4*T-F(13)*(FAL(1,13)
     *                     +Z9D2*T*(Z8D9-Z2*R-T)) )
          FAL(2, 7)=       Z0
          FAL(3, 7)= F( 7)*(Z4*R-F(13)*(FAL(3,13)
     *                     +Z9D2*R*(Z8D9-Z2*T-R)) )
          FAL(1, 8)= F( 8)*(Z4*(U-R)-F(14)*(FAL(1,14)
     *                     +Z9D2*(U-R)*(Z8D9-R-U)) )
          FAL(2, 8)=-F( 8)*(Z4*R+F(14)*(FAL(2,14)
     *                     -Z9D2*R*(Z8D9-Z2*U-R)) )
          FAL(3, 8)= FAL(2, 8)
          FAL(1, 9)=-F( 9)*(Z4*S+F(15)*(FAL(1,15)
     *                     -Z9D2*S*(Z8D9-Z2*U-S)) )
          FAL(2, 9)= F( 9)*(Z4*(U-S)-F(15)*(FAL(2,15)
     *                     +Z9D2*(U-S)*(Z8D9-S-U)) )
          FAL(3, 9)= FAL(1, 9)
          FAL(1,10)=-F(10)*(Z4*T+F(16)*(FAL(1,16)
     *                     -Z9D2*T*(Z8D9-Z2*U-T)) )
          FAL(2,10)= FAL(1,10)
          FAL(3,10)= F(10)*(Z4*(U-T)-F(16)*(FAL(3,16)
     *                     +Z9D2*(U-T)*(Z8D9-T-U)) )
          FAL(1, 1)= Z1-F( 5)*Z2*S-F( 7)*Z2*T-F( 8)*Z2*(U-R)
     *            -F(11)*(Z2D9*FAL(1, 5)-Z1D9*FAL(1,11)-Z2*S*(T+U-R))
     *            -F(13)*(Z2D9*FAL(1,13)-Z1D9*FAL(1, 7)-Z2*T*(U-R+S))
     *            -F(14)*(Z2D9*FAL(1, 8)-Z1D9*FAL(1,14)-Z2*(T+S)*(U-R))
          FAL(2, 1)=   -F( 5)*Z2*R           +F( 8)*Z2*R
     *            -F(11)*(Z2D9*FAL(2, 5)-Z1D9*FAL(2,11)-Z2*R*(U-S+T))
     *            -F(14)*(Z2D9*FAL(2, 8)-Z1D9*FAL(2,14)-Z2*R*(U-S-T))
          FAL(3, 1)=              -F( 7)*Z2*R+F( 8)*Z2*R
     *            -F(13)*(Z2D9*FAL(3,13)-Z1D9*FAL(3, 7)-Z2*R*(U-T+S))
     *            -F(14)*(Z2D9*FAL(3, 8)-Z1D9*FAL(3,14)-Z2*R*(U-T-S))
          FAL(1, 2)=   -F( 5)*Z2*S           +F( 9)*Z2*S
     *            -F(11)*(Z2D9*FAL(1,11)-Z1D9*FAL(1, 5)-Z2*S*(U-R+T))
     *            -F(15)*(Z2D9*FAL(1, 9)-Z1D9*FAL(1,15)-Z2*S*(U-R-T))
          FAL(2, 2)= Z1-F( 5)*Z2*R-F( 6)*Z2*T-F( 9)*Z2*(U-S)
     *            -F(11)*(Z2D9*FAL(2,11)-Z1D9*FAL(2, 5)-Z2*R*(U-S+T))
     *            -F(12)*(Z2D9*FAL(2, 6)-Z1D9*FAL(2,12)-Z2*T*(U-S+R))
     *            -F(15)*(Z2D9*FAL(2, 9)-Z1D9*FAL(2,15)-Z2*(R+T)*(U-S))
          FAL(3, 2)=              -F( 6)*Z2*S+F( 9)*Z2*S
     *            -F(12)*(Z2D9*FAL(3, 6)-Z1D9*FAL(3,12)-Z2*S*(U-T+R))
     *            -F(15)*(Z2D9*FAL(3, 9)-Z1D9*FAL(3,15)-Z2*S*(U-T-R))
           FAL(1, 3)=              -F( 7)*Z2*T+F(10)*Z2*T
     *            -F(13)*(Z2D9*FAL(1, 7)-Z1D9*FAL(1,13)-Z2*T*(U-R+S))
     *            -F(16)*(Z2D9*FAL(1,10)-Z1D9*FAL(1,16)-Z2*T*(U-R-S))
          FAL(2, 3)=   -F( 6)*Z2*T           +F(10)*Z2*T
     *            -F(12)*(Z2D9*FAL(2,12)-Z1D9*FAL(2, 6)-Z2*T*(U-S+R))
     *            -F(16)*(Z2D9*FAL(2,10)-Z1D9*FAL(2,16)-Z2*T*(U-S-R))
          FAL(3, 3)= Z1-F( 6)*Z2*S-F( 7)*Z2*R-F(10)*Z2*(U-T)
     *            -F(12)*(Z2D9*FAL(3,12)-Z1D9*FAL(3, 6)-Z2*S*(U-T+R))
     *            -F(13)*(Z2D9*FAL(3, 7)-Z1D9*FAL(3,13)-Z2*R*(U-T+S))
     *            -F(16)*(Z2D9*FAL(3,10)-Z1D9*FAL(3,16)-Z2*(R+S)*(U-T))
          FAL(1, 4)=-Z1-F( 8)*Z2*(U-R)+F( 9)*Z2*S+F(10)*Z2*T
     *            -F(14)*(Z2D9*FAL(1,14)-Z1D9*FAL(1, 8)-Z2*(S+T)*(U-R))
     *            -F(15)*(Z2D9*FAL(1,15)-Z1D9*FAL(1, 9)-Z2*S*(U-T-R))
     *            -F(16)*(Z2D9*FAL(1,16)-Z1D9*FAL(1,10)-Z2*T*(U-R-S))
          FAL(2, 4)=-Z1+F( 8)*Z2*R-F( 9)*Z2*(U-S)+F(10)*Z2*T
     *            -F(14)*(Z2D9*FAL(2,14)-Z1D9*FAL(2, 8)-Z2*R*(U-S-T))
     *            -F(15)*(Z2D9*FAL(2,15)-Z1D9*FAL(2, 9)-Z2*(R+T)*(U-S))
     *            -F(16)*(Z2D9*FAL(2,16)-Z1D9*FAL(2,10)-Z2*T*(U-R-S))
          FAL(3, 4)=-Z1+F( 8)*Z2*R+F( 9)*Z2*S-F(10)*Z2*(U-T)
     *            -F(14)*(Z2D9*FAL(3,14)-Z1D9*FAL(3, 8)-Z2*R*(U-T-S))
     *            -F(15)*(Z2D9*FAL(3,15)-Z1D9*FAL(3, 9)-Z2*S*(U-T-R))
     *            -F(16)*(Z2D9*FAL(3,16)-Z1D9*FAL(3,10)-Z2*(R+S)*(U-T))
        END IF
      END IF
      R E T U R N
      END
