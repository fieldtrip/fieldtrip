      SUBROUTINE LOKNKO(IELTYP,LKNEL,XLOK,IERR)                         NEU
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
C LOKNKO:
C =======
C
C BERECHNUNG DER LOKALEN KOORDINATEN DER KNOTEN DES ELEMENTES MIT DEM
C ELEMENTTYP IELTYP
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
C LKNEL     I   I  (MKNPEP)       KNOTEN-ELEMENT-ZUORDNUNG UNCOMPRESSED
C XLOK      D   O  (3,MKNPEP)     LOKALE KNOTENKOORDINATEN R,S,T
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
      DIMENSION  LKNEL(  MKNPEP),XLOK  (3,MKNPEP)
C
C +++++
C ELEMENTTYP DEKODIEREN
C +++++
      WRITE(CHTYP,'(I5)')IELTYP
      READ(CHTYP,'(BN,2X,3I1)')IGLDIM,IELNUM,ILODIM
C
      MECKE=MKNEL (IELNUM+1,ILODIM)
      IQANF=LQZANF(IELNUM+1,ILODIM)
      IQEND=LQZEND(IELNUM+1,ILODIM)
      ISUMQ=0
      ISUMK=0
      DO 10 IKNEL=IQANF,IQEND
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMQ=ISUMQ+1
      ELSE
        F(IKNEL)=Z0
      END IF
   10 CONTINUE
      DO 20 IKNEL=IQEND+1,MECKE
      IF(LKNEL(IKNEL).NE.0)THEN
        F(IKNEL)=Z1
        ISUMK=ISUMK+1
      ELSE
        F(IKNEL)=Z0
      END IF
   20 CONTINUE
      IF(ILODIM.EQ.1)THEN
C +++++
C 1-DIMENSIONALE ELEMENTE
C +++++
        IF(IELNUM.EQ.0)THEN
C +++++
C STAB
C +++++
          DO 30 IECKE=1,MECKE
          XLOK(1,IECKE)=DLKSTA(IECKE)
   30     CONTINUE
          IF(ISUMQ.NE.0.AND.ISUMK.NE.0)THEN
            DO 40 IECKE=1,MECKE-IQEND
            JECKE=IQEND+IECKE
            IQECKE=LQKSTA(IECKE)
            XLOK(1,IQECKE)=XLOK(1,IQECKE)+F(JECKE)*DQKSTA(IECKE)
   40       CONTINUE
          END IF
        ELSE
          GOTO 1000
        END IF
      ELSE IF(ILODIM.EQ.2)THEN
C +++++
C 2-DIMENSIONALE ELEMENTE
C +++++
        IF(IELNUM.EQ.0)THEN
C +++++
C DREIECK
C +++++
          DO 50 IECKE=1,MECKE
          XLOK(1,IECKE)=DLKDRE(1,IECKE)
          XLOK(2,IECKE)=DLKDRE(2,IECKE)
   50     CONTINUE
          IF(ISUMQ.NE.0.AND.ISUMK.NE.0)THEN
            DO 60 IECKE=1,MECKE-IQEND
            JECKE=IQEND+IECKE
            IQECKE=LQKDRE(IECKE)
            XLOK(1,IQECKE)=XLOK(1,IQECKE)+F(JECKE)*DQKDRE(1,IECKE)
            XLOK(2,IQECKE)=XLOK(2,IQECKE)+F(JECKE)*DQKDRE(2,IECKE)
   60       CONTINUE
          END IF
        ELSE IF(IELNUM.EQ.1)THEN
C +++++
C VIERECK
C +++++
          DO 70 IECKE=1,MECKE
          XLOK(1,IECKE)=DLKVIE(1,IECKE)
          XLOK(2,IECKE)=DLKVIE(2,IECKE)
   70     CONTINUE
          IF(ISUMQ.NE.0.AND.ISUMK.NE.0)THEN
            DO 80 IECKE=1,MECKE-IQEND
            JECKE=IQEND+IECKE
            IQECKE=LQKVIE(IECKE)
            XLOK(1,IQECKE)=XLOK(1,IQECKE)+F(JECKE)*DQKVIE(1,IECKE)
            XLOK(2,IQECKE)=XLOK(2,IQECKE)+F(JECKE)*DQKVIE(2,IECKE)
   80       CONTINUE
          END IF
        ELSE
          GOTO 1000
        END IF
      ELSE IF(ILODIM.EQ.3)THEN
C +++++
C 3-DIMENSIONALE ELEMENTE
C +++++
        IF(IELNUM.EQ.0)THEN
C +++++
C TETRAEDER
C +++++
          DO 90 IECKE=1,MECKE
          XLOK(1,IECKE)=DLKTET(1,IECKE)
          XLOK(2,IECKE)=DLKTET(2,IECKE)
          XLOK(3,IECKE)=DLKTET(3,IECKE)
   90     CONTINUE
          IF(ISUMQ.NE.0.AND.ISUMK.NE.0)THEN
            DO 100 IECKE=1,MECKE-IQEND
            JECKE=IQEND+IECKE
            IQECKE=LQKTET(IECKE)
            XLOK(1,IQECKE)=XLOK(1,IQECKE)+F(JECKE)*DQKTET(1,IECKE)
            XLOK(2,IQECKE)=XLOK(2,IQECKE)+F(JECKE)*DQKTET(2,IECKE)
            XLOK(3,IQECKE)=XLOK(3,IQECKE)+F(JECKE)*DQKTET(3,IECKE)
  100       CONTINUE
          END IF
        ELSE IF(IELNUM.EQ.1)THEN
C +++++
C KEIL
C +++++
          DO 110 IECKE=1,MECKE
          XLOK(1,IECKE)=DLKKEI(1,IECKE)
          XLOK(2,IECKE)=DLKKEI(2,IECKE)
          XLOK(3,IECKE)=DLKKEI(3,IECKE)
  110     CONTINUE
          IF(ISUMQ.NE.0.AND.ISUMK.NE.0)THEN
            DO 120 IECKE=1,MECKE-IQEND
            JECKE=IQEND+IECKE
            IQECKE=LQKKEI(IECKE)
            XLOK(1,IQECKE)=XLOK(1,IQECKE)+F(JECKE)*DQKKEI(1,IECKE)
            XLOK(2,IQECKE)=XLOK(2,IQECKE)+F(JECKE)*DQKKEI(2,IECKE)
            XLOK(3,IQECKE)=XLOK(3,IQECKE)+F(JECKE)*DQKKEI(3,IECKE)
  120       CONTINUE
          END IF
        ELSE IF(IELNUM.EQ.2)THEN
C +++++
C QUADER
C +++++
          DO 130 IECKE=1,MECKE
          XLOK(1,IECKE)=DLKQUA(1,IECKE)
          XLOK(2,IECKE)=DLKQUA(2,IECKE)
          XLOK(3,IECKE)=DLKQUA(3,IECKE)
  130     CONTINUE
          IF(ISUMQ.NE.0.AND.ISUMK.NE.0)THEN
            DO 140 IECKE=1,MECKE-IQEND
            JECKE=IQEND+IECKE
            IQECKE=LQKQUA(IECKE)
            XLOK(1,IQECKE)=XLOK(1,IQECKE)+F(JECKE)*DQKQUA(1,IECKE)
            XLOK(2,IQECKE)=XLOK(2,IQECKE)+F(JECKE)*DQKQUA(2,IECKE)
            XLOK(3,IQECKE)=XLOK(3,IQECKE)+F(JECKE)*DQKQUA(3,IECKE)
  140       CONTINUE
          END IF
        ELSE
          GOTO 1000
        END IF
      ELSE
        GOTO 1000
      END IF
      IERR=0
      R E T U R N
 1000 IERR=1
      R E T U R N
      END
