      SUBROUTINE ITPANZ(IELTYP,NECKE,INGRD,NINPKT,IERR)                 NEU
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
C ITPANZ:
C =======
C
C BERECHNUNG DER ANZAHL DER INTEGRATIONSPUNKTE FUER EINEN BESTIMMTEN
C ELEMENTTYP UND EINEN BESTIMMTEN INTEGRATIONSGRAD
C
C BEI REDLIN ODER REDINT WIRD DER INTEGRATIONSGRAD KORRIGIERT
C
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
C NECKE     I   I                 HOECHSTE BESETZTE KNOTENNUMMER DES ELEMENTES
C INGRD     I  I/O                INTEGRATIONSGRAD
C                                 INPUT:
C                                 GEWUENSCHTER INTEGRATIONSGRAD
C                                 OUTPUT:
C                                 GEW. INTGRD ODER DURCH REDLIN BZW. REDINT
C                                 KORRIGIERTER INTEGRATIONSGRAD.
C NINPKT    I   O                 ANZAHL DER INTEGRATIONSPUNKTE
C IERR      I   O                 FEHLERPARAMETER
C                                 0: ALLES KLAR
C                                 1: FALSCHER ELEMENTTYP
C                                 2: FALSCHER INTEGRATIONSGRAD
C#
C
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
C +++++
C PRUEFE UND KOORIGIERE GGF. DEN INTEGRATIONSGRAD
C +++++
      IF(REDLIN)THEN
C +++++
C LINEARE ELEMENTE
C +++++
        IF(NECKE.EQ.LQZANF(IELNUM+1,ILODIM)-1)THEN
          IF(IELNUM.EQ.0.AND.ILODIM.EQ.3)THEN
            INGRD=IVOTET(INGRD)
          ELSE IF((IELNUM.EQ.0.AND.ILODIM.EQ.2).OR.
     *            (IELNUM.EQ.1.AND.ILODIM.EQ.3))THEN
            INGRD=IVODRE(INGRD)
          ELSE
            INGRD=IVOGAU(INGRD)
          END IF
        END IF
      END IF
      IF(REDINT)THEN
C +++++
C REDUZIERTE INTEGRATION
C NICHT FUER KUBISCHE ELEMENTE DURCHFUEHREN
C +++++
        NQUA=LQZEND(IELNUM+1,ILODIM)-LQZANF(IELNUM+1,ILODIM)+1
        NECMIN=LQZANF(IELNUM+1,ILODIM)+NQUA
        IF(IGLDIM.EQ.ILODIM.AND.NECKE.LT.NECMIN)THEN
          IF(IELNUM.EQ.0.AND.ILODIM.EQ.3)THEN
            INGRD=IVOTET(INGRD)
          ELSE IF((IELNUM.EQ.0.AND.ILODIM.EQ.2).OR.
     *            (IELNUM.EQ.1.AND.ILODIM.EQ.3))THEN
            INGRD=IVODRE(INGRD)
          ELSE
            INGRD=IVOGAU(INGRD)
          END IF
        END IF
      END IF
C
      IF(ILODIM.EQ.1)THEN
C +++++
C 1-DIMENSIONALE ELEMENTE
C +++++
        IF(IELNUM.EQ.0)THEN
C +++++
C STAB
C +++++
          NINPKT=INGRD
        ELSE
          GOTO 100
        END IF
      ELSE IF(ILODIM.EQ.2)THEN
C +++++
C 2-DIMENSIONALE ELEMENTE
C +++++
        IF(IELNUM.EQ.0)THEN
C +++++
C DREIECK
C +++++
          IF(INGRD.EQ.1)THEN
            NINPKT=1
          ELSE IF(INGRD.EQ.2.OR.INGRD.EQ.3)THEN
            NINPKT=3
          ELSE IF(INGRD.EQ.4)THEN
            NINPKT=6
          ELSE IF(INGRD.EQ.5)THEN
            NINPKT=7
          ELSE
            NINPKT=(INGRD-5)*(INGRD-5)
          END IF
        ELSE IF(IELNUM.EQ.1)THEN
C +++++
C VIERECK
C +++++
          NINPKT=INGRD*INGRD
        ELSE
          GOTO 100
        END IF
      ELSE IF(ILODIM.EQ.3)THEN
C +++++
C 3-DIMENSIONALE ELEMENTE
C +++++
        IF(IELNUM.EQ.0)THEN
C +++++
C TETRAEDER
C +++++
          IF(INGRD.EQ.1)THEN
            NINPKT=1
          ELSE IF(INGRD.EQ.2)THEN
            NINPKT=4
          ELSE IF(INGRD.EQ.3)THEN
            NINPKT=5
          ELSE
            GOTO 200
          END IF
        ELSE IF(IELNUM.EQ.1)THEN
C +++++
C KEIL
C +++++
          IF(INGRD.EQ.1)THEN
            NINPKT=1
          ELSE IF(INGRD.EQ.2.OR.INGRD.EQ.3)THEN
            NINPKT=6
          ELSE IF(INGRD.EQ.4)THEN
            NINPKT=12
          ELSE IF(INGRD.EQ.5)THEN
            NINPKT=21
          ELSE
            NINPKT=(INGRD-5)*(INGRD-5)*(INGRD-5)
          END IF
        ELSE IF(IELNUM.EQ.2)THEN
C +++++
C QUADER
C +++++
          NINPKT=INGRD*INGRD*INGRD
        ELSE
          GOTO 100
        END IF
      ELSE
        GOTO 100
      END IF
      IERR=0
      R E T U R N
  100 IERR=1
      R E T U R N
  200 IERR=2
      R E T U R N
      END
