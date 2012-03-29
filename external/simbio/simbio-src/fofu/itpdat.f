      SUBROUTINE ITPDAT(IELTYP,INGRD,IINPKT,DINPKT,IERR)                NEU
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
C ITPDAT:
C =======
C
C BERECHNUNG DER DATEN DER INTEGRATIONSPUNKTE ( LAGE, WICHTUNG ) FUER
C EINEN BESTIMMTEN ELEMENTTYP UND EINEN BESTIMMTEN INTEGRATIONSGRAD
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
C INGRD     I   I                 INTEGRATIONSGRAD
C IINPKT    I   I                 NUMMER DES INTEGRATIONSPUNKTES
C DINPKT    D   O (4)             DATEN DES INTEGRATIONSPUNKTES
C                                 (1) : 1. LOKALE KOORDINATE
C                                 (2) : 2. LOKALE KOORDINATE (FALLS DA)
C                                 (3) : 3. LOKALE KOORDINATE (FALLS DA)
C                                 (4) : WICHTUNG
C IERR      I   O                 FEHLERPARAMETER
C                                 0: ALLES KLAR
C                                 1: FALSCHER ELEMENTTYP
C                                 2: FALSCHER INTEGRATIONSGRAD
C                                 3: FALSCHE INTP-NUMMER
C
C
C#
C
      INCLUDE 'fofulib.inc'
C
      CHARACTER*6 CHTYP
C
      DIMENSION DINPKT(4)
C
C +++++
C ELEMENTTYP DEKODIEREN
C +++++
      WRITE(CHTYP,'(I5)')IELTYP
      READ(CHTYP,'(BN,2X,3I1)')IGLDIM,IELNUM,ILODIM
C
      IF(ILODIM.EQ.1)THEN
C +++++
C 1-DIMENSIONALE ELEMENTE
C +++++
        IF(IELNUM.EQ.0)THEN
C +++++
C STAB
C +++++
C
C          NUMERISCHE    STUETZ-      INTEGRATIONS-
C INGRD     ORDNUNG      PUNKTZAHL       ART
C --------------------------------------------------
C  1           1             1          GAUSS
C  2           3             2          GAUSS
C  3           5             3          GAUSS
C  4           7             4          GAUSS
C  5           9             5          GAUSS
C  6          11             6          GAUSS
C  7          13             7          GAUSS
C  8          15             8          GAUSS
C  9          17             9          GAUSS
C
          DINPKT(1)=VIEKOO(IINPKT,INGRD)
          DINPKT(4)=VIEWIC(IINPKT,INGRD)
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
C
C          NUMERISCHE    STUETZ-      INTEGRATIONS-
C INGRD     ORDNUNG      PUNKTZAHL       ART
C --------------------------------------------------
C  1         1              1           DIREKT
C  2         2              3           DIREKT
C  3         2              3           DIREKT
C  4         4              6           DIREKT
C  5         5              7           DIREKT
C  6       1 X 1          1 X 1       GAUSS-RADAU
C  7       3 X 3          2 X 2       GAUSS-RADAU
C  8       5 X 5          3 X 3       GAUSS-RADAU
C  9       7 X 7          4 X 4       GAUSS-RADAU
C 10       9 X 9          5 X 5       GAUSS-RADAU
C
C BEI DER DIREKTEN INTEGRATION NACH HAMMER MUSS FUER DIE EXAKTE
C INTEGRATION DES POLYNOMS XSI**I * ETA**J GELTEN:
C          (I+J)   KLEINER GLEICH   (NUMERISCHE ORDNUNG)
C
          IF(INGRD.EQ.1)THEN
            IANF=0
            IDIR=1
          ELSE IF(INGRD.EQ.2)THEN
            IANF=1
            IDIR=1
          ELSE IF(INGRD.EQ.3)THEN
            IANF=4
            IDIR=1
          ELSE IF(INGRD.EQ.4)THEN
            IANF=7
            IDIR=1
          ELSE IF(INGRD.EQ.5)THEN
            IANF=13
            IDIR=1
          ELSE IF(INGRD.GE.6.AND.INGRD.LE.10)THEN
            IDIR=0
          ELSE
          GOTO 200
          END IF
          IF(IDIR.EQ.1)THEN
C +++++
C DIREKTE INTEGRATION NACH HAMMER
C +++++
            DINPKT(1)=DR1KOO(IANF+IINPKT,1)
            DINPKT(2)=DR1KOO(IANF+IINPKT,2)
            DINPKT(4)=DR1WIC(IANF+IINPKT  )
          ELSE
C +++++
C INTEGRATION NACH GAUSS-RADAU
C +++++
            IOR=INGRD-5
            IPOS=(IINPKT-1)/IOR+1
            JPOS=IINPKT-(IPOS-1)*IOR
            IANF=(IOR*(IOR-1))/2
            DINPKT(1)=DR2AI(IANF+IPOS)
                X1          =Z1-DINPKT(1)
            DINPKT(2)=DR2AJ(IANF+JPOS)*X1
            DINPKT(4)=DR2AS(IANF+IPOS)*DR2H(IANF+JPOS)*X1
          END IF
        ELSE IF(IELNUM.EQ.1)THEN
C +++++
C VIERECK
C +++++
C
C          NUMERISCHE    STUETZ-      INTEGRATIONS-
C INGRD     ORDNUNG      PUNKTZAHL       ART
C --------------------------------------------------
C  1         1 X  1       1 X 1         GAUSS
C  2         3 X  3       2 X 2         GAUSS
C  3         5 X  5       3 X 3         GAUSS
C  4         7 X  7       4 X 4         GAUSS
C  5         9 X  9       5 X 5         GAUSS
C  6        11 X 11       6 X 6         GAUSS
C  7        13 X 13       7 X 7         GAUSS
C  8        15 X 15       8 X 8         GAUSS
C  9        17 X 17       9 X 9         GAUSS
C
          IPOS=(IINPKT-1)/INGRD+1
          JPOS=IINPKT-(IPOS-1)*INGRD
          DINPKT(1)=VIEKOO(IPOS,INGRD)
          DINPKT(2)=VIEKOO(JPOS,INGRD)
          DINPKT(4)=VIEWIC(IPOS,INGRD)*VIEWIC(JPOS,INGRD)
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
C
C          NUMERISCHE       STUETZ-        INTEGRATIONS-
C INGRD     ORDNUNG         PUNKTZAHL         ART
C ------------------------------------------------------
C   1         1                  1           HAMMER
C   2         2                  4           HAMMER
C   3         3                  5           HAMMER
C
          IF(INGRD.EQ.1)THEN
            IANF=0
          ELSE IF(INGRD.EQ.2)THEN
            IANF=1
          ELSE IF(INGRD.EQ.3)THEN
            IANF=5
          ELSE
            GOTO 200
          END IF
          DINPKT(1)=TETKOO(IANF+IINPKT,1)
          DINPKT(2)=TETKOO(IANF+IINPKT,2)
          DINPKT(3)=TETKOO(IANF+IINPKT,3)
          DINPKT(4)=TETWIC(IANF+IINPKT  )
        ELSE IF(IELNUM.EQ.1)THEN
C +++++
C KEIL
C +++++
C
C          NUM.ORDN.     STUETZPKT.         INTEGRATIONSART
C INGRD   DRE.  VIE.    DRE.     VIER.       DREI.    VIER.
C ---------------------------------------------------------
C  1       1     1       1         1       DIREKT    GAUSS
C  2       2     3       3         2       DIREKT    GAUSS
C  3       2     3       3         2       DIREKT    GAUSS
C  4       4     3       6         2       DIREKT    GAUSS
C  5       5     5       7         3       DIREKT    GAUSS
C  6     1 X 1   1     1 X 1       1     GAUSS-RADAU GAUSS
C  7     3 X 3   3     2 X 2       2     GAUSS-RADAU GAUSS
C  8     5 X 5   5     3 X 3       3     GAUSS-RADAU GAUSS
C  9     7 X 7   7     4 X 4       4     GAUSS-RADAU GAUSS
C 10     9 X 9   9     5 X 5       5     GAUSS-RADAU GAUSS
C
C BEI DER DIREKTEN INTEGRATION NACH HAMMER MUSS FUER DIE EXAKTE
C INTEGRATION DES POLYNOMS XSI**I * ETA**J GELTEN:
C          (I+J)   KLEINER GLEICH   (NUMERISCHE ORDNUNG)
C
C
          IF(INGRD.EQ.1)THEN
            IANF = 0
            IORVI= 1
            IORDR= 1
            IDIR = 1
          ELSE IF(INGRD.EQ.2)THEN
            IANF = 1
            IORVI= 2
            IORDR= 3
            IDIR = 1
          ELSE IF(INGRD.EQ.3)THEN
            IANF = 4
            IORVI= 2
            IORDR= 3
            IDIR = 1
          ELSE IF(INGRD.EQ.4)THEN
            IANF = 7
            IORVI= 2
            IORDR= 6
            IDIR = 1
          ELSE IF(INGRD.EQ.5)THEN
            IANF =13
            IORVI= 3
            IORDR= 7
            IDIR = 1
          ELSE IF(INGRD.EQ.6)THEN
            IANF = 0
            IORVI= 1
            IORDR= 1
            IDIR = 0
          ELSE IF(INGRD.EQ.7)THEN
            IANF = 1
            IORVI= 2
            IORDR= 2
            IDIR = 0
          ELSE IF(INGRD.EQ.8)THEN
            IANF = 3
            IORVI= 3
            IORDR= 3
            IDIR = 0
          ELSE IF(INGRD.EQ.9)THEN
            IANF = 6
            IORVI= 4
            IORDR= 4
            IDIR = 0
          ELSE IF(INGRD.EQ.10)THEN
            IANF =10
            IORVI= 5
            IORDR= 5
            IDIR = 0
          ELSE
            GOTO 200
          END IF
          IF(IDIR.EQ.1)THEN
C +++++
C DIREKTE INTEGRATION NACH HAMMER PLUS GAUSS-QUADRATUR
C +++++
            KPOS=(IINPKT-1)/IORDR+1
            JINPKT=IINPKT-(KPOS-1)*IORDR
            DINPKT(1)=DR1KOO(IANF+JINPKT,1    )
            DINPKT(2)=DR1KOO(IANF+JINPKT,2    )
            DINPKT(3)=VIEKOO(KPOS       ,IORVI)
            DINPKT(4)=DR1WIC(IANF+JINPKT      )*VIEWIC(KPOS,IORVI)
          ELSE
C +++++
C INTEGRATION NACH GAUSS-RADAU PLUS GAUSS-QUADRATUR
C +++++
            IORDQ=IORDR*IORDR
            KPOS=(IINPKT-1)/IORDQ+1
            JINPKT=IINPKT-(KPOS-1)*IORDQ
            IPOS=(JINPKT-1)/IORDR+1
            JPOS=JINPKT-(IPOS-1)*IORDR
            DINPKT(1)=DR2AI(IANF+IPOS)
                X1   =Z1-DINPKT(1)
            DINPKT(2)=DR2AJ(IANF+JPOS)*X1
            DINPKT(3)=VIEKOO(KPOS,IORVI)
            DINPKT(4)=DR2AS(IANF+IPOS)*DR2H(IANF+JPOS)*X1
     *               *VIEWIC(KPOS,IORVI)
   20       CONTINUE
          END IF
        ELSE IF(IELNUM.EQ.2)THEN
C +++++
C QUADER
C +++++
C          NUMERISCHE       STUETZ-        INTEGRATIONS-
C INGRD     ORDNUNG         PUNKTZAHL         ART
C ------------------------------------------------------
C  1       1 X  1 X  1      1 X 1 X 1        GAUSS
C  2       3 X  3 X  3      2 X 2 X 2        GAUSS
C  3       5 X  5 X  5      3 X 3 X 3        GAUSS
C  4       7 X  7 X  7      4 X 4 X 4        GAUSS
C  5       9 X  9 X  9      5 X 5 X 5        GAUSS
C  6      11 X 11 X 11      6 X 6 X 6        GAUSS
C  7      13 X 13 X 13      7 X 7 X 7        GAUSS
C  8      15 X 15 X 15      8 X 8 X 8        GAUSS
C  9      17 X 17 X 17      9 X 9 X 9        GAUSS
C
          IF(INGRD.GT.9)THEN
C +++++
C INTEGRATIONSGRAD NICHT IMPLEMENTIERT, FEHLER
C +++++
            GOTO 200
          END IF
          INGRQ=INGRD*INGRD
          KPOS=(IINPKT-1)/INGRQ+1
          JINPKT=IINPKT-(KPOS-1)*INGRQ
          IPOS=(JINPKT-1)/INGRD+1
          JPOS=JINPKT-(IPOS-1)*INGRD
          DINPKT(1)=VIEKOO(IPOS,INGRD)
          DINPKT(2)=VIEKOO(JPOS,INGRD)
          DINPKT(3)=VIEKOO(KPOS,INGRD)
          DINPKT(4)=VIEWIC(IPOS,INGRD)*VIEWIC(JPOS,INGRD)
     *             *VIEWIC(KPOS,INGRD)
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
