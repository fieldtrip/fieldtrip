      SUBROUTINE FOFAGN(IELTYP,R,S,T,LKNEL ,XLK,FWE,FAG,DETJ,IERR)      NEU
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
C FOFAGN:
C =======
C
C BERECHNUNG DER GOLBALEN ABLEITUNGEN DER FORMFUNKTIONEN UND DER DETERMINANTE
C DER JAKOBISCHEN MATRIX.
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
C R         D   I                 1. LOKALE KOORDINATE DES AUSWERTEPUNKTES
C S         D   I                 2. LOKALE KOORDINATE DES AUSWERTEPUNKTES
C T         D   I                 3. LOKALE KOORDINATE DES AUSWERTEPUNKTES
C LKNEL     I   I  (MKNPEP)       KNOTEN-ELEMENT-ZUORDNUNG UNCOMPRESSED
C XLK       D   I  (3,MKNPEP)     KNOTENKOORDINATEN IN LOKALER NUMERIERUNG
C FWE       D   O  (MKNPEP)       WERTE DER FORMFUNKTIONEN
C FAG       D   O  (3,MKNPEP)     WERTE DER GLOBALEN ABL. DER FORMF.
C IERR      I   O                 FEHLERPARAMETER
C                                 0: ALLES KLAR
C                                 1: FALSCHER ELEMENTTYP
C                                 2: DETJ KLEINER ODER GLEICH NULL
C
C#
C
      INCLUDE 'fofulib.inc'
C
      CHARACTER*6 CHTYP
C
      DIMENSION LKNEL (  MKNPEP),XLK   (3,MKNPEP),FWE   (  MKNPEP),
     *          FAG   (3,MKNPEP),FAL   (3,MKNPEP)
C
C +++++
C ELEMENTTYP DEKODIEREN
C +++++
      WRITE(CHTYP,'(I5)')IELTYP
      READ(CHTYP,'(BN,2X,3I1)')IGLDIM,IELNUM,ILODIM
      IF(ILODIM.EQ.1)THEN
C +++++
C EINDIMENSIONALE ELEMENTE
C +++++
        IF(IELNUM.EQ.0)THEN
C +++++
C STAB
C +++++
          CALL FOFU01(1,1,R,LKNEL,FWE,FAL)
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
          CALL FOFU02(1,1,R,S,LKNEL,FWE,FAL)
        ELSE IF(IELNUM.EQ.1)THEN
C +++++
C VIERECK
C +++++
          CALL FOFU12(1,1,R,S,LKNEL,FWE,FAL)
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
          CALL FOFU03(1,1,R,S,T,LKNEL,FWE,FAL)
        ELSE IF(IELNUM.EQ.1)THEN
C +++++
C KEIL
C +++++
          CALL FOFU13(1,1,R,S,T,LKNEL,FWE,FAL)
        ELSE IF(IELNUM.EQ.2)THEN
C +++++
C QUADER
C +++++
          CALL FOFU23(1,1,R,S,T,LKNEL,FWE,FAL)
        ELSE
          GOTO 100
        END IF
      ELSE
        GOTO 100
      END IF
C +++++
C IKNMAX : HOECHSTE BESETZTE KNOTENNUMMER DES ELEMENTES
C +++++
      MKNANZ=MKNEL(IELNUM+1,ILODIM)
      IKNMAX=0
      DO 10 IKNANZ=1,MKNANZ
      IF(LKNEL(IKNANZ).NE.0)IKNMAX=IKNANZ
   10 CONTINUE
      IF(IGLDIM.EQ.ILODIM)THEN
C +++++
C
C NORMALE ELEMENTE
C ----------------
C
C +++++
C +++++
C BERECHNUNG DER JOKOBISCHEN MATRIX, DEREN DETERMINANTE UND DER
C GLOBALEN ABLEITUNGEN
C Z.ZT. SIND NUR NORMALE ELEMENTE MIT EINER QUADRATISCHEN JAKOBI-
C MATRIX ZUGELASSEN
C +++++
        IF(ILODIM.EQ.1)THEN
C +++++
C 1-DIMENSIONALE ELEMENTE
C +++++
          JAC11=Z0
          DO 20 K=1,IKNMAX
          IF(LKNEL(K).NE.0)THEN
            JAC11=JAC11+FAL(1,K)*XLK(1,K)
          END IF
   20     CONTINUE
          DETJ=JAC11
          IF(DETJ.LE.Z0)GOTO 200
          JACI11= Z1/DETJ
          DO 30 I=1,IKNMAX
          FAG(1,I)=JACI11*FAL(1,I)
   30     CONTINUE
        ELSE IF(ILODIM.EQ.2)THEN
C +++++
C 2-DIMENSIONALE ELEMENTE
C +++++
          JAC11=Z0
          JAC12=Z0
          JAC21=Z0
          JAC22=Z0
          DO 40 K=1,IKNMAX
          IF(LKNEL(K).NE.0)THEN
            JAC11=JAC11+FAL(1,K)*XLK(1,K)
            JAC12=JAC12+FAL(1,K)*XLK(2,K)
            JAC21=JAC21+FAL(2,K)*XLK(1,K)
            JAC22=JAC22+FAL(2,K)*XLK(2,K)
          END IF
   40     CONTINUE
          DETJ=JAC11*JAC22-JAC21*JAC12
          IF(DETJ.LE.Z0)GOTO 200
          JACI11= JAC22/DETJ
          JACI12=-JAC12/DETJ
          JACI21=-JAC21/DETJ
          JACI22= JAC11/DETJ
          DO 50 I=1,IKNMAX
          FAG(1,I)=JACI11*FAL(1,I)+JACI12*FAL(2,I)
          FAG(2,I)=JACI21*FAL(1,I)+JACI22*FAL(2,I)
   50     CONTINUE
        ELSE IF(ILODIM.EQ.3)THEN
C +++++
C 3-DIMENSIONALE ELEMENTE
C +++++
          JAC11=Z0
          JAC12=Z0
          JAC13=Z0
          JAC21=Z0
          JAC22=Z0
          JAC23=Z0
          JAC31=Z0
          JAC32=Z0
          JAC33=Z0
          DO 60 K=1,IKNMAX
          IF(LKNEL(K).NE.0)THEN
            JAC11=JAC11+FAL(1,K)*XLK(1,K)
            JAC12=JAC12+FAL(1,K)*XLK(2,K)
            JAC13=JAC13+FAL(1,K)*XLK(3,K)
            JAC21=JAC21+FAL(2,K)*XLK(1,K)
            JAC22=JAC22+FAL(2,K)*XLK(2,K)
            JAC23=JAC23+FAL(2,K)*XLK(3,K)
            JAC31=JAC31+FAL(3,K)*XLK(1,K)
            JAC32=JAC32+FAL(3,K)*XLK(2,K)
            JAC33=JAC33+FAL(3,K)*XLK(3,K)
          END IF
   60     CONTINUE
          ADJ11= JAC22*JAC33-JAC32*JAC23
          ADJ12=-JAC21*JAC33+JAC31*JAC23
          ADJ13= JAC21*JAC32-JAC31*JAC22
          ADJ21=-JAC12*JAC33+JAC32*JAC13
          ADJ22= JAC11*JAC33-JAC31*JAC13
          ADJ23=-JAC11*JAC32+JAC31*JAC12
          ADJ31= JAC12*JAC23-JAC22*JAC13
          ADJ32=-JAC11*JAC23+JAC21*JAC13
          ADJ33= JAC11*JAC22-JAC21*JAC12
          DETJ=JAC11*ADJ11+JAC12*ADJ12+JAC13*ADJ13
          IF(DETJ.LE.Z0)GOTO 200
          JACI11=ADJ11/DETJ
          JACI12=ADJ21/DETJ
          JACI13=ADJ31/DETJ
          JACI21=ADJ12/DETJ
          JACI22=ADJ22/DETJ
          JACI23=ADJ32/DETJ
          JACI31=ADJ13/DETJ
          JACI32=ADJ23/DETJ
          JACI33=ADJ33/DETJ
          DO 70 I=1,IKNMAX
          FAG(1,I)=JACI11*FAL(1,I)+JACI12*FAL(2,I)+JACI13*FAL(3,I)
          FAG(2,I)=JACI21*FAL(1,I)+JACI22*FAL(2,I)+JACI23*FAL(3,I)
          FAG(3,I)=JACI31*FAL(1,I)+JACI32*FAL(2,I)+JACI33*FAL(3,I)
   70     CONTINUE
        END IF
      ELSE IF(IGLDIM-ILODIM.EQ.1)THEN
C +++++
C
C OBERFLAECHEN-ELEMENTE
C
C +++++
        IF(IGLDIM.EQ.2)THEN
C +++++
C RANDLINIEN IN DER EBENE
C +++++
          JAC11=Z0
          JAC12=Z0
          DO 300 K=1,IKNMAX
          IF(LKNEL(K).NE.0)THEN
            JAC11=JAC11+FAL(1,K)*XLK(1,K)
            JAC12=JAC12+FAL(1,K)*XLK(2,K)
          END IF
  300     CONTINUE
c$$$          ADJ21=-JAC12
c$$$          ADJ22= JAC11
          ADJ21= JAC12
          ADJ22=-JAC11
          DETJ=SQRT(ADJ21*ADJ21+ADJ22*ADJ22)
          if(detj.le.z0) goto 200
c
c---neu (!!!)
c---globale ableitungen der formfunktionen (c) Adrian Rien"acker, Okt 1994
          uxilen=z1/(jac11*jac11+jac12*jac12)
          DO 280 I=1,IKNMAX
             FAG(1,I)=FAL(1,I)*jac11*uxilen
             FAG(2,I)=FAL(1,I)*jac12*uxilen
  280     CONTINUE
c---ende neu (!!!)
        ELSE IF(IGLDIM.EQ.3)THEN
C +++++
C FLAECHEN ALS OBERFLAECHE EINER 3D-STRUKTUR
C +++++
          JAC11=Z0
          JAC12=Z0
          JAC13=Z0
          JAC21=Z0
          JAC22=Z0
          JAC23=Z0
          DO 310 K=1,IKNMAX
          IF(LKNEL(K).NE.0)THEN
            JAC11=JAC11+FAL(1,K)*XLK(1,K)
            JAC12=JAC12+FAL(1,K)*XLK(2,K)
            JAC13=JAC13+FAL(1,K)*XLK(3,K)
            JAC21=JAC21+FAL(2,K)*XLK(1,K)
            JAC22=JAC22+FAL(2,K)*XLK(2,K)
            JAC23=JAC23+FAL(2,K)*XLK(3,K)
          END IF
  310     CONTINUE
          ADJ31= JAC12*JAC23-JAC22*JAC13
          ADJ32=-JAC11*JAC23+JAC21*JAC13
          ADJ33= JAC11*JAC22-JAC21*JAC12
          DETJ=SQRT(ADJ31*ADJ31+ADJ32*ADJ32+ADJ33*ADJ33)
          if(detj.le.z0) goto 200
c
c---neu (!!!) hier werden die globalen ableitungen entlang der Fl"ache,
c   projiziert auf die kartesischen
c   Richtungen berechnet!!! Adrian Rienaecker, Okt. 1994
          uxilen=z1/(jac11*jac11+jac12*jac12+jac13*jac13)
          vetlen=z1/(jac21*jac21+jac22*jac22+jac23*jac23)
          DO 270 I=1,IKNMAX
             FAG(1,I)=FAL(1,I)*jac11*uxilen+fal(2,i)*jac21*vetlen
             FAG(2,I)=FAL(1,I)*jac12*uxilen+fal(2,i)*jac22*vetlen
             FAG(3,I)=FAL(1,I)*jac13*uxilen+fal(2,i)*jac23*vetlen
  270     CONTINUE
c---ende neu (!!!)
        END IF
      else if (igldim-ilodim.eq.2) then
c
c---neu (!!!) stabelemente im 3D (c) Adrian Rienaecker, Okt. 1994
          JAC11=Z0
          JAC12=Z0
          JAC13=Z0
          DO 260 K=1,IKNMAX
          IF(LKNEL(K).NE.0)THEN
            JAC11=JAC11+FAL(1,K)*XLK(1,K)
            JAC12=JAC12+FAL(1,K)*XLK(2,K)
            JAC13=JAC13+FAL(1,K)*XLK(3,K)
          END IF
  260     CONTINUE
          DETJ=SQRT(jac11*jac11+jac12*jac12+jac13*jac13)
          if(detj.le.z0) goto 200
          uxilen=z1/(detj*detj)
          DO 290 I=1,IKNMAX
             FAG(1,I)=FAL(1,I)*jac11*uxilen
             FAG(2,I)=FAL(1,I)*jac12*uxilen
             FAG(3,I)=FAL(1,I)*jac13*uxilen
  290     CONTINUE
c
c--- ende neu (!!!)
      END IF
      IERR=0
      R E T U R N
  100 IERR=1
C +++++
C 100 FEHLER: FLASCHER ELEMENTTYP
C +++++
      R E T U R N
  200 IERR=2
C +++++
C 200 FEHLER: DETJ GLEINER ODER GLEICH NULL
C +++++
      R E T U R N
      END
