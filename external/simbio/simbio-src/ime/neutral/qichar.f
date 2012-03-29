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
      CHARACTER*80 FUNCTION QICHAR (TEXT,IRERAN,CDFWE,LTI,IERR)
C     Liest einen Character-String ein
C
C     ****************
C     * STANDARD F77 *
C     ****************
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C TEXT      C   I  *(*)           Auszugebender Text
C IRERAN    I   I                 Ausrichtung des Textes
C                                     0 : linksbuendig
C                                 sonst : rechtsbuendig in Spalte IRERAN
C                                 mit: IREMIN < IRERAN < IREMAX
C CDFWE     C   I  *(*)           Defaultwert der Antwort
C                                 Wird CHAR(0) uebergeben, so wird
C                                 kein Defaultwert ausgegeben, aber es muss
C                                 eine Eingabe erfolgen (KEIN Blank), bei
C                                 CHAR(1) wird kein Doppelpunkt und kein 
C                                 Defaultwert ausgegeben. Es muss keine
C                                 Eingabe erfolgen (nur RETURN reicht).
C LTI       I   O                 Laenge des eingelesenen String QICHAR in Bu.
C IERR      I   O                 Fehlercode  0=kein Fehler
C                                 1=Fehler; 2=Benutzerabbruch der Eingabe
C                                 -1=Formatueberschreitung (Warnung)
      INCLUDE 'imelibh'
C
      INTEGER QCLEN
C
      CHARACTER*60  ZEILE
      CHARACTER*80  BUFFER
      CHARACTER*(*) TEXT, CDFWE
      INTEGER LETE, IRERAN, IREMIN, IREMAX, LDFWE, LTI,
     $     I, J, NBLK, NGESAM, IERR, LDF, LE, ICLIPT, ICLIPD
C
      IERR=0
      IREMIN = 10
      IREMAX = 60
      ICLIPT = 0
      ICLIPD = 0
      NBLK = 0
      ZEILE=' '
C
C   Abfrage, ob Default-Text vorhanden oder nicht,
C   und ggf. Bestimmen der Laenge
C
      IF (CDFWE .EQ. CHAR(0)) THEN
         LDFWE=0
      ELSEIF (CDFWE .EQ. CHAR(1)) THEN
         LDFWE=0
      ELSE
         LDFWE = QCLEN(CDFWE)
      ENDIF
      LDF = LDFWE
C
C     Bestimmen der Laenge des auszugebenden Textes
C
      LETE = QCLEN (TEXT)
      LE = LETE
C
      IF (CDFWE .EQ. CHAR(0)) THEN
C     ohne Default (mit Blank und Doppelpunkt)
         NGESAM = LETE + 2
      ELSEIF (CDFWE .EQ. CHAR(1)) THEN
C     ohne Default (ohne Blank und Doppelpunkt)
         NGESAM = LETE
      ELSE
C     mit Default (mit Blank, 2 mal Klammer und Doppelpunkt)
         NGESAM = LETE + 2 + LDFWE + 2
      ENDIF
C
C     Berechne die Anzahl der Blanks, die als erstes in der auszugebenden
C     Zeile ZEILE stehen muessen, damit der Text TEXT
C     evtl. gefolgt von der Angabe des Defaultwertes
C     rechtsbuendig in der Spalte IRERAN endet.
C
C   Kontrolle von IRERAN
      IF (IRERAN .NE. 0) THEN
         IF (IRERAN .LT. IREMIN) IERR = -1
         IF (IRERAN .GT. IREMAX) IERR = -1
      ENDIF
C
C  1. Text nicht ausrichten:
      IF (IRERAN.EQ.0) THEN
         IF (NGESAM.GT.IREMAX) THEN
            IERR = -1
C     Klippen
C     LE, LDF   Laenge des ggf. abgeschnittenen Textes
C               bzw. Default-Textes
C
C     Zuerst Default-Text klippen, Platz fuer 3 Punkte reservieren
            IF (NGESAM-LDFWE+3.LE.IREMAX) THEN
               LDF = IREMAX -LETE -2 -3
               ICLIPD = 1
C     Text auch klippen, weil Default allein nicht reicht
            ELSE
               LDF = MIN (3, LDFWE)
               LE  = IREMAX -LDF -2 -3
               ICLIPD = 1
               ICLIPT = 1
            ENDIF
         ENDIF
C
C  2. Text ausrichten, Berechnung der Anzahl der fuehrenden Blanks
      ELSE
         IF (NGESAM.GT.IRERAN) THEN
            IERR = -1
C
C     Klippen
C     zuerst Default-Text klippen, Platz fuer 3 Punkte reservieren
            IF (NGESAM-LDFWE+1+3.LE.IRERAN) THEN
               LDF = IRERAN -LETE+1 -2 -3
               ICLIPD = 1
C     Text auch klippen, weil Default allein nicht reicht
            ELSE
               LDF = MIN (3, LDFWE)
               LE  = IRERAN -LDF -2 -3
               if (ldf.ne.ldfwe) ICLIPD = 1
               ICLIPT = 1
            ENDIF
         ENDIF

         IF (CDFWE .EQ. CHAR(0)) THEN
C     ohne Default (mit Blank und Doppelpunkt)
            NGESAM = LE + 2
         ELSEIF (CDFWE .EQ. CHAR(1)) THEN
C     ohne Default (ohne Blank und Doppelpunkt)
            NGESAM = LE
         ELSE
C     mit Default (mit Blank, 2 mal Klammer und Doppelpunkt)
            NGESAM = LE + 2 + LDF + 2
         ENDIF


C
C     Berechnung der fuehrenden Blanks
         NBLK = IRERAN - NGESAM
C
      ENDIF
C
C     Ausgabezeile zusammenstellen
C
      IF (IERR .EQ. 1) THEN
C
C     Fehler, gebe Fehlermeldung aus:
         ZEILE=' >> Fehler in QICHAR <<'
         RETURN
      ELSE
C
C     Korrekte Ausgabezeile moeglich:
C
         I=NBLK+1
         J=I+LE-1
         ZEILE(I:J)=TEXT(1:LE)
         IF (ICLIPT.EQ.1) ZEILE(J-2:J)= '...'
         IF (CDFWE .EQ. CHAR(1)) THEN
            ZEILE(J+1:J+1)=' '
         ELSEIF (CDFWE. EQ. CHAR(0)) THEN
            ZEILE(J+1:J+2)=' :'
         ELSE
C     
C     mit Default
C     
            I=J+1
            J=I+1
            ZEILE(I:J)=' ['
            I=J+1
            J=I+LDF-1
            ZEILE(I:J)=CDFWE(1:LDF)
            IF (ICLIPD.EQ.1) ZEILE(J-2:J)='...'
            ZEILE(J+1:J+2)=']:'
         END IF
      END IF

c Kommentar in den Protokollfile schreiben
         
      call qfhist(zeile,1)
      if (ntprot.ne.0) write(ntprot,'(''#'',a)') zeile(1:qclen(zeile))
      
 5    CONTINUE
      IF(TAST.AND.TERM)THEN

c Lesen von Tastatur und Schreiben aufs Terminal

 10      CONTINUE
         BUFFER=' '
         CALL QDTEXT(ZEILE, 1)
         READ(NTTAST,12,ERR=10,END=8020) BUFFER
 12      FORMAT (A)
      ELSE
         BUFFER=' '
         READ(NTTAST,12,ERR=8010,END=8020)BUFFER
         ihash=index(buffer,'#')
         if (ihash.eq.1) goto 5
         if (ihash.gt.0) buffer(ihash:)=' '
         CALL QDTEXT(ZEILE(1:qclen(zeile))//buffer, 0)
      END IF

 20   CONTINUE
      LTI = QCLEN(BUFFER)
      call qfhist(BUFFER(1:MAX(1,LTI)),2)
      IF (NTPROT.NE.0) WRITE(NTPROT,'(A)') BUFFER(1:MAX(1,LTI))
      
      txtalt=buffer
      call qvfind(buffer)

      IF(LTI.EQ.0.AND.CDFWE.EQ.CHAR(0).AND.LETE.GE.0) GOTO 5
      IF(LTI.EQ.0.AND.CDFWE.EQ.CHAR(1).AND.LETE.GE.0) THEN
C     Blank zurueckgeben
         QICHAR=' '
         RETURN
      ENDIF
      IF(LTI.EQ.0.AND.LDFWE.NE.0)THEN
C
C     Setze den Default-Wert ein
C
         LTI=LDFWE
         BUFFER=' '
         BUFFER(1:LTI)=CDFWE(1:LTI)
         txtalt=buffer
         call qvfind(buffer)
      END IF
      IF(.NOT.TERM)THEN
         WRITE(NTTERM,9020)ZEILE(1:NGESAM),BUFFER(1:LTI)
 9020 FORMAT(1X,A,A)
      END IF
C
C     Kopieren des Textes vom Buffer in die Function-Variable
      QICHAR=BUFFER
C
C     Wenn Abbruch eingegeben:
      IF (BUFFER(1:LABBRU) .EQ. CABBRU(1:LABBRU)) IERR=2
C
      RETURN
C
C     Bei Fehlern im Batchbetrieb Programmstop:
 8010 STOP 'Fehler beim Lesen in QICHAR'
 8020 CONTINUE
      ierr=2
      return
      END

