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
      PROGRAM LBTEST
C     -----------------------------------------------
C     Testet die maschinenabhaengige Library des IME.
C     -----------------------------------------------
C
C     Erweiterte Version von TESLIB
C
C      Das Programm ruft menuegesteuert drei Unterprogramme auf:
C      a) System-Routinen
C      b) Dialog-Routinen
C      c) File-Management-Routinen
C
      CHARACTER*80 ZEILE
      CHARACTER*20 CABBRU
      CHARACTER*75 PKTMEN(15)
      CHARACTER*2 ANTMEN, MENUE2
      CHARACTER*1 WEITER
      CHARACTER*16 QLBVER
      INTEGER  NTTAST, NTTERM, NTERRO, LABBRU, NERR
C
      WRITE (*, 9020)
 9020 FORMAT(5(/),5X,
     $'   ---------------------------------------------------    '/5X,
     $'         Programm zum Testen der IME-Library              '/5X,
     $'   ---------------------------------------------------    '/5X,
     $'                                                          '/5X,
     $'   Ein Programm des     INSTITUTS FUER MASCHINENELEMENTE  '/5X,
     $'                        UND MASCHINENGESTALTUNG           '/5X,
     $'                        DER RWTH AACHEN                   '/5X,
     $'                                                          '/5X,
     $'  in Zusammenarbeit mit FORSCHUNGSVEREINIGUNG             '/5X,
     $'                        ANTRIEBSTECHNIK E.V.  (FVA)       '/5X,
     $'                        Lyoner Strasse 18                 '/5X,
     $'                        6000 Frankfurt/Main               ')
C
C
C     Zuallererst: Initialisieren der Library
C
      LABBRU = 0
      CALL QLBINI (NTTAST, NTTERM, NTERRO, CABBRU, LABBRU)
C
      WRITE  (*,9021)
 9021 FORMAT (/,1X,'(QLBINI:) Library initialisiert mit:',/)
      WRITE  (*,9022) NTTAST, NTTERM, NTERRO, CABBRU, LABBRU,
     $     QLBVER()
 9022 FORMAT (1X,'NTTAST (Standard-Eingabe-Unit)          ',I2,/,
     $        1X,'NTTERM (Standard-Ausgabe-Unit)          ',I2,/,
     $        1X,'NTERRO (Standard-Fehler-Unit)           ',I2,/,
     $        1X,'CABBRU (Programm-Abbruch moeglich mit:)  ',A, /,
     $        1X,'LABBRU (Laenge des Abbruch-Strings)     ',I2,/,
     $        1X,'Versionsnummer                           ',A)
C
      WRITE (*,'(A)') ' Weiter mit RETURN...'
      READ  (*,4711) WEITER
 4711 FORMAT (A1)
C
C      Erzeugen des Hauptmenues
C
 100  PKTMEN(1)='HAUPTMENUE des Testprogramms'
      PKTMEN(2)='SY...System-Routinen'
      PKTMEN(3)='DI...Dialog-Routinen'
      PKTMEN(4)='FM...File-Management'
      PKTMEN(5)='EN...ENDE'
      ANTMEN=MENUE2('EN', 4, PKTMEN)
C
C     Verzweigen in die Menues
C
      IF (ANTMEN.EQ.'SY') CALL SYTEST (NTTAST, NTTERM, CABBRU, LABBRU)
      IF (ANTMEN.EQ.'DI') CALL DITEST (NTTAST, NTTERM)
      IF (ANTMEN.EQ.'FM') CALL FMTEST (NTTAST, NTTERM)
      IF (ANTMEN.EQ.'EN') STOP 'ENDE Testprogramm'
      GOTO 100
C
      END
C
C     ############################################################
C
      SUBROUTINE SYTEST (NTTAST, NTTERM, CABBRU, LABBRU)
C
C     Testen der System-Routinen
C
C     QSUSER, QSDATI, QSCPUT, QSECHT, QSBZZZ
C     QSPREA, QSPDBL
C
      CHARACTER*8  DAT, UHRZ
      CHARACTER*80 LOGNA, FILNAM
      CHARACTER*20 CABBRU
      CHARACTER*1 WEITER
      INTEGER  LELO, NTTAST, NTTERM, LABBRU, NERR
      REAL QSPREA, REAMIN, REAMAX, REAGEN
      DOUBLE PRECISION  DUMMY, T1, T2, QSCPUT, QSECHT
      DOUBLE PRECISION  QSPDBL, DBLMIN, DBLMAX, DBLGEN
C
      CALL LINIE0
      CALL LINIE2
      WRITE(NTTERM, 8000)
 8000 FORMAT(25X,' System-Informationen')
      CALL LINIE1
      CALL LINIE2
C
      WRITE (NTTERM,'(A)') ' Test von QSUSER'
      CALL QSUSER (LOGNA, LELO)
      WRITE (NTTERM, 8001) LOGNA, LELO
 8001 FORMAT(1X,'Benutzername  : ',A,/,
     *       1X,'mit der Laenge: ',I2,/)
C
      WRITE (NTTERM,'(A)') ' Test von QSDATI'
      CALL QSDATI (DAT, UHRZ)
      WRITE (NTTERM, 8002) DAT, UHRZ
 8002 FORMAT(1X,'Datum  : ',A,/,
     *       1X,'Uhrzeit: ',A,/)
C
      WRITE (NTTERM,'(A)') ' Test von QSCPUT und QSECHT'
      DUMMY = QSCPUT (1, 0, NERR)
      IF (NERR.NE.0) THEN
      WRITE (NTTERM, 8910)
 8910 FORMAT(1X,/,
     * 1X,'SUBROUTINE QSCPUT liefert einen Fehler',/,
     * 1X,'beim Starten der Stoppuhr.')
      STOP
      ENDIF
      DUMMY = QSECHT (1, 0, NERR)
      IF (NERR.NE.0) THEN
      WRITE (NTTERM, 8920)
 8920 FORMAT(1X,/,
     * 1X,'SUBROUTINE QSECHT liefert einen Fehler',/,
     * 1X,'beim Starten der Stoppuhr.')
      STOP
      ENDIF
      Y = 0.0D0
      DO 10 I=1, 50000
         Y = Y+SIN(DBLE(I))
   10 CONTINUE
      T1 = QSCPUT (1, 1, NERR)
      IF (NERR.NE.0) THEN
      WRITE (NTTERM, 8911)
 8911 FORMAT(1X,/,
     * 1X,'SUBROUTINE QSCPUT liefert einen Fehler',/,
     * 1X,'beim Anhalten der Stoppuhr.')
      STOP
      ENDIF
      T2 = QSECHT (1, 1, NERR)
      IF (NERR.NE.0) THEN
      WRITE (NTTERM, 8921)
 8921 FORMAT(1X,/,
     * 1X,'SUBROUTINE QSECHT liefert einen Fehler',/,
     * 1X,'beim Anhalten der Stoppuhr.')
      STOP
      ENDIF
C
      WRITE (NTTERM, 8004) T1, T2
 8004 FORMAT(
     * 1X,'mit 50000 mal Sinus:',/,
     * 1X,'  CPU-Zeit: ',D14.7,/,
     * 1X,'Echte Zeit: ',D14.7)
C
 4711 FORMAT (A1)
C
C  2. Lauf mit der Stoppuhr
C
      WRITE (NTTERM,*)
      WRITE (NTTERM,'(A)') ' Test von QSBZZZ'
      WRITE (NTTERM, 8030)
 8030 FORMAT(1X,'mit einer Pause von 5 Sekunden:')
      DUMMY = QSCPUT (1, 0, NERR)
      IF (NERR.NE.0) THEN
      WRITE (NTTERM, 8912)
 8912 FORMAT(1X,/,
     * 1X,'SUBROUTINE QSCPUT liefert einen Fehler',/,
     * 1X,'beim 2. Starten der Stoppuhr.')
      STOP
      ENDIF
      DUMMY = QSECHT (1, 0, NERR)
      IF (NERR.NE.0) THEN
      WRITE (NTTERM, 8922)
 8922 FORMAT(1X,/,
     * 1X,'SUBROUTINE QSECHT liefert einen Fehler',/,
     * 1X,'beim 2. Starten der Stoppuhr.')
      STOP
      ENDIF
C
      CALL QSBZZZ (5)
      T1 = QSCPUT (1, 1, NERR)
      IF (NERR.NE.0) THEN
      WRITE (NTTERM, 8913)
 8913 FORMAT(1X,/,
     * 1X,'SUBROUTINE QSCPUT liefert einen Fehler',/,
     * 1X,'beim 2. Anhalten der Stoppuhr.')
      STOP
      ENDIF
      T2 = QSECHT (1,1,NERR)
      IF (NERR.NE.0) THEN
      WRITE (NTTERM, 8923)
 8923 FORMAT(1X,/,
     * 1X,'SUBROUTINE QSECHT liefert einen Fehler',/,
     * 1X,'beim 2. Anhalten der Stoppuhr.')
      STOP
      ENDIF
C
      WRITE (NTTERM, 8003) T1, T2
 8003 FORMAT(1X,'  CPU-Zeit: ',D14.7,/,
     *       1X,'Echte Zeit: ',D14.7)
C
      WRITE (NTTERM,*)
      WRITE (NTTERM,'(A)') ' Weiter mit RETURN...'
      READ  (NTTAST,4711) WEITER
C
C      Maschinengenauigkeit
C      QSPDBL und QSPREA
C
      WRITE (NTTERM,'(A)') ' Test von QSPDBL und QSPREA'
      WRITE (NTTERM,*)
 2    DBLMIN = QSPDBL(2)
      DBLMAX = QSPDBL(3)
      DBLGEN = QSPDBL(1)
      WRITE (NTTERM,8151) DBLMIN, DBLMAX, DBLGEN
 8151 FORMAT (1X,'Kleinste DOUBLE-Zahl :',D20.12,/,
     *        1X,'groesste DOUBLE-Zahl :',D20.12,/,
     *        1X,'Genauigkeit          :',D20.12,/)
C
      REAMIN = QSPREA(2)
      REAMAX = QSPREA(3)
      REAGEN = QSPREA(1)
      WRITE (NTTERM,8152) REAMIN, REAMAX, REAGEN
 8152 FORMAT (1X,'Kleinste REAL-Zahl :',E15.7,/,
     *        1X,'groesste REAL-Zahl :',E15.7,/,
     *        1X,'Genauigkeit        :',E15.7,/)
C
      WRITE(NTTERM, '(A)') ' Test von QSCLR :'
      WRITE (NTTERM,*)
      WRITE (NTTERM,'(A)') ' Nach RETURN erfolgt ein Loeschen '//
     $     'des Bildschirms'
      READ  (NTTAST,4711) WEITER
      CALL QSCLR
C
      WRITE (NTTERM,'(A)') ' Weiter mit RETURN...'
      READ  (NTTAST,4711) WEITER
C
C     Systemshell aufrufen:
      WRITE (NTTERM,'(A)') ' Test von QSHELL'
      WRITE (NTTERM,*)
      CALL QSHELL(NERR)
      IF (NERR .GT. 0) THEN
         WRITE (NTTERM, '(A)') ' Fehler in QSHELL'
         STOP 'Fehler'
      ELSEIF (NERR .EQ. -2) THEN
         WRITE (NTTERM, '(A)') ' QSHELL auf diesem Rechner nicht '//
     $        'verfuegbar.'
         WRITE (NTTERM,'(A)') ' Weiter mit RETURN...'
         READ  (NTTAST,4711) WEITER
      ENDIF
C
C     Druckroutine testen:
      WRITE(NTTERM, '(A)') ' Test von QSPRNT'
      WRITE(NTTERM, '(A)') ' Geben Sie den Namen der Datei ein, die'//
     $     ' gedruckt werden soll:'
      READ (NTTAST, 8200) FILNAM
 8200 FORMAT (A)
      IF (FILNAM(1:LABBRU) .NE. CABBRU(1:LABBRU)) THEN
         CALL QSPRNT(FILNAM, NERR)
         IF (NERR .GT. 0) THEN
            WRITE (NTTERM, '(A)') ' Fehler in QSPRNT'
            STOP 'Fehler'
         ELSEIF (NERR .EQ. -2) THEN
            WRITE (NTTERM, '(A)') ' QSPRNT auf diesem Rechner nicht '//
     $           'verfuegbar.'
         ELSEIF (NERR .EQ. -1) THEN
            WRITE (NTTERM, '(A)') ' Warnung: Datei nicht gefunden.'
         ENDIF
         WRITE (NTTERM,'(A)') ' Weiter mit RETURN...'
         READ  (NTTAST,4711) WEITER
      ENDIF
      RETURN
      END
C
C     ############################################################
C
      SUBROUTINE DITEST (NTTAST, NTTERM)
C
C     Testen der Dialog-Routinen
C
C      QICHAR, QCTOUP,QCTOLO,
C      QINTE, QIREAL, QIDBL, QCITOA, QCRTOA, QCDTOA,
C      QCLEN, QIYENO, QDTEXT
C
      LOGICAL QIYENO
      CHARACTER*10  TE
      CHARACTER*80  ZEILE, QICHAR
      CHARACTER*20  CABBRU
      INTEGER  LE, I, NN, NE, LZEILE, NTTAST, NTTERM,
     $     LABBRU, NERR, QINTE, QCLEN
      REAL R, R0, QIREAL
      PARAMETER (R0=0.0)
      DOUBLE PRECISION  D, D0, QIDBL
      PARAMETER (D0=0.0D0)
C
      NN=0
      NE=0
C
      CALL LINIE0
      CALL LINIE2
      WRITE (NTTERM, 8005)
 8005 FORMAT(30X,'Test der Dialog-Routinen')
      CALL LINIE1
C
      CALL LINIE2
      CALL QDTEXT ('Test von QDTEXT :',1)
      CALL QDTEXT ('Linksbuendiger Text.',0)
      CALL QDTEXT ('Dieser Text sollte daneben stehen.',1)
      CALL QDTEXT ('Dieser Text sollte darunter stehen.',0)
C
      WRITE(NTTERM, *)
      WRITE(NTTERM, '(A)') ' Test von QCLEN. '
      WRITE(NTTERM, '(A)') ' Bitte geben Sie einen Text ein :'
      READ (NTTAST,'(A)')  ZEILE
      LZEILE = QCLEN (ZEILE)
      IF (LZEILE .EQ. 0) LZEILE=1
      WRITE (NTTERM, 6201) ZEILE(1:LZEILE), LZEILE
 6201 FORMAT (1X,'Der Text >',A,'<',/,
     $        1X,'hat die Laenge ',I2,/)
C
      WRITE(NTTERM, '(A)') ' Test von QICHAR :'
      ZEILE=QICHAR ('Bitte einen Text eingeben', 0, CHAR(0),
     $     LZEILE, NERR)
      IF (NERR .EQ. -1) THEN
         WRITE (NTTERM,'(A)')
     $        ' Warnung: Formatueberschreitung in der Routine QICHAR'
      ENDIF
      IF (NERR .EQ. 1) THEN
         WRITE (NTTERM,'(A)') ' FEHLER in der Routine QICHAR'
         STOP
      ENDIF
      IF (NERR .EQ. 2) THEN
         WRITE(NTTERM, 8020)
 8020 FORMAT(1X,'Abbruch der Eingabe durch Benutzer!')
      RETURN
      ELSE
         WRITE (NTTERM, 8006) ZEILE (1:LZEILE), LZEILE
 8006 FORMAT(1X,'Gelesener Text: ',A,/,
     *       1X,'Mit der Laenge: ',I3,/)
      ENDIF
C
      ZEILE=QICHAR ('Bitte nochmal einen Text eingeben', 0,
     $     'Dieser Defaultwert ist eindeutig zu lang',
     $     LZEILE, NERR)
      IF (NERR .EQ. -1) THEN
         WRITE (NTTERM,'(A)')
     $        ' Warnung: Formatueberschreitung in der Routine QICHAR'
      ENDIF
      IF (NERR .EQ. 1) THEN
         WRITE (NTTERM,'(A)') 'FEHLER in der Routine QICHAR'
         STOP
      ENDIF
      IF (NERR .EQ. 2) THEN
         WRITE(NTTERM, 8020)
      RETURN
      ELSE
         WRITE (NTTERM, 8006) ZEILE (1:LZEILE), LZEILE
      ENDIF
C
      WRITE (NTTERM,'(A)') ' Test von QCTOUP und QCTOLO:'
      CALL QCTOUP (ZEILE)
      WRITE (NTTERM, 8016) ZEILE (1:LEN(ZEILE))
 8016 FORMAT(1X, 'Text in Grossbuchstaben: ', A)
C
      CALL QCTOLO (ZEILE)
      WRITE (NTTERM, 8026) ZEILE (1:LEN(ZEILE))
 8026 FORMAT(1X, 'Text in Kleinbuchstaben: ', A, /)
C
      ZEILE=QICHAR('Druecken Sie nur RETURN', 0, CHAR(1), LZEILE, NERR)
C
      WRITE (NTTERM,*)
      WRITE (NTTERM, '(A)') ' Test von QINTE :'
      I=QINTE ('Bitte eine Integerzahl eingeben', 0, 0, 'NO',
     $     0, 0, 0, NERR)
      IF (NERR .EQ. -1) THEN
         WRITE (NTTERM,'(A)')
     $        ' Warnung: Formatueberschreitung in der Routine QINTE'
      ENDIF
      IF (NERR .EQ. 1) THEN
         WRITE (NTTERM,'(A)') ' FEHLER in der Routine QINTE'
         STOP
      ENDIF
      IF (NERR .EQ. 2) THEN
         WRITE(NTTERM, 8020)
         RETURN
         ENDIF
      WRITE (NTTERM, 8007) I
 8007 FORMAT(1X,'Gelesene Zahl: ', I10, /)
C
      WRITE(NTTERM, '(A)') ' Test von QCITOA :'
      CALL QCITOA (I, TE, LE)
      WRITE (NTTERM, 8017) TE
 8017 FORMAT(1X,'Zahl in CHARACTER-Darstellung: ',A,/)
C
      WRITE(NTTERM, '(A)') ' Test von QIREAL :'
      R=QIREAL ('Bitte eine REAL-Zahl eingeben',  0, R0, 'NO',
     $     0, R0, R0, NERR)
      IF (NERR .EQ. -1) THEN
         WRITE (NTTERM,'(A)')
     $        ' Warnung: Formatueberschreitung in der Routine QIREAL'
      ENDIF
      IF (NERR .EQ. 2) THEN
         WRITE(NTTERM, 8020)
         RETURN
         ENDIF
      WRITE (NTTERM, 8108) R
 8108 FORMAT(1X,'Gelesene Zahl : ',F12.5, /)
C
      WRITE(NTTERM, '(A)') ' Test von QCRTOA :'
      CALL QCRTOA (R, NN, NE, TE, LE)
      WRITE (NTTERM, 8018) TE(1:LE)
 8018 FORMAT(1X,'Zahl in CHARACTER-Darstellung: ',A,/)
C
      WRITE(NTTERM, '(A)') ' Test von QIDBL :'
      D = QIDBL ('Bitte eine Doublezahl eingeben', 0, Z0, 'NO',
     $     0, Z0, Z0, NERR)
      IF (NERR .EQ. -1) THEN
         WRITE (NTTERM,'(A)')
     $        ' Warnung: Formatueberschreitung in der Routine QIDBL'
      ENDIF
      IF (NERR .EQ. 2) THEN
         WRITE(NTTERM, 8020)
         RETURN
         ENDIF
      WRITE (NTTERM, 8008) D
 8008 FORMAT(1X,'Gelesene Zahl: ',D12.5, /)
C
      WRITE(NTTERM, '(A)') ' Test von QCDTOA :'
      CALL QCDTOA (D, NN, NE, TE, LE)
      WRITE (NTTERM, 8019) TE(1:LE)
 8019 FORMAT(1X, 'Zahl in CHARACTER-Darstellung: ',A,/)
C
      WRITE(NTTERM, '(A)') ' Test von QIYENO :'
 1111 IF (QIYENO ('Gefaellt Ihnen dieses Programm? (J/N)',
     $             0, 'N')) THEN
         WRITE (NTTERM,'(A)') ' Diese Antwort gefaellt uns.'
      ELSE
         WRITE (NTTERM,'(A)') ' Die Antwort NEIN lassen wir '//
     $        'nicht gelten !'
         GOTO 1111
      ENDIF
      WRITE (NTTERM,*)
      WRITE (NTTERM,'(A)') ' Ende Test Dialog-Routinen'
      WRITE (NTTERM,'(A)') ' Weiter mit RETURN...'
      READ  (NTTAST,4711) WEITER
 4711 FORMAT (A1)
C
      RETURN
      END
C
C     ############################################################
C
      SUBROUTINE FMTEST (NTTAST, NTTERM)
C
C     Testen von File-Management-Funktionen
C
C      QFILDA, QFREEU, QFRECL, QFSIZE, QFLEGN, QFLEGP,
C      QSDEFI, QSDEFO,
C      QFILAO, QFOPSE, QFOPDI, QFCRSE, QFCRDI,
C      QFDEL,  QFCLOS, QFCLAL,
C      QFWRDI, QFREDI
C
      CHARACTER*80  ZEILE, FILNA, FILNA1, FILNA2
      CHARACTER*20  CABBRU
      CHARACTER*80  CZAHL
      LOGICAL  TAST, TERM, YENO, QIYENO
      INTEGER NDWEGP, NRWEGP, NDWEG, NRWEG, L, LE,
     $     I, LZEILE, NT, J, LR, IRWEG, IDWEG, NTTAST, NTTERM,
     $     LABBRU, NERR
      PARAMETER ( NDWEGP=100, NRWEGP=100)
      DOUBLE PRECISION DWEG(NDWEGP), DSOLL, DUMMY
C
C      Deklaration der benutzten Library-Funktionen
C
      CHARACTER*80   QICHAR
      CHARACTER*1    WEITER
      LOGICAL QFILDA, QFLEGP, QFLEGN
      INTEGER IRELE, QCLEN, QFREEU, QFRECL
C
C
      NDWEG=NDWEGP
      NRWEG=NRWEGP
C
      WRITE (NTTERM, 8009)
 8009 FORMAT(1X,//,
     * 1X,' Datei-Routinen',/,
     * 1X,'----------------',/)
 8000 ZEILE=QICHAR ('Name einer NICHT-existierenden Datei',
     $     40, CHAR(0), L, NERR)
      IF (NERR .EQ. -1) THEN
         WRITE (NTTERM,'(A)')
     $        ' Warnung: Formatueberschreitung in der Routine QICHAR'
      ENDIF
      IF (NERR .EQ. 1) STOP 'Fehler in QICHAR'
      IF (NERR .EQ. 2) THEN
         WRITE(NTTERM, 8020)
 8020    FORMAT(1X,'Abbruch der Eingabe durch Benutzer !')
         RETURN
      ENDIF
C
C   Untersuchung von Dateinamen und Pfad
C
      IF (QFLEGN(ZEILE)) THEN
         WRITE (NTTERM, 8050) ZEILE(1:L)
 8050    FORMAT (1X,'QFLEGN: Dateiname >',A,'< zulaessig')
         ELSE
         WRITE (NTTERM, 8051) ZEILE(1:L)
 8051    FORMAT (1X,'QFLEGN: Dateiname >',A,'< unzulaessig')
         GO TO 8000
      ENDIF
C     QFLEGP (Pfad) kommt spaeter...
      IF (QFILDA (ZEILE)) THEN
         WRITE (NTTERM, 8904) ZEILE(1:L)
 8904    FORMAT(1X,/,
     *        1X,'QFILDA liefert den falschen Wert .TRUE.',/,
     *        1X,' Datei >',A,'< bereits vorhanden.')
         STOP 'Fehler'
      ELSE
         WRITE (NTTERM,6004) ZEILE(1:L)
 6004    FORMAT (1X,/,1X,'Bestaetigung von QFILDA:',/,
     $     1X,'Datei >',A,'< nicht vorhanden, Programm macht weiter',/)
      ENDIF
      NT=QFREEU()
      WRITE (NTTERM,6005) NT
 6005 FORMAT (1X,'Von QFREEU gewaehlte Unit :',I3,/)
C
      CALL QFOPSE (NT, ZEILE, 'NE', 'FO', NERR)
      IF (NERR .NE. 0) THEN
         WRITE (NTTERM, 8905)
 8905    FORMAT(1X,/,
     *        1X,'QFOPSE liefert einen Fehler.')
         STOP 'Fehler'
      ELSE
         WRITE (NTTERM, 6905)
 6905    FORMAT (1X,'QFOPSE fehlerfrei.',/)
      ENDIF
      IF ( QFILDA (ZEILE)) THEN
         WRITE (NTTERM, '(A)') ' QFILDA liefert wie erwartet .TRUE.'
      ELSE
         WRITE (NTTERM, '(A)') ' QFILDA liefert den falschen '//
     $        'Wert .FALSE.'
         STOP 'Fehler in QFILDA'
      ENDIF
C
C     Beschreiben der Datei
      DO 20 I=1, 20
         WRITE (NT, 8010, ERR=9007) I
 8010    FORMAT (I6)
 20   CONTINUE
C
      REWIND NT
C
C     Lesen der Datei
      DO 30 I=1, 20
         READ (NT, 8010, ERR=9008, END=9008) J
         IF (I.NE.J) THEN
            WRITE (NTTERM, 8909) I, J
 8909 FORMAT(1X,/,
     * 1X,'Vergleich beim formatierten sequentiellen Lesen falsch.',/,
     * 1X,'Geschrieben: ',I10,/,
     * 1X,'    Gelesen: ',I10)
            STOP 'Fehler'
         ENDIF
 30   CONTINUE
C
      CALL QFCLOS (NT, 0)
      WRITE (NTTERM, 7002) NT
 7002 FORMAT (1X, 'Unit ',I3,' von QFCLOS geschlossen.')
C
C     Nochmal aufmachen
      NT=QFREEU()
      WRITE (NTTERM,6006) NT
 6006 FORMAT (1X,'Neue von QFREEU gewaehlte Unit :',I3,/)
C
      CALL QFCRSE (NT, ZEILE, 'UN', NERR)
      IF (NERR .NE. 0) THEN
         WRITE (NTTERM, '(A)') ' QFCRSE liefert einen Fehler.'
         STOP 'Fehler'
      ELSE
         WRITE (NTTERM, '(A)') ' QFCRSE fehlerfrei.'
      ENDIF
      IF (.NOT.QFILDA (ZEILE)) THEN
         WRITE (NTTERM, '(A)') ' QFILDA liefert den falschen '//
     $        'Wert .FALSE.'
         WRITE (NTTERM, '(A)') ZEILE
         STOP 'Fehler'
      ENDIF
      DO 40 I=1, 20
         WRITE (NT, ERR=9010) I
 40   CONTINUE
      REWIND NT
      DO 50 I=1, 20
         READ (NT, ERR=9011, END=9011) J
         IF (I.NE.J) THEN
            WRITE (NTTERM, 8912) I, J
 8912 FORMAT(1X,/,
     * 1X,'Vergleich beim binaeren sequentiellen Lesen falsch.',/,
     * 1X,'Geschrieben: ',I10,/,
     * 1X,'    Gelesen: ',I10)
            STOP 'Fehler'
         ENDIF
 50   CONTINUE
      CALL QFCLOS (NT, 0)
      NT=QFREEU()
      LR = QFRECL (0, NDWEG, NRWEG)
      WRITE (NTTERM, 7914) LR
 7914 FORMAT (1X, 'LR (Recordlaenge von QFRECL) :',I7,/)
C
      WRITE (NTTERM,*)
      WRITE (NTTERM,'(A)') ' Weiter mit RETURN...'
      READ  (NTTAST,4711) WEITER
 4711 FORMAT (A1)
C
C     Direct Access unformatiert
C     (vorher alles schliessen)
C
      CALL QFCLAL
C
      WRITE (NTTERM,'(A)') ' Test von QFCRDI, '//
     $     'Direct Access unformatiert...'
      CALL QFCRDI (NT, ZEILE, 'UN', LR, NERR)
      IF (NERR .GT. 0) THEN
         WRITE (NTTERM, 8104)
 8104    FORMAT (1X,'QFCRDI liefert einen Fehler aus QFDEL,',/,
     *           1X,'Datei existiert schon und ist offen.')
         STOP 'Fehler'
      ENDIF
C
      WRITE (NTTERM, 8011) NT, LR, NDWEG, NRWEG
 8011 FORMAT(1X,/,
     *     1X,'Parameter-Werte fuer DIRECT ACCESS, unformatiert:',/,
     *     1X,'           UNIT: ',I6,/,
     *     1X,'   RECORDLAENGE: ',I6,/,
     *     1X,'DBLE PRO RECORD: ',I6,/,
     *     1X,' ANZAHL RECORDS: ',I6)
            WRITE (NTTERM, 8975)
 8975       FORMAT(1X,/,
     *           1X,'Schreiben mit QFWRDI.')
      DO 70 IRWEG=1, NRWEG
         DO 60 IDWEG=1, NDWEG
            DWEG(IDWEG)=100*IRWEG+IDWEG
 60      CONTINUE
         CALL QFWRDI (NT, IRWEG, DWEG, LR, NERR)
         IF (NERR .NE. 0) THEN
            WRITE (NTTERM, 8915)
 8915       FORMAT(1X,/,
     *           1X,'QFWRDI liefert einen Fehler.')
            STOP 'Fehler'
         ENDIF
 70   CONTINUE
      WRITE (NTTERM, 7060)
 7060 FORMAT (1X,'QFWRDI in Ordnung',/)
C
            WRITE (NTTERM, 8976)
 8976       FORMAT(1X,/,
     *           1X,'Lesen mit QFREDI.')
      DO 90 IRWEG=1, NRWEG
         CALL QFREDI (NT, IRWEG, DWEG, LR, NERR)
         IF (NERR .GT. 0) THEN
            WRITE (NTTERM, 8916)
 8916       FORMAT(1X,/,
     *           1X,'QFREDI liefert einen Fehler.')
            STOP 'Fehler'
         ENDIF
         DO 80 IDWEG=1, NDWEG
            DSOLL=100*IRWEG+IDWEG
            IF (DWEG(IDWEG).NE.DSOLL) THEN
               WRITE (NTTERM, 8917) DSOLL, DWEG(IDWEG)
 8917          FORMAT(1X,/,
     *              1X,'Vergleich beim Direct-Access-Lesen falsch.',/,
     *              1X,'Geschrieben: ',D12.5,/,
     *              1X,'    Gelesen: ',D12.5)
               STOP 'Fehler'
            ENDIF
 80      CONTINUE
 90   CONTINUE
      WRITE (NTTERM, 7061)
 7061 FORMAT (1X,'QFREDI in Ordnung',/)
C
C
      IF (QIYENO ('Soll die Datei geloescht werden ? (J/N)',
     $     0, 'J')) THEN
         CALL QFCLOS (NT, -1)
         WRITE (NTTERM,7062) NT
 7062    FORMAT (1X, 'Unit ',I3, ' geschlossen und geloescht',/)
      ELSE
         CALL QFCLOS (NT, 0)
         WRITE (NTTERM,7063) NT
 7063    FORMAT (1X, 'Unit ',I3, ' geschlossen',/)
      ENDIF
C
C   Test von QFILAO, QFSIZE, QSDEFI, QSDEFO
C
C     Aufruf von QSDEFI und QSDEFO
      WRITE (NTTERM,7260)
 7260 FORMAT (1X, 'Aufruf von QFILAO, QSDEFI und QSDEFO',/)
C
      NTI=QFREEU()
      CALL QFSIZE (1000)
      CALL QFILAO (NTI, FILNA1, 'CR', 'FO', 0,
     *     ' Wie soll die Eingabedatei heissen ?',
     *     0, CHAR(0), 1, NERR)
C
      IF (NERR.GT.0) WRITE (NTTERM,'(A)') ' Fehler in QFILAO'
      IF (NERR.LT.0) WRITE (NTTERM,'(A)') ' Warning in QFILAO'
C
      NTO=QFREEU()
      FILNA2=QICHAR ('Wie soll die Ausgabedatei heissen?',
     $                0, CHAR(0), L, NERR)
C
      CALL QSDEFI (NTI, FILNA1, NERR)
      IF (NERR.NE.0) WRITE (*,'(A)') ' Fehler in QSDEFI'
      WRITE (*,'(A)') ' QSDEFI fertig'
      CALL QSDEFO (NTO, FILNA2, NERR)
      IF (NERR.NE.0) WRITE (*,'(A)') ' Fehler in QSDEFO'
      WRITE (*,'(A)') ' QSDEFO fertig'
         WRITE (*,7160)
 7160    FORMAT (1X, 'Tastatur mit QSDEFI substituiert ')
         WRITE (*,7161)
 7161    FORMAT (1X,'Bildschirm mit QSDEFO substituiert ')
C
C      Zuruecksetzen auf Terminal und Bildschirm
C
      CALL QSDEFO (NTTERM, CHAR(0), NERR)
      CALL QSDEFI (NTTAST, CHAR(0), NERR)
C
C   Oeffnen mehrerer Dateien (als Scratch-File, ohne Dateinamen)
C   fuer den Aufruf von QFCLAL
      NT=QFREEU()
      OPEN (NT,STATUS='SCRATCH')
      NT=QFREEU()
      OPEN (NT,STATUS='SCRATCH')
      NT=QFREEU()
      OPEN (NT,STATUS='SCRATCH')
      NT=QFREEU()
      OPEN (NT,STATUS='SCRATCH')
C
      WRITE (NTTERM,7064)
 7064 FORMAT (1X,/,1X, 'Schliessen aller geoeffneten Dateien')
      CALL QFCLAL
      WRITE (NTTERM,7065)
 7065 FORMAT (1X, 'QFCLAL - alle geoeffneten Dateien geschlossen',/)
C     
      WRITE (NTTERM,'(A)') ' Weiter mit RETURN...'
      READ  (NTTAST,4711) WEITER
C
      RETURN
C
 9011 WRITE (NTTERM, 8911)
      GOTO 9100
 9010 WRITE (NTTERM, 8910)
      GOTO 9100
 9008 WRITE (NTTERM, 8908)
      GOTO 9100
 9007 WRITE (NTTERM, 8907)
      GOTO 9100
 9100 STOP 'FEHLER'
 8907 FORMAT(1X,/,
     * 1X,'ERR beim formatierten, sequentiellen Schreiben.')
 8908 FORMAT(1X,/,
     * 1X,'ERR oder END beim formatierten, sequentiellen Lesen.')
 8910 FORMAT(1X,/,
     * 1X,'ERR beim binaeren, sequentiellen Schreiben.')
 8911 FORMAT(1X,/,
     * 1X,'ERR oder END beim binaeren, sequentiellen Lesen.')
C
      END
C
C     ############################################################
C
      FUNCTION MENUE2 (VORGAB, AZMENP, MENPKT)
C     --------------------------------------
C     Subroutine fuer FVA-Programm DRESP
C     hier benutzt fuer das Programm LibTest
C     --------------------------------------
C
C     Unterprogramm zur Ausgabe eines Menues und Erfragen der Antwort
C
      INCLUDE 'imelibh'
C
      INTEGER AZMENP, LENANT, NERR, I, CRTLIN
      CHARACTER*2 VORGAB, ANTTXT, MENUE2
      CHARACTER*75 MENPKT(AZMENP+1)
      CHARACTER*80 QICHAR
C
      CRTLIN = 24
      CALL LINIE2
      CALL LINIE2
      CALL LINIE2
C
C     Erzeugen des Menues
C
      CALL LINIE0
      CALL LINIE2
C     Titel schreiben
      WRITE (NTTERM, 10) MENPKT(1)
 10   FORMAT (9X, A)
      CALL LINIE1
      CALL LINIE2
C     Alle Menuepunkte schreiben
      DO 100 I=1, AZMENP
         WRITE(NTTERM, 20) MENPKT(I+1)
 20      FORMAT(4X, A)
 100  CONTINUE
      CALL LINIE1
C     Leerzeilen schreiben:
      DO 110 I=1, (CRTLIN -9 -AZMENP)
         WRITE(NTTERM, *)
 110  CONTINUE
C     Antwort einlesen...
      ANTTXT=QICHAR(' Geben Sie bitte den gewuenschten Menuepunkt ein',
     &     0, VORGAB, LENANT, NERR)
C     Bei Fehler in Eingaberoutine nichts zurueckgeben
      IF (NERR .NE. 0) THEN
         ANTTXT=' '
      ENDIF
      CALL QCTOUP (ANTTXT)
C     ... und zurueckgeben:
      MENUE2=ANTTXT
C     Vorher aber noch zwei Leerzeilen schreiben:
      WRITE (NTTERM, *)
      WRITE (NTTERM, *)
      RETURN
C     Fertig!
      END
C
C     ############################################################
C
      SUBROUTINE LINIE0
C     -----------------------------------
C     Subroutine fuer FVA-Programm DRESP
C     -----------------------------------
      INCLUDE 'imelibh'
C
      WRITE(NTTERM,1)
    1 FORMAT(2X, 58('='),' LIBTEST 1.05')
      RETURN
      END
C
C     ############################################################
C
      SUBROUTINE LINIE1
C     -----------------------------------
C     Subroutine fuer FVA-Programm DRESP
C     -----------------------------------
      INCLUDE 'imelibh'
C
      WRITE(NTTERM,1)
    1 FORMAT(2X, 76('='))
      RETURN
      END
C
C     ############################################################
C
      SUBROUTINE LINIE2
C     -----------------------------------
C     Subroutine fuer FVA-Programm DRESP
C     -----------------------------------
      INCLUDE 'imelibh'
C
      WRITE(NTTERM,1)
    1 FORMAT(2X,' ')
      RETURN
      END
