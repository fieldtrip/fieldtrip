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
      SUBROUTINE QLBINI (NTASZU, NTERZU, NERRZU, CABB, LABB)
C     Initialisierung fuer IMELIB
C
C     ****************
C     * STANDARD F77 *
C     ****************
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C NTASZU    I   O                 FORTRAN-Unit der Tastatur
C NTERZU    I   O                 FORTRAN-Unit des Bildschirms
C NERRZU    I   O                 FORTRAN-Unit fuer Fehlermeldungen
C CABB      C  I/O    *(*)        Abbruch-String
C LABB      I  I/O                Laenge des Abbruchstrings
C                                 Wird fuer LABB der Wert 0 uebergeben,
C                                 so nimmt die IMELIB den Defaultwert,
C                                 ansonsten wird der String genommen der
C                                 mit CABB uebergeben wird. ACHTUNG uebergeben
C                                 Sie NIEMALS Konstanten fuer LABB und CABB
C                                 sondern immer nur Variablen
C
      INCLUDE 'imelibh'
C
      INTEGER I, NTASZU, NTERZU, NERRZU, LABB, ILISIN
      integer qfreeu
      CHARACTER *(*) CABB
      SAVE ILISIN
      DATA ILISIN /0/
C
C     Markiere alle Fortran-Units als frei.
C     Sollen irgendwelche Units fuer das Programm unerreichbar sein,
C     so sind die entsprechenden Komponenten des Vektors IFIKEN
C     mit 3 zu besetzen.
C     0: Unit ist besetzt
C     1: Unit ist reserviert
C     2: Unit ist frei
C     3: Unit ist gesperrt
C
C     Alle Units unterhalb IFIMIP sperren (in diesem Bereich liegen auch
C     die Units fuer Tastatur und Terminal):
C
      IF (ILISIN.EQ.0) THEN
         CALL QSPARA
         IF (MANFIP.GT.1000) THEN
            STOP 'Bitte Feld IFIKEN in imelibh vergroessern.'
         endif
         ILISIN=1
         NTPROT=0
         DO 10 I= 1, (IFIMIP-1)
            IFIKEN(I)=3
 10      CONTINUE
C
C     Alle Units ab IFIMIP bis MANFIP sollen frei sein:
C
         DO 20 I= IFIMIP, MANFIP
            IFIKEN(I)=2
 20      CONTINUE
C
C
C     Zur Eingabe:
C     a) Eingabe ueber Tastatur
C     b) Fortran-Unit der Tastatur ist 5
         TAST   = .TRUE.
         NTTAST = ISTIN
C
C     Zur Ausgabe:
C     a) Ausgabe auf Terminal
C     b) Fortran-Unit des Terminals ist 6
         TERM   = .TRUE.
         NTTERM = ISTOUT
C     
C     Fehlermeldungen auf das Terminal:
C     Fortran-Unit fuer Fehlermeldungen ist 6 
         NTERRO = ISTERR
C
C     Hier wird der String definiert, den der Benutzer fuer
C     einen Abbruch der Eingaben eingeben darf.
C
         IF (LABB .EQ. 0) THEN
            LABB=2
            CABB='$q'
         ENDIF
C     Sonst den String nehmen, der uebergeben wurde:
         CABBRU=CABB
         LABBRU=LABB
      ELSE
         CABB=CABBRU
         LABB=LABBRU
      ENDIF
C
C     Umspeichern zu Uebergabe:
C
      NTASZU=NTTAST
      NTERZU=NTTERM
      NERRZU=NTERRO

C Variablen loeschen

      do i=1,maxvar
         lnvar(i)=0
         ltvar(i)=0
      enddo

c History-Datei oeffen

      nthist=0
      iwhist=0

      RETURN
      END
