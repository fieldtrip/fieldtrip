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
      SUBROUTINE QDTEXT (TEXT, IFORT)
C     Starten eines Dialogtextes. Schreibt den Text TEXT auf das
C     Terminal. Es wird der Zeilenvorschub so gesteuert, dass die
C     von der Subroutine DLGTXS produzierten Texte direkt unterein-
C     ander stehen, wenn jeweils zwischen zwei Aufrufen von DLGTXS
C     ein Read-Statement steht. Bei der PRIME wird durch betaetigen
C     der RETURN-Taste ein CR und LF ausgeloest, so dass DLGTXS nur
C     den reinen Text ausgeben muss, ohne irgendwelche Vorschubsteue-
C     rung. Da einige Terminals das erste Zeichen als Druckersteuer-
c     zeichen interpretieren, gibt DLGTXS ggf. vor dem eigentlichen
C     Text ein Blank aus.
C
C     ***********************
C     * RECHNER - ABHAENGIG *
C     ***********************
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C TEXT      C   I  *(*)           Auszugebender Text
C IFORT     I   I                 Flag, ob Ausgabe gestartet (IFORT=0)
C                                 oder fortgesetzt werden soll (IFORT=1)
C                                 oder beendet werden soll (IFORT=2)
C
C
      INCLUDE 'imelibh'
C
      CHARACTER TEXT*(*)
      INTEGER LETE, IFORT, QCLEN
C
C
C     Check, ob Textlaenge OK
      LETE = QCLEN (TEXT)
      IF (LETE .GT. 80) LETE=80
C
C     TEXT ausgeben:
C
C     Abfrage, ob neue Zeile
C
      IF (LETE .EQ. 0) LETE=1
      IF (IFORT.EQ.1) THEN
         WRITE (NTTERM,8000) TEXT(1:LETE)
 8000    FORMAT (A,$)
      ELSE if (ifort.eq.0) then
         WRITE (NTTERM,8001) TEXT(1:LETE)
 8001    FORMAT (/,A,$)
      ELSE if (ifort.eq.2) then
         WRITE (NTTERM,8002) TEXT(1:LETE)
 8002    FORMAT (A)
      ENDIF
C
      RETURN
      END
