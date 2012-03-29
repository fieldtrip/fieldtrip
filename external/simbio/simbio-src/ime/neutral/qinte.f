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
      INTEGER FUNCTION QINTE
     * (TEXT,IRERAN,IDFWE,FMDFWE,IRA,IRAU,IRAO,IERR)
C     Liest eine Integerzahl ein
C
C     ****************
C     * STANDARD F77 *
C     ****************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C TEXT      C   I  *(*)           Auszugebender Text
C IRERAN    I   I                 Ausrichtung des Textes
C                                     0 : linksbuendig
C                                 sonst : rechtsbuendig in Spalte IRERAN
C IDFWE     I   I                 Defaultwert der Antwort
C FMDFWE    C   I  *80            'NO'  : kein Default zulaessig
C                                 sonst : Format, in dem IDFWE ausgegeben
C                                         werden soll (z.B.: '(I5)')
C IRA       I   I                 Range-Kennung
C                                 -1 : erlaubt IRAU<=QDINTI
C                                  0 : keine Grenzen da
C                                  1 : erlaubt       QDINTI<=IRAO
C                                  2 : erlaubt IRAU<=QDINTI<=IRAO
C IRAU      I   I                 untere Grenze
C IRAO      I   I                 obere Grenze
C IERR      I   O                 Fehlervariable
C                                 =1: Fehler in QICHAR
C                                 =2: Benutzerabbruch
C
      INCLUDE 'imelibh'
C
      CHARACTER*(*) TEXT
      CHARACTER*80 FMDFWE,CDFWE,QICHAR,LESTEX
      INTEGER IRERAN, IDFWE, IRA, IRAU, IRAO,
     $     IBUFF, LTI, IERR
      REAL BUFFER
      LOGICAL OKEY
C
C     Default-Wert in Text umwandeln, um QICHAR aufzurufen.
C     Wenn kein Default:
      IF (FMDFWE(1:2) .EQ. 'NO') THEN
         CDFWE=CHAR(0)
      ELSE
C     Default ausgeben
         CDFWE=' '
         WRITE(CDFWE,FMDFWE,ERR=10)IDFWE
         GOTO 20
 10      CDFWE='>> FEHLER in QINTE <<'
      END IF
C
C     Lese-Routine QICHAR aufrufen
 20   LESTEX=QICHAR(TEXT,IRERAN,CDFWE,LTI, IERR)
C
C     Abfrage, ob kein Fehler bzw. Benutzerabbruch:
      IF (IERR .EQ. 1) THEN
cc         WRITE (NTTERM,*) 'FEHLER in der Routine QICHAR'
cc         STOP
         RETURN
      ENDIF
C     Bei Abbruch durch Benutzer:
      IF (IERR .EQ. 2) THEN
cc         WRITE (NTTERM,*) 'Benutzerabbruch in QINTE '
         QINTE=0
         RETURN
      ENDIF
C
C     Wert von Text in Zahl umformen:
C
cc      READ(LESTEX,'(BN,F30.10)',ERR=20) BUFFER
      READ(LESTEX,'(BN,F80.0)',ERR=20) BUFFER
      IBUFF=INT(BUFFER)
C
C     Check, ob IBUFF im Range-Bereich ist, falls RANGE gesetzt
      OKEY=.TRUE.
      IF (IRA.EQ.-1) THEN
C
C     Untere Grenze ist gegeben
C
         IF (IBUFF.LT.IRAU) OKEY=.FALSE.
      ELSE IF(IRA.EQ.1)THEN
C
C     Obere Grenze ist gegeben
C
         IF (IBUFF.GT.IRAO) OKEY=.FALSE.
      ELSE IF (IRA.EQ.2) THEN
C
C     Intervall ist gegeben
C
         IF (IBUFF.LT.IRAU.OR.IBUFF.GT.IRAO) OKEY=.FALSE.
      END IF
      IF (.NOT.OKEY) THEN
C
C     Range ist nicht okey, Fehlermeldung
C
        IF (IRA.EQ.-1) THEN
           WRITE (NTTERM,8010) IBUFF, IRAU
        ELSE IF (IRA.EQ.1) THEN
           WRITE(NTTERM, 8020) IBUFF, IRAO
        ELSE IF(IRA.EQ.2)THEN
           WRITE (NTTERM, 8030) IBUFF, IRAU, IRAO
        END IF
        GOTO 20
      END IF
C
C     Umspeichern auf Function-Variable
C
      IERR=0
      QINTE = IBUFF
      R E T U R N
 8010 FORMAT(/,1X,'Die gelesene Zahl: ',I20,/,
     *         1X,'       muss ueber: ',I20,/,
     *         1X,'liegen.',/)
 8020 FORMAT(/,1X,'Die gelesene Zahl: ',I20,/,
     *         1X,'       muss unter: ',I20,/,
     *         1X,'liegen.',/)
 8030 FORMAT(/,1X,'Die gelesene Zahl: ',I20,/,
     *         1X,'   muss  zwischen: ',I20,/,
     *         1X,'              und: ',I20,/,
     *         1X,'liegen.',/)
      END
