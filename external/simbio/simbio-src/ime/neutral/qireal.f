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
      REAL FUNCTION QIREAL (TEXT,
     $     IRERAN, RDFWE, FMDFWE, IRA, RRAU, RRAO, IERR)
C     Liest eine REAL-Zahl ein
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
C RDFWE     D   I                 Defaultwert der Antwort
C FMDFWE    C   I  *80            'NO'  : kein Default zulaessig
C                                 sonst : Format, in dem IDFWE ausgegeben
C                                         werden soll (z.B.: '(E12.5)'
C IRA       I   I                 Range-Kennung
C                                 -1 : erlaubt DRAU<=DBLINP
C                                  0 : keine Grenzen da
C                                  1 : erlaubt       DBLINP<=DRAO
C                                  2 : erlaubt DRAU<=DBLINP<=DRAO
C RRAU      D   I                 Untere Grenze
C RRAO      D   I                 Obere Grenze
C IERR      I   O                 Fehlercode
C
C
      INCLUDE 'imelibh'
C
cc      INTEGER LAEVRE
C
      CHARACTER*(*) TEXT
      CHARACTER*80 FMDFWE, QICHAR, LESTEX, CDFWE
      LOGICAL OKEY
      INTEGER IRA, IRERAN, LTI, IERR
      REAL RDFWE, RRAU, RRAO, BUFFER
C
      IF(FMDFWE(1:2).EQ.'NO')THEN
         CDFWE=CHAR(0)
      ELSE
         CDFWE=' '
         WRITE (CDFWE, FMDFWE, ERR=10) RDFWE
         GOTO 20
 10      CDFWE='>> FEHLER in QIREAL <<'
      END IF
 20   LESTEX=QICHAR (TEXT,IRERAN,CDFWE,LTI,IERR)
      IF (IERR .EQ. 1) THEN
         WRITE (NTTERM,*) 'FEHLER in der Routine QICHAR'
         STOP
      ENDIF
C     Bei Abbruch durch Benutzer:
      IF (IERR .EQ. 2) THEN
         QIREAL=0.0
         RETURN
      ENDIF
C
      READ(LESTEX,'(BN,F80.0)',ERR=20)BUFFER
      OKEY=.TRUE.
      IF(IRA.EQ.-1)THEN
         IF(BUFFER.LT.RRAU)OKEY=.FALSE.
      ELSE IF(IRA.EQ.1)THEN
         IF(BUFFER.GT.RRAO)OKEY=.FALSE.
      ELSE IF(IRA.EQ.2)THEN
         IF(BUFFER.LT.RRAU.OR.BUFFER.GT.RRAO)OKEY=.FALSE.
      END IF
      IF(.NOT.OKEY)THEN
         IF(IRA.EQ.-1)THEN
            WRITE(NTTERM,8010)BUFFER,RRAU
         ELSE IF(IRA.EQ.1)THEN
            WRITE(NTTERM,8020)BUFFER,RRAO
         ELSE
            WRITE(NTTERM,8030)BUFFER,RRAU,RRAO
         END IF
         GOTO 20
      END IF
      QIREAL=BUFFER
      IERR=0
      R E T U R N
 8010 FORMAT(/,1X,'Die gelesene Zahl: ',G20.5,/,
     *         1X,'       muss ueber: ',G20.5,/,
     *         1X,'liegen.',/)
 8020 FORMAT(/,1X,'Die gelesene Zahl: ',G20.5,/,
     *         1X,'       muss unter: ',G20.5,/,
     *         1X,'liegen.',/)
 8030 FORMAT(/,1X,'Die gelesene Zahl: ',G20.5,/,
     *         1X,'   muss  zwischen: ',G20.5,/,
     *         1X,'              und: ',G20.5,/,
     *         1X,'liegen.',/)
      END
