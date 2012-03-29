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
      DOUBLE PRECISION FUNCTION QIDBL (TEXT,
     $     IRERAN, DDFWE, FMDFWE, IRA, DRAU, DRAO, IERR)
C     Liest eine Doublezahl ein
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
C                                 sonst : rechtsbuendig in spalte IRERAN
C DDFWE     D   I                 Defaultwert der Antwort
C FMDFWE    C   I  *80            'NO'  : kein Default zulaessig
C                                 sonst : Format, in dem IDFWE ausgegeben
C                                         werden soll (z.B.: '(D12.5)'
C IRA       I   I                 Range-Kennung
C                                 -1 : erlaubt DRAU<=QDBLIN
C                                  0 : keine Grenzen da
C                                  1 : erlaubt       QDBLIN<=DRAO
C                                  2 : erlaubt DRAU<=QDBLIN<=DRAO
C DRAU      D   I                 untere Grenze
C DRAO      D   I                 obere Grenze
C IERR      I   O                 Fehlercode
C
C
      INCLUDE 'imelibh'
C
      CHARACTER*(*) TEXT
      CHARACTER*80 FMDFWE, CDFWE, QICHAR, LESTEX
      LOGICAL OKEY
      INTEGER IRA, IRERAN, LTI, IERR
      DOUBLE PRECISION DDFWE, DRAU, DRAO, BUFFER
C
      IF(FMDFWE(1:2).EQ.'NO')THEN
         CDFWE=CHAR(0)
      ELSE
         CDFWE=' '
         WRITE (CDFWE, FMDFWE, ERR=10) DDFWE
         GOTO 20
 10      CDFWE='>> Fehler in QIDBL <<'
      END IF
   20 LESTEX= QICHAR (TEXT, IRERAN, CDFWE, LTI, IERR)
C     Bei Abbruch durch Benutzer:
      IF (IERR .EQ. 2) THEN
         QIDBL=0.0D0
         RETURN
      ENDIF
C
      READ(LESTEX,'(BN,D80.0)',ERR=20)BUFFER
      OKEY=.TRUE.
      IF(IRA.EQ.-1)THEN
         IF(BUFFER.LT.DRAU)OKEY=.FALSE.
      ELSE IF(IRA.EQ.1)THEN
         IF(BUFFER.GT.DRAO)OKEY=.FALSE.
      ELSE IF(IRA.EQ.2)THEN
         IF(BUFFER.LT.DRAU.OR.BUFFER.GT.DRAO)OKEY=.FALSE.
      END IF
      IF(.NOT.OKEY)THEN
         IF(IRA.EQ.-1)THEN
            WRITE(NTTERM,8010)BUFFER,DRAU
         ELSE IF(IRA.EQ.1)THEN
            WRITE(NTTERM,8020)BUFFER,DRAO
         ELSE
            WRITE(NTTERM,8030)BUFFER,DRAU,DRAO
         END IF
         GOTO 20
      END IF
      QIDBL=BUFFER
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
