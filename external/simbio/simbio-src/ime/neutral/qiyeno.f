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
      LOGICAL FUNCTION QIYENO (TEXT, IRERAN, CDFWE)
C     Bringt eine Ja/Nein Frage auf das Terminal, liest die Antwort
C     ein und formt diese in einen Logical-Wert um.
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
C CDFWE     C   I  *(*)           Defaultwert der Antwort
C                                 Wenn CDFWE=CHAR(0), dann kein Default
C
      INCLUDE 'imelibh'
C
      INTEGER IRERAN, LTI, IERR
      CHARACTER*(*) TEXT, CDFWE
      CHARACTER*80  QICHAR, BUFFER
C
   10 CONTINUE
      BUFFER=QICHAR (TEXT,IRERAN,CDFWE,LTI,IERR)
      IF(LTI.EQ.0)GOTO 10
      CALL QCTOUP(BUFFER)
      IF(LTI.GT.4)GOTO 10
      IF(BUFFER(1:4).EQ.'Y   '.OR.BUFFER(1:4).EQ.'YE  '.OR.
     *   BUFFER(1:4).EQ.'YES '.OR.BUFFER(1:4).EQ.'J   '.OR.
     *   BUFFER(1:4).EQ.'JA  ')THEN
         QIYENO=.TRUE.
      ELSE IF(BUFFER(1:4).EQ.'N   '.OR.BUFFER(1:4).EQ.'NO  '.OR.
     *        BUFFER(1:4).EQ.'NE  '.OR.BUFFER(1:4).EQ.'NEI '.OR.
     *        BUFFER(1:4).EQ.'NEIN')THEN
         QIYENO=.FALSE.
      ELSE
         GOTO 10
      END IF
      RETURN
      END
