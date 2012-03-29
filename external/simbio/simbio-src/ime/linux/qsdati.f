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
      SUBROUTINE QSDATI (DATUM, UHRZ)
C     Bestimmt das Datum und die Uhrzeit des Aufrufes
C
C     ***********************
C     * RECHNER - ABHAENGIG *
C     ***********************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C DATUM     C   O  *8             Datum des Aufrufes dieser Routine
C                                 TT.MM.JJ
C UHRZ      C   O  *8             Uhrzeit des Aufrufes dieser Routine
C                                 HH:MM:SS
C
      INCLUDE 'imelibh'
C
      INTEGER TDATE(3)
      INTEGER IDAY, IMON, IYEAR, I
      CHARACTER*(*) DATUM, UHRZ
      CHARACTER*9 DUMMY
      CHARACTER*2 CH                                                    AIX
C
      TDATE(1)=IMON
      TDATE(2)=IDAY
      TDATE(3)=IYEAR
      CALL IDATE(TDATE)
      DATUM='  .  .  '                                                  AIX
      WRITE(DATUM(1:2),'(I2)')IDAY                                      AIX
      WRITE(DATUM(4:5),'(I2)')IMON                                      AIX
      WRITE(DATUM(7:8),'(I2)')IYEAR                                     AIX
      DO 10 I=1,8                                                       AIX
         IF (DATUM(I:I) .EQ. ' ') THEN                                  AIX
            DATUM(I:I)='0'                                              AIX
         END IF                                                         AIX
   10 CONTINUE                                                          AIX
      CALL AIXTIME(DUMMY)
      UHRZ=DUMMY(1:8)                                                   AIX

      if (len(datum).gt.8) datum(9:)=' '
      if (len(uhrz).gt.8)  uhrz(9:)=' '

      R E T U R N
      END
