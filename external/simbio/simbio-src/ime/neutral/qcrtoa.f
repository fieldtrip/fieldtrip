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
      SUBROUTINE QCRTOA (RA, NN, NE, TEX, LE)
C     Wandelt die REAL-Zahl RA in einen entsprechenden Textstring TEX
C     der Laenge LE Buchstaben um.
C
C     ****************
C     * STANDARD F77 *
C     ****************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C RA        D   I                 Umzuwandelnde Zahl
C NN        I   I                 Anzahl der Nachkommastellen in der Text-
C                                 darstellung
C NE        I   I                 10-hoch-Exponent, durch den RA vor der
C                                 Umwandlung geteilt werden soll.
C TEX       C   O *(*)            Text-Darstellung der Zahl RA
C LE        I   O                 Laenge von TEX in Buchstaben
C
C
      INCLUDE 'imelibh'
C
      CHARACTER*(*) TEX
      CHARACTER*80  BUF
      CHARACTER*10  FMT, CH
      INTEGER NN, NE, LE, LCH, IP, I
      REAL RA, RAH, R10
C
      PARAMETER (R10=10.E0)
C
      RAH=RA/R10**NE
      FMT=' '
      BUF=' '
      TEX=' '
      FMT(1:5)='(F80.'
      CALL QCITOA(NN,CH,LCH)
      if (lch .gt. 4) then
         lch=1
         ch=' '
      endif
      IP=5
      FMT(IP+1:IP+LCH)=CH(1:LCH)
      IP=IP+LCH
      FMT(IP+1:IP+1)=')'
      WRITE(BUF,FMT)RAH
      DO 1 I=1,80
         IF(BUF(I:I).NE.' ')THEN
            IP=I
            GOTO 2
         END IF
    1 CONTINUE
      IP=80
      BUF(IP:IP)='0'
    2 LE=80-IP+1
      IF (BUF(80:80) .EQ. '.') LE=LE-1
      TEX(1:LE)=BUF(IP:IP+LE-1)
      R E T U R N
      END

