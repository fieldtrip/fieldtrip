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
      SUBROUTINE QSUSER (LOGNA, LELO)
C     Gibt den Login-Namen des Users zurueck
C
C     ***********************
C     * RECHNER - ABHAENGIG *
C     ***********************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C LOGNA     C   O  *80            LOGIN-Name des Users
C LELO      I   O                 Laenge von LOGNA in Buchstaben
C
      INCLUDE 'imelibh'
      EXTERNAL CUSERID
C
      INTEGER LELO, LMAX, IPO, ICHAR
      CHARACTER*(*) LOGNA
      CHARACTER*2  CH                                                   AIX
      CHARACTER*81 TMPSTR                                               AIX
C
      CALL CUSERID(TMPSTR)                                              AIX
      DO 10 IPO=1,80                                                    AIX
         IF(ICHAR(TMPSTR(IPO:IPO)).EQ.0) GOTO 20                        AIX
   10 CONTINUE                                                          AIX
      IPO=81                                                            AIX
   20 CONTINUE                                                          AIX
      LELO=IPO-1                                                        AIX
      LOGNA=TMPSTR(1:LELO)                                              AIX
      R E T U R N
      END
