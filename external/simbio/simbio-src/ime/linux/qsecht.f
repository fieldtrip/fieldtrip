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
      DOUBLE PRECISION FUNCTION QSECHT (I, J, IERR)
C
C     Ermittlung der echten Zeit in Sekunden
C
C     ***********************
C     * RECHNER - ABHAENGIG *
C     ***********************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C I         I   I                 NUMMER DER STOP-UHR (MAX. IMAX)
C J         I   I                 0 : STARTEN DER STOP-UHR (INIT)
C                                 1 : ABFRAGEN DES STANDES DER STOP-UHR
C IERR      I   O                 FEHLERPARAMETER
C
C
C
      INCLUDE 'imelibh'
      EXTERNAL TIME
C
      INTEGER IMAX
      PARAMETER (IMAX=20)                  
      DOUBLE PRECISION  TE(IMAX), TA(IMAX) 
      INTEGER I, J, IERR, INI(IMAX)
      INTEGER *4 TIME
C
      SAVE INI,TA,TE                    
C                                       
      DATA INI/IMAX*0/                  
C                                       
      QSECHT=0.0D0                      
      IF(I.LT.0.OR.I.GT.IMAX)GOTO 9000  
      IF(J.NE.1.AND.J.NE.0)  GOTO 9000  
      IF(J.EQ.0)THEN                    
        TA(I)=TIME()
        INI(I)=1                        
      ELSE                              
        IF(INI(I).EQ.0)GOTO 9000        
        TE(I)=TIME()
        QSECHT=TE(I)-TA(I)              
      END IF                            
      IERR=0                            
      R E T U R N                       
 9000 IERR=1                            
      R E T U R N                       
      END
