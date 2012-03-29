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
      SUBROUTINE QFOPDI (IUNIT ,FILNAM, STAT, FM, IR, IERR)
C     Oeffnet einen Direct-Access-File und ordnet ihn
C     der Fortran-Unit IUNIT zu.
C
C     ***********************
C     * RECHNER - ABHAENGIG *
C     ***********************
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C IUNIT     I   I                 FORTRAN-Unit, der der File zugeordn. wird
C FILNAM    C   I *(*)            Name des Files
C STAT      C   I *(*)            STATUS des Files
C                                 'OL D'
C                                 'NE W'
C                                 'UN KNOWN'
C                                 'SC RATCH'
C FM        C   I *(*)            Formatierung der Datei
C                                 'FO rmatted'
C                                 'UN formatted'
C IR        I   I                 Recordlaenge, die in der RECL-Option des
C                                 OPEN-Statements stehen muss
C IERR      I   O                 Fehlercode
C
C
      INCLUDE 'imelibh'
C
      INTEGER IUNIT, IR, IERR, IOERR
      CHARACTER*8   STA
      CHARACTER*11  FMT
      CHARACTER*(*) STAT, FM, FILNAM
C
C
      STA=' '
      FMT=' '
      IOERR=0
      IF(STAT(1:1).EQ.'O'.or.stat(1:1).eq.'o')THEN
        STA(1:3)='OLD'                                                  CHK
      ELSE IF(STAT(1:1).EQ.'N'.or.stat(1:1).eq.'n')THEN
        STA(1:3)='NEW'                                                  CHK
      ELSE IF(STAT(1:1).EQ.'S'.or.stat(1:1).eq.'s')THEN
        STA(1:7)='SCRATCH'                                              CHK
      ELSE IF(STAT(1:1).EQ.'U'.or.stat(1:1).eq.'u')THEN
        STA(1:7)='UNKNOWN'                                              CHK
      END IF
C
      IF(FM(1:2).EQ.'FO'.or.fm(1:2).eq.'fo')THEN
         FMT(1:9)='FORMATTED'                                           CHK
      ELSE IF(FM(1:2).EQ.'UN'.or.fm(1:2).eq.'un')THEN
         FMT(1:11)='UNFORMATTED'                                        CHK
      END IF
      IF(STA(1:7).EQ.'SCRATCH')THEN
        OPEN(IUNIT                        ,                             CHK
     *       STATUS      = STA            ,                             CHK
     *       ACCESS      = 'DIRECT'       ,                             CHK
     *       FORM        = FMT            ,                             CHK
     *       RECL        = IR             ,                             CHK
     *       IOSTAT      = IOERR          )                             CHK
      ELSE
        OPEN(IUNIT                        ,                             CHK
     *       FILE        = FILNAM         ,                             CHK
     *       STATUS      = STA            ,                             CHK
     *       ACCESS      = 'DIRECT'       ,                             CHK
     *       FORM        = FMT            ,                             CHK
     *       RECL        = IR             ,                             CHK
     *       IOSTAT      = IOERR          )                             CHK
      END IF
C
      IF (IOERR .NE. 0) THEN
         WRITE(NTERRO, 200) IOERR
 200     FORMAT (1X, '(IMELIB-QFOPDI:) F E H L E R: IOSTAT= ', I3)
         IERR=1
         R E T U R N
      ENDIF
C
      IF (IFIKEN(IUNIT) .NE. 3) IFIKEN(IUNIT)=0
      IERR=0
      R E T U R N
      END
