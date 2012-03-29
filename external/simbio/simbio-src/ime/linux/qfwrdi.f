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
      SUBROUTINE QFWRDI (IUNIT, IREC, IVEC, IR, IERR)
C
C     Schreibt logischen Record auf Direct-Access-File
C
C     ***********************
C     * RECHNER - ABHAENGIG *
C     ***********************
C
C
C ACHTUNG:
C ========
C Die Variable IVEC stellt einen Speicherbereich des Rechners dar.
C Es ist daher unbedingt erforderlich, dass die Laenge des Datentypes
C dieses Vektors mit der Einheit der Recordlaenge im OPEN-Statement
C uebereinstimmt.
C
C Beispiele:
C PRIME 2450 unter PRIMOS:
C OPEN-Statement verwendet als Einheit 16-BIT-Woerter
C -> Datentyp IVEC ist INTEGER*2
C
C IBM3090 unter VM/XA:
C OPEN-Statement verwendet als Einheit BYTES
C -> Datentyp IVEC ist LOGICAL*1
C
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C IUNIT     I   I                 FORTRAN-UNIT DES FILES
C IREC      I   I                 RECORDNUMMER, IN DEN GESCHR. WERDEN SOLL
C                                 0 : ES ERFOLGT EIN BINAERES SEQUENTIELLES
C                                     WRITE
C IVEC      I   I  IR             ZU SPEICHERNDE WERTE
C                                 IVEC MUSS EINEN DATENTYP HABEN, DESSEN
C                                 LAENGE GLEICH DER EINHEIT DER RECL-OPTION
C                                 IST. (PRIME = 16 BIT)
C IR        I   I                 RECORDLAENGE DES FILES
C IERR      I   O                 FEHLERCODE
C
C
      INCLUDE 'imelibh'
C
      INTEGER   IREC, IR, IUNIT, I, IERR
c$$$      LOGICAL*1 IVEC(IR)
      CHARACTER*1 IVEC(IR)      ! War frueher mal LOGICAL*1
C
C
C     Schreiben mit Standard F77-Routinen
C     Voraussetzung:
C     Daten-Laenge IVEC = Einheit der RECL-Option!!
C
      IERR=0
C
      IF (IREC.EQ.0) THEN
         WRITE (IUNIT,ERR=9000) (IVEC(I),I=1,IR)
      ELSE
         WRITE (IUNIT,REC=IREC,ERR=9000) (IVEC(I),I=1,IR)
      END IF
      R E T U R N
C
 9000 IERR=1
      R E T U R N
      END
 
 
 
