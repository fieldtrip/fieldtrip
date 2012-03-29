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
C
C             FILE: FOFULIB.F77
C            AUTOR: HEINZ WALTERMANN
C         REVISION: 1.0
C LETZTE AENDERUNG: 01.06.1991
C
C
C
C FOFULIB:
C ========
C
C IN DIESER LIBRARY SIND ALLE ROUTINEN GESAMMELT, DIE SICH AUF
C DIE FEM-FORMFUNKTIONEN ODER IHRE ABLEITUNGEN BEZIEHEN.
C
C
C *******************************
C * K O M P L E T  STANDARD F77 *
C *******************************
C
C
      SUBROUTINE FOFINI(RELIN,REINT)                                    NEU
C
C#
C
C FOFINI:
C =======
C
C INITIALISIERT DIE FOFULIB
C
C VARIABLE TYP ART DIMENSION      ERKLAERUNG
C ----------------------------------------------------------------------
C RELIN     L   I                 .T.
C                                 LINEARE ELEMENTE WERDEN GEGENUEBER DEM
C                                 IN DER SUBROUTINE ITPANZ GEFORDERTEN
C                                 INTEGRATIONSGRAD UM EINE ORDNUNG ZURUECK
C                                 GESTUFT. DER UEBERGABEPARAMETER WIRD DANN
C                                 ENTSPRECHEND KORRIGIERT.
C REINT     L   I                 .T.
C                                 ALLE NORMALEN ELEMENTE WERDEN UM EINE
C                                 STUFE IM INTEGRATIONSGRAD GEAENDERT.
C
C#
C
      INCLUDE 'fofulib.inc'
      LOGICAL RELIN,REINT

      REDLIN=RELIN
      REDINT=REINT
      R E T U R N
      END
