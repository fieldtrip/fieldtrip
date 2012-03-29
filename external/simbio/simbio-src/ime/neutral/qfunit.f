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
c@======================================================================
c Sperren bzw. Freigeben von Fortran-Units
c
c csteu*1  I      = 'k' : Kleinste Fortran-Unit in IUNIT zurueckgeben
c                 = 'g' : groesste Fortran-Unit in IUNIT zurueckgeben
c                 = 's' : IUNIT sperren
c                 = 'f' : IUNIT freigeben
c iunit    I/O    siehe oben
c ierr     O      = 0 : Alles klar
c                 = 1 : IUNIT < IFIMIP oder IUNIT > MANFIP
c                 =-1 : zu sperrende (freizugebende) Unit war vorher
c                       nicht gesperrt (freigegeben)
c
c=======================================================================
      subroutine qfunit(csteu,iunit,ierr)
      include 'imelibh'
      character*1 csteu
      integer iunit,ierr

      ierr=0
      if (csteu.eq.'k' .or. csteu.eq.'K') then
         iunit=ifimip
         return
      endif
      if (csteu.eq.'g' .or. csteu.eq.'G') then
         iunit=manfip
         return
      endif
      if (csteu.eq.'s' .or. csteu.eq.'S') then
         if (iunit.ge.ifimip .and. iunit.le.manfip) then
            if (ifiken(iunit).ne.2) ierr=-1
            ifiken(iunit)=3
         else
            ierr=1
         endif
         return
      endif
      if (csteu.eq.'f' .or. csteu.eq.'F') then
         if (iunit.ge.ifimip .and. iunit.le.manfip) then
            if (ifiken(iunit).ne.3) ierr=-1
            ifiken(iunit)=2
         else
            ierr=1
         endif
         return
      endif

      return
      end
