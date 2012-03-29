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
c Datei und Pfadnamen zusammensetzen, doppelte // werden bearbeitet
c=======================================================================
      subroutine qsmnam(name,path,base)
      include 'imelibh'
      character*(*) name,path,base
      character*(maxlen) cpath,cbase
      integer qclen

      cbase=base
      cpath=path
      lb=qclen(cbase)
      lp=qclen(cpath)
      ln=lp
      name=cpath
      if (lp.gt.0) then
         ln=ln+1
         name(ln:ln)=patsep
      endif
      name(ln+1:)=cbase
      ln=ln+lb
      call delspa(name(1:ln),ln)
      do i=1,ln-1
         if (name(i:i).eq.patsep .and. name(i+1:i+1).eq.patsep) then
            name(1:i)=' '
         endif
      enddo
      call delspa(name(1:ln),ln)

      return
      end

