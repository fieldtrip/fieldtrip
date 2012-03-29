!  $1 21/03/2003  Anwander A.  released for SIMBIO

!---->---1---------2---------3---------4---------5---------6---------7--<
!
!     findobfl.f :
!                              -------------------
!     begin                : Mai 2000
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

      subroutine findobfl(inelem,inelty,inelne,knobfl,knobty,knobne,
     *                    knobsc,nextob,lastob,neknob,knobbe,knoban,
     *                    istack,ialles, nelem, npoin,nobfl ,maobfl,
     *                    mxknin)
c
c
c
c findobfl:
c =======
c
c ermittelt die oberflaechen einer struktur-struktur (ialles=0)
c oder
c alle in der struktur vorkommenden flaechen (ialles=1)
c
c
c variable typ art dimension      erklaerung
c ----------------------------------------------------------------------
c inelem    i   i  mxkn,nelem     node-element order
c inelty    i   i  nelem          type of the elements
c inelne    i   i  nelem          necke of the elements
c knobfl    i   o  mknpfp,maobfl  knoten-element-zuordnung der oberflaechen
c knobty    i   o  maobfl         typen der oberflaechen
c knobne    i   o  maobfl         necke der oberflaechen
c knobsc    i   o  maobfl         0 : flaeche kommt aus einem volumenelement
c                                 1 : flaeche kommt aus einem schalenelement
c nextob    i   o  maobfl         nextob(i) ist die nummer der auf die ober-
c                                 flaeche i folgenden oberflaeche
c lastob    i   o  maobfl         lastob(i) ist die nummer der der oberflaeche
c                                 i vorausgehenden oberflaeche
c neknob    i   o  maobfl         nextob, nur bezogen auf den anfangsknoten
c                                 der flaechen
c knobbe    i   o  npoin          position in knobfl, an der die erste flaeche
c                                 mit anfangsknoten ipoin steht
c knoban    i   o  npoin          anzahl der flaechen mit anfangsknoten ipoin
c
c istack    i   w  npoin          work-feld, das von der routine flglei ge-
c                                 braucht wird, um die gleichheit zweier
c                                 flaechen zu ermitteln
c ialles    i   i                 =0 : nur oberflaechen suchen
c                                 =1 : alle flaechen suchen
c nelem     i   i                 number of elements in the volume
c npoin     i   i                 number of nodes in the volume
c nobfl     i   o                 number of elements in the surface structure 
c maobfl    i   i                 max number of surfaces (size of the fields)
c mxknin    i   i                 max nodes per element in the input element list
c
c
      character*6 ch
c
      parameter (mknpep = 32   ,
     *           meltyp =  6   ,
     *           mknpfp = 12   ,
     *           mkelem = mknpep+2,
     *           mflpep =   6     ,
     *           mgelfp =  17     )

      dimension 
     * lfelnu(10    , 3), nknpel(meltyp), nhkpel(meltyp),
     * nkapel(   meltyp), nflpel(meltyp), lfltyp(mflpep,meltyp),
     * loffse(   meltyp), knelfl(mknpfp,mgelfp)

      dimension
     * inelem(mxknin,nelem ),inelty(nelem ),inelne(nelem ),
     * knobfl(mknpfp,maobfl),knobty(maobfl),knobne(maobfl),
     * nextob(       maobfl),lastob(maobfl),neknob(maobfl),
     * knobbe(       npoin ),knoban(npoin ),kn    (mkelem),
     * knfl  (       mknpfp),knobsc(maobfl),istack(npoin )
  
      data lfelnu/
     *  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     *  2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     *  4, 5, 6, 0, 0, 0, 0, 0, 0, 0/

      data nknpel/
     *  4, 9,12,16,24,32/
      data nhkpel/
     *  2, 3, 4, 4, 6, 8/
      data nkapel/
     *  1, 3, 4, 6, 9,12/

c anzahl der flaechen in abhaengigkeit der laufenden elementnummer

      data nflpel/
     *  0, 1, 1, 4, 5, 6/

c element-typen der einzelnen flaechen

      data lfltyp/
     *  0, 0, 0, 0, 0, 0,
     *  2, 0, 0, 0, 0, 0,
     * 12, 0, 0, 0, 0, 0,
     *  2, 2, 2, 2, 0, 0,
     *  2, 2,12,12,12, 0,
     * 12,12,12,12,12,12/

c offset der flaechen-definitionen in knelfl-matrix fuer die elementtypen

      data loffse/
     *  0, 0, 1, 2, 6,11/

c knoten-element-zuordnung der element-flaechen

      data knelfl/
     *  1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0,
     *  1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,
     *  1, 2, 3, 5, 6, 7,11,12,13, 0, 0, 0,
     *  2, 4, 3, 9,16,12,15,10, 6, 0, 0, 0,
     *  4, 1, 3,14,13,10, 8, 7,16, 0, 0, 0,
     *  1, 4, 2, 8,15,11,14, 9, 5, 0, 0, 0,
     *  1, 3, 2,18,17,16, 9, 8, 7, 0, 0, 0,
     *  4, 5, 6,10,11,12,19,20,21, 0, 0, 0,
     *  1, 4, 6, 3,13,21,24, 9,22,12,15,18,
     *  3, 6, 5, 2,15,20,23, 8,24,11,14,17,
     *  1, 2, 5, 4, 7,14,19,22,16,23,10,13,
     *  1, 4, 3, 2,24,23,22,21,12,11,10, 9,
     *  4, 8, 7, 3,20,27,31,11,32,15,19,23,
     *  8, 5, 6, 7,16,13,14,15,28,25,26,27,
     *  5, 1, 2, 6,29, 9,18,25,17,21,30,13,
     *  1, 5, 8, 4,17,28,32,12,29,16,20,24,
     *  2, 3, 7, 6,10,19,26,30,22,31,14,18/
c


c initialisierung

      nobfl=0
      minobf=1
      maxobf=0
      minloc=1
c
      do 50 iobfl=1,maobfl
        nextob(iobfl)=0
        lastob(iobfl)=0
        neknob(iobfl)=0
        knobsc(iobfl)=0
   50 continue
      do 60 ipoin=1,npoin
        knobbe(ipoin)=0
        knoban(ipoin)=0
        istack(ipoin)=0
   60 continue
c

c main-loop ueber alle elemente

      do 210 ielem=1,nelem

c knoten-element-zuordnung, typ und anzahl der knoten einlesen

      do 61 iecke=1,mkelem
         kn(iecke)=0
   61 continue
      do 62 iecke=1,inelne(ielem)
         kn(iecke)=inelem(iecke,ielem)
   62 continue
  
      ieltyp=inelty(ielem)
      necke =inelne(ielem)
      write(ch,'(    I6)')ieltyp
      read (ch,'(BN,3X,3I1)')igldim,ilonum,ilodim
      launum=lfelnu(ilonum+1,ilodim)

c bestimme nun die anzahl der flaechen dieses elementes

      nelfl=nflpel(launum)

c stelle fest, ob dies eine schale ist (nelfl=1), oder nicht

      if(nelfl.eq.1)then
        ischal=1
      else
        ischal=0
      end if

c loop ueber alle flaechen des elementes

      do 200 ielfl=1,nelfl

c knoten-element-zuordnung, typ und anzahl der knoten der flaeche einlesen
c dazu:
c position der knoten-element-zuordnung in knelfl bestimmen

      ipos=loffse(launum)+ielfl

c typ der flaeche bestimmen

      ifltyp=lfltyp(ielfl,launum)
      write(ch,'(    I6)')ifltyp
      read (ch,'(BN,3X,3I1)')igldim,ilonum,ilodim
      lfnuf=lfelnu(ilonum+1,ilodim)

c es werden nur flaechige flaechen beachtet

      if(ilodim.ne.2)goto 200

c knoten-element-zuordnung der flaeche bestimmen

      do 70 iecke=1,nknpel(lfnuf)
      knfl(iecke)=kn(knelfl(iecke,ipos))
   70 continue
      do 80 iecke=nknpel(lfnuf)+1,mknpfp
      knfl(iecke)=0
   80 continue

c bei kanten, die zwar quadratische, aber keine kubischen zwischen-
c knoten enthalten, die positionen korrigieren

      ioffse=nhkpel(lfnuf)
      ikan  =nkapel(lfnuf)
      do 90 i=ioffse+1,ioffse+ikan
      if(knfl(i).eq.0.and.knfl(i+ikan).ne.0)then
        knfl(i     )=knfl(i+ikan)
        knfl(i+ikan)=0
      end if
   90 continue

c bestimme nun die hoechste besetzte knotennummer

      do 100 iecke=nknpel(lfnuf),1,-1
      if(knfl(iecke).ne.0)goto 110
  100 continue
      goto 9002
  110 neckfl=iecke

c ermittlung der nummer des kleinsten hauptknotens

      if(ilonum.eq.0)then

c dreiecke

        knoten=min(knfl(1),knfl(2),knfl(3))
      else if(ilonum.eq.1)then

c vierecke

        knoten=min(knfl(1),knfl(2),knfl(3),knfl(4))
      else
        goto 9003
      end if

c iflda = 1 eine flaeche mit diesem anfangsknoten ist bereits da

      iflda=knobbe(knoten)
      if(iflda.eq.0)then

c es gibt noch keine flaeche mit diesem anfangsknoten

        if(nobfl.eq.maobfl)then

c zu viele flaechen

          goto 9004
        else if(nobfl.eq.1)then

c es ist erst eine flaeche vorhanden

          lastlo=minloc-1
          if(minloc.eq.1)lastlo=maxobf
          lastob(minloc)=lastlo
          lastob(lastlo)=minloc
          nextob(lastlo)=minloc
          nextob(minloc)=lastlo
        else if(nobfl.gt.1)then

c es sind schon mehrere flaechen da

          lastlo=minloc-1
          if(minloc.eq.1)lastlo=maxobf
          nextlo=nextob(lastlo)
          nextob(minloc)=nextob(lastlo)
          nextob(lastlo)=minloc
          lastob(minloc)=lastlo
          lastob(nextlo)=minloc
        end if

c daten der flaeche speichern

        do 120 ikn=1,mknpfp
           knobfl(ikn,minloc)=knfl(ikn)
  120   continue
        knobty(minloc)=ifltyp
        knobne(minloc)=neckfl
        knobsc(minloc)=ischal
        knobbe(knoten)=minloc
        knoban(knoten)=1

c position der letzten flaeche aktualisieren

        if(minloc.gt.maxobf)maxobf=minloc

c position der erste flaeche aktualisieren

        if(minloc.lt.minobf)minobf=minloc

c minloc aktualisieren

        do 130 i=minloc+1,maobfl
           if(nextob(i).eq.0)goto 140
  130   continue

c alle positionen sind besetzt, setze minloc auf null

        i=0

c bei i wurde ein freie position gefunden

  140   continue
        minloc=i

c anzahl der flaechen aktualisieren

        nobfl=nobfl+1
      else

c es gibt bereits mind. eine flaeche mit diesem anfangsknoten
c nfla : anzahl der bereits gespeicherten flaechen

        if(ischal.eq.1)then

c bei schalen wird kein test auf gleichheit durchgefuehrt

          igl=0
          ilocsa=0
          iloc=knobbe(knoten)
          nfla=knoban(knoten)
          do 141 ifla=1,nfla
             ilocsa=iloc
             iloc=neknob(iloc)
  141     continue
        else
          nfla=knoban(knoten)
          ilocsa=0
          iloc=knobbe(knoten)
          igl=0
          do 144 ifla=1,nfla
            if(knobsc(iloc).eq.1)then

c schalen sind immer oberflaechen

            igl=0
          else
            call sameplane(  knfl     ,ifltyp      ,neckfl      ,
     *                  knobfl(1,iloc),knobty(iloc),knobne(iloc),
     *                  istack        ,ifla        ,nfla        ,
     *                  igl           )
          end if
          if(igl.eq.1)goto 145
          ilocsa=iloc
          iloc=neknob(iloc)
  144     continue

c die beiden flaechen sind nicht gleich, setze igl auf 0

          igl=0
  145     continue
        end if
        if(igl.eq.0)then

c dies ist eine neue flaeche mit einem alten anfangsknoten

          if(nobfl.eq.maobfl)then

c zu viele flaechen

            goto 9004
          else if(nobfl.eq.1)then

c es ist erst eine flaeche vorhanden

            lastlo=ilocsa
            lastob(minloc)=lastlo
            lastob(lastlo)=minloc
            nextob(lastlo)=minloc
            nextob(minloc)=lastlo
          else if(nobfl.gt.1)then

c es sind schon mehrere flaechen da

            lastlo=ilocsa
            nextlo=nextob(lastlo)
            nextob(minloc)=nextob(lastlo)
            nextob(lastlo)=minloc
            lastob(minloc)=lastlo
            lastob(nextlo)=minloc
          end if

c daten der flaeche speichern

          do 150 ikn=1,mknpfp
          knobfl(ikn,minloc)=knfl(ikn)
  150     continue
          knobty(minloc)=ifltyp
          knobne(minloc)=neckfl
          knobsc(minloc)=ischal
          knoban(knoten)=knoban(knoten)+1
          neknob(ilocsa)=minloc

c position der letzten flaeche aktualisieren

          if(minloc.gt.maxobf)maxobf=minloc

c position der erste flaeche aktualisieren

          if(minloc.lt.minobf)minobf=minloc

c minloc aktualisieren

          do 160 i=minloc+1,maobfl
          if(nextob(i).eq.0)goto 170
  160     continue

c alle positionen sind besetzt, setze minloc auf null

          i=0

c bei i wurde ein freie position gefunden

  170     continue
          minloc=i

c anzahl der flaechen aktualisieren

          nobfl=nobfl+1
        else if(ialles.eq.0)then

c die beiden flaechen sind gleich, loesche die bereits gespeicherte
c flaeche, falls nur oberflaechen gesucht werden.

          if(nobfl.eq.1)then

c es ist nur eine flaeche gespeichert, loesche diese

            nextob(iloc)=0
            lastob(iloc)=0
            minobf=1
            maxobf=0
            minloc=1
            nobfl =0
            knobbe(knoten)=0
            knoban(knoten)=0
            neknob(knoten)=0
          else

c es sind noch mind. zwei flaechen im speicher

            lastlo=lastob(iloc)
            nextlo=nextob(iloc)
            lastob(iloc)=0
            nextob(iloc)=0
            lastob(nextlo)=lastlo
            nextob(lastlo)=nextlo
            nobfl=nobfl-1
            if(ilocsa.eq.0)then

c loesche die erste flaeche fuer diesen knoten

              iwo=knobbe(knoten)
              knobbe(knoten)=neknob(iwo)
              neknob(   iwo)=0
              knoban(knoten)=knoban(knoten)-1
            else if(ifla.eq.nfla)then

c loesche die letzte flaeche fuer diesen knoten

              neknob(ilocsa)=0
              knoban(knoten)=knoban(knoten)-1
            else

c loesche eine flaeche mitten raus

              neknob(ilocsa)=neknob(iloc)
              neknob(  iloc)=0
              knoban(knoten)=knoban(knoten)-1
            end if
          end if

c position der letzten flaeche aktualisieren

          if(iloc.eq.maxobf)then
            if(nobfl.eq.0)then
              maxobf=0
            else
  180         maxobf=maxobf-1
              if(nextob(maxobf).eq.0)goto 180
            end if
          end if

c position der ersten flaeche aktualisieren

          if(iloc.eq.minobf)then
            if(nobfl.eq.0)then
              minobf=1
            else
  190         minobf=minobf+1
              if(nextob(minobf).eq.0)goto 190
            end if
          end if

c minloc aktualisieren

          if(iloc.lt.minloc)minloc=iloc
        end if
      end if

c nehme nun die naechste flaeche dieses elementes

  200 continue

c nehme nun das naechste element

  210 continue

c komprimiere nun die gefundenen flaechen

      iloc=minobf-1
      do 250 iobfl=1,nobfl
  230 continue
      iloc=iloc+1
      if(nextob(iloc).eq.0)goto 230
      do 240 ikn=1,mknpfp
      knobfl(ikn,iobfl)=knobfl(ikn,iloc)
  240 continue
      knobty(iobfl)=knobty(iloc)
      knobne(iobfl)=knobne(iloc)
      knobsc(iobfl)=knobsc(iloc)
  250 continue

c add the spatial dimension to the type
      do 260 iobfl=1,nobfl
         ityp=300+knobty(iobfl)
         knobty(iobfl)=ityp
  260 continue


      r e t u r n
c
 9002 write(7,8902)
      goto 9100
 9003 write(7,8903)ifltyp
      goto 9100
 9004 write(7,8904)
      goto 9100
 9100 r e t u r n
 8902 format(1x,'ES WURDE EINE NULL-FLAECHE GEFUNDEN!')
 8903 format(1x,'ILLEGALER FLAECHENTYP: ',i6)
 8904 format(1x,'ZU VIELE FLAECHEN (PARAMETER MIWRKP ERHOEHEN)')
      end

      
  
      subroutine sameplane(knfl1 ,iflty1,neckf1,knfl2 ,iflty2,  
     *                     neckf2,istack,itest ,ntest ,iglei ) 
c
c#
c
c sameplane:
c =======
c
c testet 2 flaechen auf gleichheit
c
c variable typ art dimension      erklaerung
c ----------------------------------------------------------------------
c knfl1     i   i  mknpfp         knoten-element-zuordnung flaeche 1
c iflty1    i   i                 typ der flaeche 1
c neckf1    i   i                 hoechste besetzte knotennummer flaeche 1
c knfl2     i   i                                                        2
c iflty2    i   i                                                        2
c neckf2    i   i                                                        2
c istack    i  i/w npoin          istack(i) ist eins fuer einen knoten, der
c                                 in der flaeche 1 enthalten ist, sonst null
c                                 bei der uebergabe muss istack komplett null
c                                 sein, falls itest=1 ist. istack wird nach
c                                 der bearbeitung (itest=ntest) an den ver-
c                                 wendeten stellen genullt.
c itest     i   i                 laufende testnummer fuer diese referenz-
c                                 flaeche (knfl1)
c ntest     i   i                 anzahl der tests mit dieser knf1-flaeche
c iglei     i   o                 0 : flaechen sind nicht gleich
c                                 1 : flaechen sind gleich
c
c#
c
c
      character*6 ch
c
      parameter (meltyp =  6 , mknpfp = 12 )
c
      dimension lfelnu(10 , 3), nhkpel(meltyp)
c
      dimension knfl1(mknpfp),knfl2(mknpfp),istack(*)
c
      save lfnu,nhaukn
c
      data lfelnu/
     *  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     *  2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     *  4, 5, 6, 0, 0, 0, 0, 0, 0, 0/

      data nhkpel/
     *  2, 3, 4, 4, 6, 8/
c

c setze kennung auf ungleich

      iglei=0
      if(itest.eq.1)then

c erster test: initialisierung
c laufende nummer des typs und anzahl hauptknoten bestimmen
c istack besetzen

        write(ch,'(    I6)')iflty1
        read (ch,'(4X,2I1)')ilonum,ilodim
        lfnu=lfelnu(ilonum+1,ilodim)
        nhaukn=nhkpel(lfnu)
        do 10 ihaukn=1,nhaukn
        istack(knfl1(ihaukn))=1
   10   continue
      end if

c gleiche flaechen haben gleichen typ und gleiche hoechste knotennummer

      if(iflty1.eq.iflty2.and.neckf1.eq.neckf2)then
        isum=0
        do 20 ihaukn=1,nhaukn
        isum=isum+istack(knfl2(ihaukn))
   20   continue
        if(isum.eq.nhaukn)iglei=1
      end if
      if(itest.eq.ntest.or.iglei.eq.1)then

c letzter test oder flaechen gleich, istack freigeben

        do 30 ihaukn=1,nhaukn
        istack(knfl1(ihaukn))=0
   30   continue
      end if
      r e t u r n
      end

