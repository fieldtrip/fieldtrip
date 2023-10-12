function [transform] = ft_affinecoordinates(original, target, varargin)

% FT_AFFINECOORDINATES returns the affine coordinate transformation matrix that
% converts FROM a specific head coordinate TO a specific head coordinate system.
%
% Use as
%   [transform] = ft_affinecoordinates(from, to)
%
% Note that translations are expressed in millimeters, therefore the geometrical data
% to which this coordinate transformation is applied must also be specified in
% millimeters.
%
% See also FT_CONVERT_COORDSYS, FT_CONVERT_UNITS, FT_HEADCOORDINATES, FT_WARP_APPLY

% Copyright (C) 2005-2022, Robert Oostenveld & Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are the 48 generic axis orientation triplets, these specify the axes but no origin
%   a = anterior
%   p = posterior
%   l = left
%   r = right
%   s = superior
%   i = inferior

% these are for speeding up subsequent calls
persistent initialized
persistent ctf2ctf ctf2bti ctf2fourd ctf2yokogawa ctf2eeglab ctf2neuromag ctf2itab ctf2acpc ctf2spm ctf2mni ctf2fsaverage ctf2tal ctf2scanras ctf2scanlps ctf2dicom ctf2als ctf2ali ctf2ars ctf2ari ctf2pls ctf2pli ctf2prs ctf2pri ctf2las ctf2lai ctf2ras ctf2rai ctf2lps ctf2lpi ctf2rps ctf2rpi ctf2asl ctf2ail ctf2asr ctf2air ctf2psl ctf2pil ctf2psr ctf2pir ctf2sal ctf2ial ctf2sar ctf2iar ctf2spl ctf2ipl ctf2spr ctf2ipr ctf2sla ctf2ila ctf2sra ctf2ira ctf2slp ctf2ilp ctf2srp ctf2irp ctf2lsa ctf2lia ctf2rsa ctf2ria ctf2lsp ctf2lip ctf2rsp ctf2rip
persistent bti2ctf bti2bti bti2fourd bti2yokogawa bti2eeglab bti2neuromag bti2itab bti2acpc bti2spm bti2mni bti2fsaverage bti2tal bti2scanras bti2scanlps bti2dicom bti2als bti2ali bti2ars bti2ari bti2pls bti2pli bti2prs bti2pri bti2las bti2lai bti2ras bti2rai bti2lps bti2lpi bti2rps bti2rpi bti2asl bti2ail bti2asr bti2air bti2psl bti2pil bti2psr bti2pir bti2sal bti2ial bti2sar bti2iar bti2spl bti2ipl bti2spr bti2ipr bti2sla bti2ila bti2sra bti2ira bti2slp bti2ilp bti2srp bti2irp bti2lsa bti2lia bti2rsa bti2ria bti2lsp bti2lip bti2rsp bti2rip
persistent fourd2ctf fourd2bti fourd2fourd fourd2yokogawa fourd2eeglab fourd2neuromag fourd2itab fourd2acpc fourd2spm fourd2mni fourd2fsaverage fourd2tal fourd2scanras fourd2scanlps fourd2dicom fourd2als fourd2ali fourd2ars fourd2ari fourd2pls fourd2pli fourd2prs fourd2pri fourd2las fourd2lai fourd2ras fourd2rai fourd2lps fourd2lpi fourd2rps fourd2rpi fourd2asl fourd2ail fourd2asr fourd2air fourd2psl fourd2pil fourd2psr fourd2pir fourd2sal fourd2ial fourd2sar fourd2iar fourd2spl fourd2ipl fourd2spr fourd2ipr fourd2sla fourd2ila fourd2sra fourd2ira fourd2slp fourd2ilp fourd2srp fourd2irp fourd2lsa fourd2lia fourd2rsa fourd2ria fourd2lsp fourd2lip fourd2rsp fourd2rip
persistent yokogawa2ctf yokogawa2bti yokogawa2fourd yokogawa2yokogawa yokogawa2eeglab yokogawa2neuromag yokogawa2itab yokogawa2acpc yokogawa2spm yokogawa2mni yokogawa2fsaverage yokogawa2tal yokogawa2scanras yokogawa2scanlps yokogawa2dicom yokogawa2als yokogawa2ali yokogawa2ars yokogawa2ari yokogawa2pls yokogawa2pli yokogawa2prs yokogawa2pri yokogawa2las yokogawa2lai yokogawa2ras yokogawa2rai yokogawa2lps yokogawa2lpi yokogawa2rps yokogawa2rpi yokogawa2asl yokogawa2ail yokogawa2asr yokogawa2air yokogawa2psl yokogawa2pil yokogawa2psr yokogawa2pir yokogawa2sal yokogawa2ial yokogawa2sar yokogawa2iar yokogawa2spl yokogawa2ipl yokogawa2spr yokogawa2ipr yokogawa2sla yokogawa2ila yokogawa2sra yokogawa2ira yokogawa2slp yokogawa2ilp yokogawa2srp yokogawa2irp yokogawa2lsa yokogawa2lia yokogawa2rsa yokogawa2ria yokogawa2lsp yokogawa2lip yokogawa2rsp yokogawa2rip
persistent eeglab2ctf eeglab2bti eeglab2fourd eeglab2yokogawa eeglab2eeglab eeglab2neuromag eeglab2itab eeglab2acpc eeglab2spm eeglab2mni eeglab2fsaverage eeglab2tal eeglab2scanras eeglab2scanlps eeglab2dicom eeglab2als eeglab2ali eeglab2ars eeglab2ari eeglab2pls eeglab2pli eeglab2prs eeglab2pri eeglab2las eeglab2lai eeglab2ras eeglab2rai eeglab2lps eeglab2lpi eeglab2rps eeglab2rpi eeglab2asl eeglab2ail eeglab2asr eeglab2air eeglab2psl eeglab2pil eeglab2psr eeglab2pir eeglab2sal eeglab2ial eeglab2sar eeglab2iar eeglab2spl eeglab2ipl eeglab2spr eeglab2ipr eeglab2sla eeglab2ila eeglab2sra eeglab2ira eeglab2slp eeglab2ilp eeglab2srp eeglab2irp eeglab2lsa eeglab2lia eeglab2rsa eeglab2ria eeglab2lsp eeglab2lip eeglab2rsp eeglab2rip
persistent neuromag2ctf neuromag2bti neuromag2fourd neuromag2yokogawa neuromag2eeglab neuromag2neuromag neuromag2itab neuromag2acpc neuromag2spm neuromag2mni neuromag2fsaverage neuromag2tal neuromag2scanras neuromag2scanlps neuromag2dicom neuromag2als neuromag2ali neuromag2ars neuromag2ari neuromag2pls neuromag2pli neuromag2prs neuromag2pri neuromag2las neuromag2lai neuromag2ras neuromag2rai neuromag2lps neuromag2lpi neuromag2rps neuromag2rpi neuromag2asl neuromag2ail neuromag2asr neuromag2air neuromag2psl neuromag2pil neuromag2psr neuromag2pir neuromag2sal neuromag2ial neuromag2sar neuromag2iar neuromag2spl neuromag2ipl neuromag2spr neuromag2ipr neuromag2sla neuromag2ila neuromag2sra neuromag2ira neuromag2slp neuromag2ilp neuromag2srp neuromag2irp neuromag2lsa neuromag2lia neuromag2rsa neuromag2ria neuromag2lsp neuromag2lip neuromag2rsp neuromag2rip
persistent itab2ctf itab2bti itab2fourd itab2yokogawa itab2eeglab itab2neuromag itab2itab itab2acpc itab2spm itab2mni itab2fsaverage itab2tal itab2scanras itab2scanlps itab2dicom itab2als itab2ali itab2ars itab2ari itab2pls itab2pli itab2prs itab2pri itab2las itab2lai itab2ras itab2rai itab2lps itab2lpi itab2rps itab2rpi itab2asl itab2ail itab2asr itab2air itab2psl itab2pil itab2psr itab2pir itab2sal itab2ial itab2sar itab2iar itab2spl itab2ipl itab2spr itab2ipr itab2sla itab2ila itab2sra itab2ira itab2slp itab2ilp itab2srp itab2irp itab2lsa itab2lia itab2rsa itab2ria itab2lsp itab2lip itab2rsp itab2rip
persistent acpc2ctf acpc2bti acpc2fourd acpc2yokogawa acpc2eeglab acpc2neuromag acpc2itab acpc2acpc acpc2spm acpc2mni acpc2fsaverage acpc2tal acpc2scanras acpc2scanlps acpc2dicom acpc2als acpc2ali acpc2ars acpc2ari acpc2pls acpc2pli acpc2prs acpc2pri acpc2las acpc2lai acpc2ras acpc2rai acpc2lps acpc2lpi acpc2rps acpc2rpi acpc2asl acpc2ail acpc2asr acpc2air acpc2psl acpc2pil acpc2psr acpc2pir acpc2sal acpc2ial acpc2sar acpc2iar acpc2spl acpc2ipl acpc2spr acpc2ipr acpc2sla acpc2ila acpc2sra acpc2ira acpc2slp acpc2ilp acpc2srp acpc2irp acpc2lsa acpc2lia acpc2rsa acpc2ria acpc2lsp acpc2lip acpc2rsp acpc2rip
persistent spm2ctf spm2bti spm2fourd spm2yokogawa spm2eeglab spm2neuromag spm2itab spm2acpc spm2spm spm2mni spm2fsaverage spm2tal spm2scanras spm2scanlps spm2dicom spm2als spm2ali spm2ars spm2ari spm2pls spm2pli spm2prs spm2pri spm2las spm2lai spm2ras spm2rai spm2lps spm2lpi spm2rps spm2rpi spm2asl spm2ail spm2asr spm2air spm2psl spm2pil spm2psr spm2pir spm2sal spm2ial spm2sar spm2iar spm2spl spm2ipl spm2spr spm2ipr spm2sla spm2ila spm2sra spm2ira spm2slp spm2ilp spm2srp spm2irp spm2lsa spm2lia spm2rsa spm2ria spm2lsp spm2lip spm2rsp spm2rip
persistent mni2ctf mni2bti mni2fourd mni2yokogawa mni2eeglab mni2neuromag mni2itab mni2acpc mni2spm mni2mni mni2fsaverage mni2tal mni2scanras mni2scanlps mni2dicom mni2als mni2ali mni2ars mni2ari mni2pls mni2pli mni2prs mni2pri mni2las mni2lai mni2ras mni2rai mni2lps mni2lpi mni2rps mni2rpi mni2asl mni2ail mni2asr mni2air mni2psl mni2pil mni2psr mni2pir mni2sal mni2ial mni2sar mni2iar mni2spl mni2ipl mni2spr mni2ipr mni2sla mni2ila mni2sra mni2ira mni2slp mni2ilp mni2srp mni2irp mni2lsa mni2lia mni2rsa mni2ria mni2lsp mni2lip mni2rsp mni2rip
persistent fsaverage2ctf fsaverage2bti fsaverage2fourd fsaverage2yokogawa fsaverage2eeglab fsaverage2neuromag fsaverage2itab fsaverage2acpc fsaverage2spm fsaverage2mni fsaverage2fsaverage fsaverage2tal fsaverage2scanras fsaverage2scanlps fsaverage2dicom fsaverage2als fsaverage2ali fsaverage2ars fsaverage2ari fsaverage2pls fsaverage2pli fsaverage2prs fsaverage2pri fsaverage2las fsaverage2lai fsaverage2ras fsaverage2rai fsaverage2lps fsaverage2lpi fsaverage2rps fsaverage2rpi fsaverage2asl fsaverage2ail fsaverage2asr fsaverage2air fsaverage2psl fsaverage2pil fsaverage2psr fsaverage2pir fsaverage2sal fsaverage2ial fsaverage2sar fsaverage2iar fsaverage2spl fsaverage2ipl fsaverage2spr fsaverage2ipr fsaverage2sla fsaverage2ila fsaverage2sra fsaverage2ira fsaverage2slp fsaverage2ilp fsaverage2srp fsaverage2irp fsaverage2lsa fsaverage2lia fsaverage2rsa fsaverage2ria fsaverage2lsp fsaverage2lip fsaverage2rsp fsaverage2rip
persistent tal2ctf tal2bti tal2fourd tal2yokogawa tal2eeglab tal2neuromag tal2itab tal2acpc tal2spm tal2mni tal2fsaverage tal2tal tal2scanras tal2scanlps tal2dicom tal2als tal2ali tal2ars tal2ari tal2pls tal2pli tal2prs tal2pri tal2las tal2lai tal2ras tal2rai tal2lps tal2lpi tal2rps tal2rpi tal2asl tal2ail tal2asr tal2air tal2psl tal2pil tal2psr tal2pir tal2sal tal2ial tal2sar tal2iar tal2spl tal2ipl tal2spr tal2ipr tal2sla tal2ila tal2sra tal2ira tal2slp tal2ilp tal2srp tal2irp tal2lsa tal2lia tal2rsa tal2ria tal2lsp tal2lip tal2rsp tal2rip
persistent scanras2ctf scanras2bti scanras2fourd scanras2yokogawa scanras2eeglab scanras2neuromag scanras2itab scanras2acpc scanras2spm scanras2mni scanras2fsaverage scanras2tal scanras2scanras scanras2scanlps scanras2dicom scanras2als scanras2ali scanras2ars scanras2ari scanras2pls scanras2pli scanras2prs scanras2pri scanras2las scanras2lai scanras2ras scanras2rai scanras2lps scanras2lpi scanras2rps scanras2rpi scanras2asl scanras2ail scanras2asr scanras2air scanras2psl scanras2pil scanras2psr scanras2pir scanras2sal scanras2ial scanras2sar scanras2iar scanras2spl scanras2ipl scanras2spr scanras2ipr scanras2sla scanras2ila scanras2sra scanras2ira scanras2slp scanras2ilp scanras2srp scanras2irp scanras2lsa scanras2lia scanras2rsa scanras2ria scanras2lsp scanras2lip scanras2rsp scanras2rip
persistent scanlps2ctf scanlps2bti scanlps2fourd scanlps2yokogawa scanlps2eeglab scanlps2neuromag scanlps2itab scanlps2acpc scanlps2spm scanlps2mni scanlps2fsaverage scanlps2tal scanlps2scanras scanlps2scanlps scanlps2dicom scanlps2als scanlps2ali scanlps2ars scanlps2ari scanlps2pls scanlps2pli scanlps2prs scanlps2pri scanlps2las scanlps2lai scanlps2ras scanlps2rai scanlps2lps scanlps2lpi scanlps2rps scanlps2rpi scanlps2asl scanlps2ail scanlps2asr scanlps2air scanlps2psl scanlps2pil scanlps2psr scanlps2pir scanlps2sal scanlps2ial scanlps2sar scanlps2iar scanlps2spl scanlps2ipl scanlps2spr scanlps2ipr scanlps2sla scanlps2ila scanlps2sra scanlps2ira scanlps2slp scanlps2ilp scanlps2srp scanlps2irp scanlps2lsa scanlps2lia scanlps2rsa scanlps2ria scanlps2lsp scanlps2lip scanlps2rsp scanlps2rip
persistent dicom2ctf dicom2bti dicom2fourd dicom2yokogawa dicom2eeglab dicom2neuromag dicom2itab dicom2acpc dicom2spm dicom2mni dicom2fsaverage dicom2tal dicom2scanras dicom2scanlps dicom2dicom dicom2als dicom2ali dicom2ars dicom2ari dicom2pls dicom2pli dicom2prs dicom2pri dicom2las dicom2lai dicom2ras dicom2rai dicom2lps dicom2lpi dicom2rps dicom2rpi dicom2asl dicom2ail dicom2asr dicom2air dicom2psl dicom2pil dicom2psr dicom2pir dicom2sal dicom2ial dicom2sar dicom2iar dicom2spl dicom2ipl dicom2spr dicom2ipr dicom2sla dicom2ila dicom2sra dicom2ira dicom2slp dicom2ilp dicom2srp dicom2irp dicom2lsa dicom2lia dicom2rsa dicom2ria dicom2lsp dicom2lip dicom2rsp dicom2rip
persistent als2ctf als2bti als2fourd als2yokogawa als2eeglab als2neuromag als2itab als2acpc als2spm als2mni als2fsaverage als2tal als2scanras als2scanlps als2dicom als2als als2ali als2ars als2ari als2pls als2pli als2prs als2pri als2las als2lai als2ras als2rai als2lps als2lpi als2rps als2rpi als2asl als2ail als2asr als2air als2psl als2pil als2psr als2pir als2sal als2ial als2sar als2iar als2spl als2ipl als2spr als2ipr als2sla als2ila als2sra als2ira als2slp als2ilp als2srp als2irp als2lsa als2lia als2rsa als2ria als2lsp als2lip als2rsp als2rip
persistent ali2ctf ali2bti ali2fourd ali2yokogawa ali2eeglab ali2neuromag ali2itab ali2acpc ali2spm ali2mni ali2fsaverage ali2tal ali2scanras ali2scanlps ali2dicom ali2als ali2ali ali2ars ali2ari ali2pls ali2pli ali2prs ali2pri ali2las ali2lai ali2ras ali2rai ali2lps ali2lpi ali2rps ali2rpi ali2asl ali2ail ali2asr ali2air ali2psl ali2pil ali2psr ali2pir ali2sal ali2ial ali2sar ali2iar ali2spl ali2ipl ali2spr ali2ipr ali2sla ali2ila ali2sra ali2ira ali2slp ali2ilp ali2srp ali2irp ali2lsa ali2lia ali2rsa ali2ria ali2lsp ali2lip ali2rsp ali2rip
persistent ars2ctf ars2bti ars2fourd ars2yokogawa ars2eeglab ars2neuromag ars2itab ars2acpc ars2spm ars2mni ars2fsaverage ars2tal ars2scanras ars2scanlps ars2dicom ars2als ars2ali ars2ars ars2ari ars2pls ars2pli ars2prs ars2pri ars2las ars2lai ars2ras ars2rai ars2lps ars2lpi ars2rps ars2rpi ars2asl ars2ail ars2asr ars2air ars2psl ars2pil ars2psr ars2pir ars2sal ars2ial ars2sar ars2iar ars2spl ars2ipl ars2spr ars2ipr ars2sla ars2ila ars2sra ars2ira ars2slp ars2ilp ars2srp ars2irp ars2lsa ars2lia ars2rsa ars2ria ars2lsp ars2lip ars2rsp ars2rip
persistent ari2ctf ari2bti ari2fourd ari2yokogawa ari2eeglab ari2neuromag ari2itab ari2acpc ari2spm ari2mni ari2fsaverage ari2tal ari2scanras ari2scanlps ari2dicom ari2als ari2ali ari2ars ari2ari ari2pls ari2pli ari2prs ari2pri ari2las ari2lai ari2ras ari2rai ari2lps ari2lpi ari2rps ari2rpi ari2asl ari2ail ari2asr ari2air ari2psl ari2pil ari2psr ari2pir ari2sal ari2ial ari2sar ari2iar ari2spl ari2ipl ari2spr ari2ipr ari2sla ari2ila ari2sra ari2ira ari2slp ari2ilp ari2srp ari2irp ari2lsa ari2lia ari2rsa ari2ria ari2lsp ari2lip ari2rsp ari2rip
persistent pls2ctf pls2bti pls2fourd pls2yokogawa pls2eeglab pls2neuromag pls2itab pls2acpc pls2spm pls2mni pls2fsaverage pls2tal pls2scanras pls2scanlps pls2dicom pls2als pls2ali pls2ars pls2ari pls2pls pls2pli pls2prs pls2pri pls2las pls2lai pls2ras pls2rai pls2lps pls2lpi pls2rps pls2rpi pls2asl pls2ail pls2asr pls2air pls2psl pls2pil pls2psr pls2pir pls2sal pls2ial pls2sar pls2iar pls2spl pls2ipl pls2spr pls2ipr pls2sla pls2ila pls2sra pls2ira pls2slp pls2ilp pls2srp pls2irp pls2lsa pls2lia pls2rsa pls2ria pls2lsp pls2lip pls2rsp pls2rip
persistent pli2ctf pli2bti pli2fourd pli2yokogawa pli2eeglab pli2neuromag pli2itab pli2acpc pli2spm pli2mni pli2fsaverage pli2tal pli2scanras pli2scanlps pli2dicom pli2als pli2ali pli2ars pli2ari pli2pls pli2pli pli2prs pli2pri pli2las pli2lai pli2ras pli2rai pli2lps pli2lpi pli2rps pli2rpi pli2asl pli2ail pli2asr pli2air pli2psl pli2pil pli2psr pli2pir pli2sal pli2ial pli2sar pli2iar pli2spl pli2ipl pli2spr pli2ipr pli2sla pli2ila pli2sra pli2ira pli2slp pli2ilp pli2srp pli2irp pli2lsa pli2lia pli2rsa pli2ria pli2lsp pli2lip pli2rsp pli2rip
persistent prs2ctf prs2bti prs2fourd prs2yokogawa prs2eeglab prs2neuromag prs2itab prs2acpc prs2spm prs2mni prs2fsaverage prs2tal prs2scanras prs2scanlps prs2dicom prs2als prs2ali prs2ars prs2ari prs2pls prs2pli prs2prs prs2pri prs2las prs2lai prs2ras prs2rai prs2lps prs2lpi prs2rps prs2rpi prs2asl prs2ail prs2asr prs2air prs2psl prs2pil prs2psr prs2pir prs2sal prs2ial prs2sar prs2iar prs2spl prs2ipl prs2spr prs2ipr prs2sla prs2ila prs2sra prs2ira prs2slp prs2ilp prs2srp prs2irp prs2lsa prs2lia prs2rsa prs2ria prs2lsp prs2lip prs2rsp prs2rip
persistent pri2ctf pri2bti pri2fourd pri2yokogawa pri2eeglab pri2neuromag pri2itab pri2acpc pri2spm pri2mni pri2fsaverage pri2tal pri2scanras pri2scanlps pri2dicom pri2als pri2ali pri2ars pri2ari pri2pls pri2pli pri2prs pri2pri pri2las pri2lai pri2ras pri2rai pri2lps pri2lpi pri2rps pri2rpi pri2asl pri2ail pri2asr pri2air pri2psl pri2pil pri2psr pri2pir pri2sal pri2ial pri2sar pri2iar pri2spl pri2ipl pri2spr pri2ipr pri2sla pri2ila pri2sra pri2ira pri2slp pri2ilp pri2srp pri2irp pri2lsa pri2lia pri2rsa pri2ria pri2lsp pri2lip pri2rsp pri2rip
persistent las2ctf las2bti las2fourd las2yokogawa las2eeglab las2neuromag las2itab las2acpc las2spm las2mni las2fsaverage las2tal las2scanras las2scanlps las2dicom las2als las2ali las2ars las2ari las2pls las2pli las2prs las2pri las2las las2lai las2ras las2rai las2lps las2lpi las2rps las2rpi las2asl las2ail las2asr las2air las2psl las2pil las2psr las2pir las2sal las2ial las2sar las2iar las2spl las2ipl las2spr las2ipr las2sla las2ila las2sra las2ira las2slp las2ilp las2srp las2irp las2lsa las2lia las2rsa las2ria las2lsp las2lip las2rsp las2rip
persistent lai2ctf lai2bti lai2fourd lai2yokogawa lai2eeglab lai2neuromag lai2itab lai2acpc lai2spm lai2mni lai2fsaverage lai2tal lai2scanras lai2scanlps lai2dicom lai2als lai2ali lai2ars lai2ari lai2pls lai2pli lai2prs lai2pri lai2las lai2lai lai2ras lai2rai lai2lps lai2lpi lai2rps lai2rpi lai2asl lai2ail lai2asr lai2air lai2psl lai2pil lai2psr lai2pir lai2sal lai2ial lai2sar lai2iar lai2spl lai2ipl lai2spr lai2ipr lai2sla lai2ila lai2sra lai2ira lai2slp lai2ilp lai2srp lai2irp lai2lsa lai2lia lai2rsa lai2ria lai2lsp lai2lip lai2rsp lai2rip
persistent ras2ctf ras2bti ras2fourd ras2yokogawa ras2eeglab ras2neuromag ras2itab ras2acpc ras2spm ras2mni ras2fsaverage ras2tal ras2scanras ras2scanlps ras2dicom ras2als ras2ali ras2ars ras2ari ras2pls ras2pli ras2prs ras2pri ras2las ras2lai ras2ras ras2rai ras2lps ras2lpi ras2rps ras2rpi ras2asl ras2ail ras2asr ras2air ras2psl ras2pil ras2psr ras2pir ras2sal ras2ial ras2sar ras2iar ras2spl ras2ipl ras2spr ras2ipr ras2sla ras2ila ras2sra ras2ira ras2slp ras2ilp ras2srp ras2irp ras2lsa ras2lia ras2rsa ras2ria ras2lsp ras2lip ras2rsp ras2rip
persistent rai2ctf rai2bti rai2fourd rai2yokogawa rai2eeglab rai2neuromag rai2itab rai2acpc rai2spm rai2mni rai2fsaverage rai2tal rai2scanras rai2scanlps rai2dicom rai2als rai2ali rai2ars rai2ari rai2pls rai2pli rai2prs rai2pri rai2las rai2lai rai2ras rai2rai rai2lps rai2lpi rai2rps rai2rpi rai2asl rai2ail rai2asr rai2air rai2psl rai2pil rai2psr rai2pir rai2sal rai2ial rai2sar rai2iar rai2spl rai2ipl rai2spr rai2ipr rai2sla rai2ila rai2sra rai2ira rai2slp rai2ilp rai2srp rai2irp rai2lsa rai2lia rai2rsa rai2ria rai2lsp rai2lip rai2rsp rai2rip
persistent lps2ctf lps2bti lps2fourd lps2yokogawa lps2eeglab lps2neuromag lps2itab lps2acpc lps2spm lps2mni lps2fsaverage lps2tal lps2scanras lps2scanlps lps2dicom lps2als lps2ali lps2ars lps2ari lps2pls lps2pli lps2prs lps2pri lps2las lps2lai lps2ras lps2rai lps2lps lps2lpi lps2rps lps2rpi lps2asl lps2ail lps2asr lps2air lps2psl lps2pil lps2psr lps2pir lps2sal lps2ial lps2sar lps2iar lps2spl lps2ipl lps2spr lps2ipr lps2sla lps2ila lps2sra lps2ira lps2slp lps2ilp lps2srp lps2irp lps2lsa lps2lia lps2rsa lps2ria lps2lsp lps2lip lps2rsp lps2rip
persistent lpi2ctf lpi2bti lpi2fourd lpi2yokogawa lpi2eeglab lpi2neuromag lpi2itab lpi2acpc lpi2spm lpi2mni lpi2fsaverage lpi2tal lpi2scanras lpi2scanlps lpi2dicom lpi2als lpi2ali lpi2ars lpi2ari lpi2pls lpi2pli lpi2prs lpi2pri lpi2las lpi2lai lpi2ras lpi2rai lpi2lps lpi2lpi lpi2rps lpi2rpi lpi2asl lpi2ail lpi2asr lpi2air lpi2psl lpi2pil lpi2psr lpi2pir lpi2sal lpi2ial lpi2sar lpi2iar lpi2spl lpi2ipl lpi2spr lpi2ipr lpi2sla lpi2ila lpi2sra lpi2ira lpi2slp lpi2ilp lpi2srp lpi2irp lpi2lsa lpi2lia lpi2rsa lpi2ria lpi2lsp lpi2lip lpi2rsp lpi2rip
persistent rps2ctf rps2bti rps2fourd rps2yokogawa rps2eeglab rps2neuromag rps2itab rps2acpc rps2spm rps2mni rps2fsaverage rps2tal rps2scanras rps2scanlps rps2dicom rps2als rps2ali rps2ars rps2ari rps2pls rps2pli rps2prs rps2pri rps2las rps2lai rps2ras rps2rai rps2lps rps2lpi rps2rps rps2rpi rps2asl rps2ail rps2asr rps2air rps2psl rps2pil rps2psr rps2pir rps2sal rps2ial rps2sar rps2iar rps2spl rps2ipl rps2spr rps2ipr rps2sla rps2ila rps2sra rps2ira rps2slp rps2ilp rps2srp rps2irp rps2lsa rps2lia rps2rsa rps2ria rps2lsp rps2lip rps2rsp rps2rip
persistent rpi2ctf rpi2bti rpi2fourd rpi2yokogawa rpi2eeglab rpi2neuromag rpi2itab rpi2acpc rpi2spm rpi2mni rpi2fsaverage rpi2tal rpi2scanras rpi2scanlps rpi2dicom rpi2als rpi2ali rpi2ars rpi2ari rpi2pls rpi2pli rpi2prs rpi2pri rpi2las rpi2lai rpi2ras rpi2rai rpi2lps rpi2lpi rpi2rps rpi2rpi rpi2asl rpi2ail rpi2asr rpi2air rpi2psl rpi2pil rpi2psr rpi2pir rpi2sal rpi2ial rpi2sar rpi2iar rpi2spl rpi2ipl rpi2spr rpi2ipr rpi2sla rpi2ila rpi2sra rpi2ira rpi2slp rpi2ilp rpi2srp rpi2irp rpi2lsa rpi2lia rpi2rsa rpi2ria rpi2lsp rpi2lip rpi2rsp rpi2rip
persistent asl2ctf asl2bti asl2fourd asl2yokogawa asl2eeglab asl2neuromag asl2itab asl2acpc asl2spm asl2mni asl2fsaverage asl2tal asl2scanras asl2scanlps asl2dicom asl2als asl2ali asl2ars asl2ari asl2pls asl2pli asl2prs asl2pri asl2las asl2lai asl2ras asl2rai asl2lps asl2lpi asl2rps asl2rpi asl2asl asl2ail asl2asr asl2air asl2psl asl2pil asl2psr asl2pir asl2sal asl2ial asl2sar asl2iar asl2spl asl2ipl asl2spr asl2ipr asl2sla asl2ila asl2sra asl2ira asl2slp asl2ilp asl2srp asl2irp asl2lsa asl2lia asl2rsa asl2ria asl2lsp asl2lip asl2rsp asl2rip
persistent ail2ctf ail2bti ail2fourd ail2yokogawa ail2eeglab ail2neuromag ail2itab ail2acpc ail2spm ail2mni ail2fsaverage ail2tal ail2scanras ail2scanlps ail2dicom ail2als ail2ali ail2ars ail2ari ail2pls ail2pli ail2prs ail2pri ail2las ail2lai ail2ras ail2rai ail2lps ail2lpi ail2rps ail2rpi ail2asl ail2ail ail2asr ail2air ail2psl ail2pil ail2psr ail2pir ail2sal ail2ial ail2sar ail2iar ail2spl ail2ipl ail2spr ail2ipr ail2sla ail2ila ail2sra ail2ira ail2slp ail2ilp ail2srp ail2irp ail2lsa ail2lia ail2rsa ail2ria ail2lsp ail2lip ail2rsp ail2rip
persistent asr2ctf asr2bti asr2fourd asr2yokogawa asr2eeglab asr2neuromag asr2itab asr2acpc asr2spm asr2mni asr2fsaverage asr2tal asr2scanras asr2scanlps asr2dicom asr2als asr2ali asr2ars asr2ari asr2pls asr2pli asr2prs asr2pri asr2las asr2lai asr2ras asr2rai asr2lps asr2lpi asr2rps asr2rpi asr2asl asr2ail asr2asr asr2air asr2psl asr2pil asr2psr asr2pir asr2sal asr2ial asr2sar asr2iar asr2spl asr2ipl asr2spr asr2ipr asr2sla asr2ila asr2sra asr2ira asr2slp asr2ilp asr2srp asr2irp asr2lsa asr2lia asr2rsa asr2ria asr2lsp asr2lip asr2rsp asr2rip
persistent air2ctf air2bti air2fourd air2yokogawa air2eeglab air2neuromag air2itab air2acpc air2spm air2mni air2fsaverage air2tal air2scanras air2scanlps air2dicom air2als air2ali air2ars air2ari air2pls air2pli air2prs air2pri air2las air2lai air2ras air2rai air2lps air2lpi air2rps air2rpi air2asl air2ail air2asr air2air air2psl air2pil air2psr air2pir air2sal air2ial air2sar air2iar air2spl air2ipl air2spr air2ipr air2sla air2ila air2sra air2ira air2slp air2ilp air2srp air2irp air2lsa air2lia air2rsa air2ria air2lsp air2lip air2rsp air2rip
persistent psl2ctf psl2bti psl2fourd psl2yokogawa psl2eeglab psl2neuromag psl2itab psl2acpc psl2spm psl2mni psl2fsaverage psl2tal psl2scanras psl2scanlps psl2dicom psl2als psl2ali psl2ars psl2ari psl2pls psl2pli psl2prs psl2pri psl2las psl2lai psl2ras psl2rai psl2lps psl2lpi psl2rps psl2rpi psl2asl psl2ail psl2asr psl2air psl2psl psl2pil psl2psr psl2pir psl2sal psl2ial psl2sar psl2iar psl2spl psl2ipl psl2spr psl2ipr psl2sla psl2ila psl2sra psl2ira psl2slp psl2ilp psl2srp psl2irp psl2lsa psl2lia psl2rsa psl2ria psl2lsp psl2lip psl2rsp psl2rip
persistent pil2ctf pil2bti pil2fourd pil2yokogawa pil2eeglab pil2neuromag pil2itab pil2acpc pil2spm pil2mni pil2fsaverage pil2tal pil2scanras pil2scanlps pil2dicom pil2als pil2ali pil2ars pil2ari pil2pls pil2pli pil2prs pil2pri pil2las pil2lai pil2ras pil2rai pil2lps pil2lpi pil2rps pil2rpi pil2asl pil2ail pil2asr pil2air pil2psl pil2pil pil2psr pil2pir pil2sal pil2ial pil2sar pil2iar pil2spl pil2ipl pil2spr pil2ipr pil2sla pil2ila pil2sra pil2ira pil2slp pil2ilp pil2srp pil2irp pil2lsa pil2lia pil2rsa pil2ria pil2lsp pil2lip pil2rsp pil2rip
persistent psr2ctf psr2bti psr2fourd psr2yokogawa psr2eeglab psr2neuromag psr2itab psr2acpc psr2spm psr2mni psr2fsaverage psr2tal psr2scanras psr2scanlps psr2dicom psr2als psr2ali psr2ars psr2ari psr2pls psr2pli psr2prs psr2pri psr2las psr2lai psr2ras psr2rai psr2lps psr2lpi psr2rps psr2rpi psr2asl psr2ail psr2asr psr2air psr2psl psr2pil psr2psr psr2pir psr2sal psr2ial psr2sar psr2iar psr2spl psr2ipl psr2spr psr2ipr psr2sla psr2ila psr2sra psr2ira psr2slp psr2ilp psr2srp psr2irp psr2lsa psr2lia psr2rsa psr2ria psr2lsp psr2lip psr2rsp psr2rip
persistent pir2ctf pir2bti pir2fourd pir2yokogawa pir2eeglab pir2neuromag pir2itab pir2acpc pir2spm pir2mni pir2fsaverage pir2tal pir2scanras pir2scanlps pir2dicom pir2als pir2ali pir2ars pir2ari pir2pls pir2pli pir2prs pir2pri pir2las pir2lai pir2ras pir2rai pir2lps pir2lpi pir2rps pir2rpi pir2asl pir2ail pir2asr pir2air pir2psl pir2pil pir2psr pir2pir pir2sal pir2ial pir2sar pir2iar pir2spl pir2ipl pir2spr pir2ipr pir2sla pir2ila pir2sra pir2ira pir2slp pir2ilp pir2srp pir2irp pir2lsa pir2lia pir2rsa pir2ria pir2lsp pir2lip pir2rsp pir2rip
persistent sal2ctf sal2bti sal2fourd sal2yokogawa sal2eeglab sal2neuromag sal2itab sal2acpc sal2spm sal2mni sal2fsaverage sal2tal sal2scanras sal2scanlps sal2dicom sal2als sal2ali sal2ars sal2ari sal2pls sal2pli sal2prs sal2pri sal2las sal2lai sal2ras sal2rai sal2lps sal2lpi sal2rps sal2rpi sal2asl sal2ail sal2asr sal2air sal2psl sal2pil sal2psr sal2pir sal2sal sal2ial sal2sar sal2iar sal2spl sal2ipl sal2spr sal2ipr sal2sla sal2ila sal2sra sal2ira sal2slp sal2ilp sal2srp sal2irp sal2lsa sal2lia sal2rsa sal2ria sal2lsp sal2lip sal2rsp sal2rip
persistent ial2ctf ial2bti ial2fourd ial2yokogawa ial2eeglab ial2neuromag ial2itab ial2acpc ial2spm ial2mni ial2fsaverage ial2tal ial2scanras ial2scanlps ial2dicom ial2als ial2ali ial2ars ial2ari ial2pls ial2pli ial2prs ial2pri ial2las ial2lai ial2ras ial2rai ial2lps ial2lpi ial2rps ial2rpi ial2asl ial2ail ial2asr ial2air ial2psl ial2pil ial2psr ial2pir ial2sal ial2ial ial2sar ial2iar ial2spl ial2ipl ial2spr ial2ipr ial2sla ial2ila ial2sra ial2ira ial2slp ial2ilp ial2srp ial2irp ial2lsa ial2lia ial2rsa ial2ria ial2lsp ial2lip ial2rsp ial2rip
persistent sar2ctf sar2bti sar2fourd sar2yokogawa sar2eeglab sar2neuromag sar2itab sar2acpc sar2spm sar2mni sar2fsaverage sar2tal sar2scanras sar2scanlps sar2dicom sar2als sar2ali sar2ars sar2ari sar2pls sar2pli sar2prs sar2pri sar2las sar2lai sar2ras sar2rai sar2lps sar2lpi sar2rps sar2rpi sar2asl sar2ail sar2asr sar2air sar2psl sar2pil sar2psr sar2pir sar2sal sar2ial sar2sar sar2iar sar2spl sar2ipl sar2spr sar2ipr sar2sla sar2ila sar2sra sar2ira sar2slp sar2ilp sar2srp sar2irp sar2lsa sar2lia sar2rsa sar2ria sar2lsp sar2lip sar2rsp sar2rip
persistent iar2ctf iar2bti iar2fourd iar2yokogawa iar2eeglab iar2neuromag iar2itab iar2acpc iar2spm iar2mni iar2fsaverage iar2tal iar2scanras iar2scanlps iar2dicom iar2als iar2ali iar2ars iar2ari iar2pls iar2pli iar2prs iar2pri iar2las iar2lai iar2ras iar2rai iar2lps iar2lpi iar2rps iar2rpi iar2asl iar2ail iar2asr iar2air iar2psl iar2pil iar2psr iar2pir iar2sal iar2ial iar2sar iar2iar iar2spl iar2ipl iar2spr iar2ipr iar2sla iar2ila iar2sra iar2ira iar2slp iar2ilp iar2srp iar2irp iar2lsa iar2lia iar2rsa iar2ria iar2lsp iar2lip iar2rsp iar2rip
persistent spl2ctf spl2bti spl2fourd spl2yokogawa spl2eeglab spl2neuromag spl2itab spl2acpc spl2spm spl2mni spl2fsaverage spl2tal spl2scanras spl2scanlps spl2dicom spl2als spl2ali spl2ars spl2ari spl2pls spl2pli spl2prs spl2pri spl2las spl2lai spl2ras spl2rai spl2lps spl2lpi spl2rps spl2rpi spl2asl spl2ail spl2asr spl2air spl2psl spl2pil spl2psr spl2pir spl2sal spl2ial spl2sar spl2iar spl2spl spl2ipl spl2spr spl2ipr spl2sla spl2ila spl2sra spl2ira spl2slp spl2ilp spl2srp spl2irp spl2lsa spl2lia spl2rsa spl2ria spl2lsp spl2lip spl2rsp spl2rip
persistent ipl2ctf ipl2bti ipl2fourd ipl2yokogawa ipl2eeglab ipl2neuromag ipl2itab ipl2acpc ipl2spm ipl2mni ipl2fsaverage ipl2tal ipl2scanras ipl2scanlps ipl2dicom ipl2als ipl2ali ipl2ars ipl2ari ipl2pls ipl2pli ipl2prs ipl2pri ipl2las ipl2lai ipl2ras ipl2rai ipl2lps ipl2lpi ipl2rps ipl2rpi ipl2asl ipl2ail ipl2asr ipl2air ipl2psl ipl2pil ipl2psr ipl2pir ipl2sal ipl2ial ipl2sar ipl2iar ipl2spl ipl2ipl ipl2spr ipl2ipr ipl2sla ipl2ila ipl2sra ipl2ira ipl2slp ipl2ilp ipl2srp ipl2irp ipl2lsa ipl2lia ipl2rsa ipl2ria ipl2lsp ipl2lip ipl2rsp ipl2rip
persistent spr2ctf spr2bti spr2fourd spr2yokogawa spr2eeglab spr2neuromag spr2itab spr2acpc spr2spm spr2mni spr2fsaverage spr2tal spr2scanras spr2scanlps spr2dicom spr2als spr2ali spr2ars spr2ari spr2pls spr2pli spr2prs spr2pri spr2las spr2lai spr2ras spr2rai spr2lps spr2lpi spr2rps spr2rpi spr2asl spr2ail spr2asr spr2air spr2psl spr2pil spr2psr spr2pir spr2sal spr2ial spr2sar spr2iar spr2spl spr2ipl spr2spr spr2ipr spr2sla spr2ila spr2sra spr2ira spr2slp spr2ilp spr2srp spr2irp spr2lsa spr2lia spr2rsa spr2ria spr2lsp spr2lip spr2rsp spr2rip
persistent ipr2ctf ipr2bti ipr2fourd ipr2yokogawa ipr2eeglab ipr2neuromag ipr2itab ipr2acpc ipr2spm ipr2mni ipr2fsaverage ipr2tal ipr2scanras ipr2scanlps ipr2dicom ipr2als ipr2ali ipr2ars ipr2ari ipr2pls ipr2pli ipr2prs ipr2pri ipr2las ipr2lai ipr2ras ipr2rai ipr2lps ipr2lpi ipr2rps ipr2rpi ipr2asl ipr2ail ipr2asr ipr2air ipr2psl ipr2pil ipr2psr ipr2pir ipr2sal ipr2ial ipr2sar ipr2iar ipr2spl ipr2ipl ipr2spr ipr2ipr ipr2sla ipr2ila ipr2sra ipr2ira ipr2slp ipr2ilp ipr2srp ipr2irp ipr2lsa ipr2lia ipr2rsa ipr2ria ipr2lsp ipr2lip ipr2rsp ipr2rip
persistent sla2ctf sla2bti sla2fourd sla2yokogawa sla2eeglab sla2neuromag sla2itab sla2acpc sla2spm sla2mni sla2fsaverage sla2tal sla2scanras sla2scanlps sla2dicom sla2als sla2ali sla2ars sla2ari sla2pls sla2pli sla2prs sla2pri sla2las sla2lai sla2ras sla2rai sla2lps sla2lpi sla2rps sla2rpi sla2asl sla2ail sla2asr sla2air sla2psl sla2pil sla2psr sla2pir sla2sal sla2ial sla2sar sla2iar sla2spl sla2ipl sla2spr sla2ipr sla2sla sla2ila sla2sra sla2ira sla2slp sla2ilp sla2srp sla2irp sla2lsa sla2lia sla2rsa sla2ria sla2lsp sla2lip sla2rsp sla2rip
persistent ila2ctf ila2bti ila2fourd ila2yokogawa ila2eeglab ila2neuromag ila2itab ila2acpc ila2spm ila2mni ila2fsaverage ila2tal ila2scanras ila2scanlps ila2dicom ila2als ila2ali ila2ars ila2ari ila2pls ila2pli ila2prs ila2pri ila2las ila2lai ila2ras ila2rai ila2lps ila2lpi ila2rps ila2rpi ila2asl ila2ail ila2asr ila2air ila2psl ila2pil ila2psr ila2pir ila2sal ila2ial ila2sar ila2iar ila2spl ila2ipl ila2spr ila2ipr ila2sla ila2ila ila2sra ila2ira ila2slp ila2ilp ila2srp ila2irp ila2lsa ila2lia ila2rsa ila2ria ila2lsp ila2lip ila2rsp ila2rip
persistent sra2ctf sra2bti sra2fourd sra2yokogawa sra2eeglab sra2neuromag sra2itab sra2acpc sra2spm sra2mni sra2fsaverage sra2tal sra2scanras sra2scanlps sra2dicom sra2als sra2ali sra2ars sra2ari sra2pls sra2pli sra2prs sra2pri sra2las sra2lai sra2ras sra2rai sra2lps sra2lpi sra2rps sra2rpi sra2asl sra2ail sra2asr sra2air sra2psl sra2pil sra2psr sra2pir sra2sal sra2ial sra2sar sra2iar sra2spl sra2ipl sra2spr sra2ipr sra2sla sra2ila sra2sra sra2ira sra2slp sra2ilp sra2srp sra2irp sra2lsa sra2lia sra2rsa sra2ria sra2lsp sra2lip sra2rsp sra2rip
persistent ira2ctf ira2bti ira2fourd ira2yokogawa ira2eeglab ira2neuromag ira2itab ira2acpc ira2spm ira2mni ira2fsaverage ira2tal ira2scanras ira2scanlps ira2dicom ira2als ira2ali ira2ars ira2ari ira2pls ira2pli ira2prs ira2pri ira2las ira2lai ira2ras ira2rai ira2lps ira2lpi ira2rps ira2rpi ira2asl ira2ail ira2asr ira2air ira2psl ira2pil ira2psr ira2pir ira2sal ira2ial ira2sar ira2iar ira2spl ira2ipl ira2spr ira2ipr ira2sla ira2ila ira2sra ira2ira ira2slp ira2ilp ira2srp ira2irp ira2lsa ira2lia ira2rsa ira2ria ira2lsp ira2lip ira2rsp ira2rip
persistent slp2ctf slp2bti slp2fourd slp2yokogawa slp2eeglab slp2neuromag slp2itab slp2acpc slp2spm slp2mni slp2fsaverage slp2tal slp2scanras slp2scanlps slp2dicom slp2als slp2ali slp2ars slp2ari slp2pls slp2pli slp2prs slp2pri slp2las slp2lai slp2ras slp2rai slp2lps slp2lpi slp2rps slp2rpi slp2asl slp2ail slp2asr slp2air slp2psl slp2pil slp2psr slp2pir slp2sal slp2ial slp2sar slp2iar slp2spl slp2ipl slp2spr slp2ipr slp2sla slp2ila slp2sra slp2ira slp2slp slp2ilp slp2srp slp2irp slp2lsa slp2lia slp2rsa slp2ria slp2lsp slp2lip slp2rsp slp2rip
persistent ilp2ctf ilp2bti ilp2fourd ilp2yokogawa ilp2eeglab ilp2neuromag ilp2itab ilp2acpc ilp2spm ilp2mni ilp2fsaverage ilp2tal ilp2scanras ilp2scanlps ilp2dicom ilp2als ilp2ali ilp2ars ilp2ari ilp2pls ilp2pli ilp2prs ilp2pri ilp2las ilp2lai ilp2ras ilp2rai ilp2lps ilp2lpi ilp2rps ilp2rpi ilp2asl ilp2ail ilp2asr ilp2air ilp2psl ilp2pil ilp2psr ilp2pir ilp2sal ilp2ial ilp2sar ilp2iar ilp2spl ilp2ipl ilp2spr ilp2ipr ilp2sla ilp2ila ilp2sra ilp2ira ilp2slp ilp2ilp ilp2srp ilp2irp ilp2lsa ilp2lia ilp2rsa ilp2ria ilp2lsp ilp2lip ilp2rsp ilp2rip
persistent srp2ctf srp2bti srp2fourd srp2yokogawa srp2eeglab srp2neuromag srp2itab srp2acpc srp2spm srp2mni srp2fsaverage srp2tal srp2scanras srp2scanlps srp2dicom srp2als srp2ali srp2ars srp2ari srp2pls srp2pli srp2prs srp2pri srp2las srp2lai srp2ras srp2rai srp2lps srp2lpi srp2rps srp2rpi srp2asl srp2ail srp2asr srp2air srp2psl srp2pil srp2psr srp2pir srp2sal srp2ial srp2sar srp2iar srp2spl srp2ipl srp2spr srp2ipr srp2sla srp2ila srp2sra srp2ira srp2slp srp2ilp srp2srp srp2irp srp2lsa srp2lia srp2rsa srp2ria srp2lsp srp2lip srp2rsp srp2rip
persistent irp2ctf irp2bti irp2fourd irp2yokogawa irp2eeglab irp2neuromag irp2itab irp2acpc irp2spm irp2mni irp2fsaverage irp2tal irp2scanras irp2scanlps irp2dicom irp2als irp2ali irp2ars irp2ari irp2pls irp2pli irp2prs irp2pri irp2las irp2lai irp2ras irp2rai irp2lps irp2lpi irp2rps irp2rpi irp2asl irp2ail irp2asr irp2air irp2psl irp2pil irp2psr irp2pir irp2sal irp2ial irp2sar irp2iar irp2spl irp2ipl irp2spr irp2ipr irp2sla irp2ila irp2sra irp2ira irp2slp irp2ilp irp2srp irp2irp irp2lsa irp2lia irp2rsa irp2ria irp2lsp irp2lip irp2rsp irp2rip
persistent lsa2ctf lsa2bti lsa2fourd lsa2yokogawa lsa2eeglab lsa2neuromag lsa2itab lsa2acpc lsa2spm lsa2mni lsa2fsaverage lsa2tal lsa2scanras lsa2scanlps lsa2dicom lsa2als lsa2ali lsa2ars lsa2ari lsa2pls lsa2pli lsa2prs lsa2pri lsa2las lsa2lai lsa2ras lsa2rai lsa2lps lsa2lpi lsa2rps lsa2rpi lsa2asl lsa2ail lsa2asr lsa2air lsa2psl lsa2pil lsa2psr lsa2pir lsa2sal lsa2ial lsa2sar lsa2iar lsa2spl lsa2ipl lsa2spr lsa2ipr lsa2sla lsa2ila lsa2sra lsa2ira lsa2slp lsa2ilp lsa2srp lsa2irp lsa2lsa lsa2lia lsa2rsa lsa2ria lsa2lsp lsa2lip lsa2rsp lsa2rip
persistent lia2ctf lia2bti lia2fourd lia2yokogawa lia2eeglab lia2neuromag lia2itab lia2acpc lia2spm lia2mni lia2fsaverage lia2tal lia2scanras lia2scanlps lia2dicom lia2als lia2ali lia2ars lia2ari lia2pls lia2pli lia2prs lia2pri lia2las lia2lai lia2ras lia2rai lia2lps lia2lpi lia2rps lia2rpi lia2asl lia2ail lia2asr lia2air lia2psl lia2pil lia2psr lia2pir lia2sal lia2ial lia2sar lia2iar lia2spl lia2ipl lia2spr lia2ipr lia2sla lia2ila lia2sra lia2ira lia2slp lia2ilp lia2srp lia2irp lia2lsa lia2lia lia2rsa lia2ria lia2lsp lia2lip lia2rsp lia2rip
persistent rsa2ctf rsa2bti rsa2fourd rsa2yokogawa rsa2eeglab rsa2neuromag rsa2itab rsa2acpc rsa2spm rsa2mni rsa2fsaverage rsa2tal rsa2scanras rsa2scanlps rsa2dicom rsa2als rsa2ali rsa2ars rsa2ari rsa2pls rsa2pli rsa2prs rsa2pri rsa2las rsa2lai rsa2ras rsa2rai rsa2lps rsa2lpi rsa2rps rsa2rpi rsa2asl rsa2ail rsa2asr rsa2air rsa2psl rsa2pil rsa2psr rsa2pir rsa2sal rsa2ial rsa2sar rsa2iar rsa2spl rsa2ipl rsa2spr rsa2ipr rsa2sla rsa2ila rsa2sra rsa2ira rsa2slp rsa2ilp rsa2srp rsa2irp rsa2lsa rsa2lia rsa2rsa rsa2ria rsa2lsp rsa2lip rsa2rsp rsa2rip
persistent ria2ctf ria2bti ria2fourd ria2yokogawa ria2eeglab ria2neuromag ria2itab ria2acpc ria2spm ria2mni ria2fsaverage ria2tal ria2scanras ria2scanlps ria2dicom ria2als ria2ali ria2ars ria2ari ria2pls ria2pli ria2prs ria2pri ria2las ria2lai ria2ras ria2rai ria2lps ria2lpi ria2rps ria2rpi ria2asl ria2ail ria2asr ria2air ria2psl ria2pil ria2psr ria2pir ria2sal ria2ial ria2sar ria2iar ria2spl ria2ipl ria2spr ria2ipr ria2sla ria2ila ria2sra ria2ira ria2slp ria2ilp ria2srp ria2irp ria2lsa ria2lia ria2rsa ria2ria ria2lsp ria2lip ria2rsp ria2rip
persistent lsp2ctf lsp2bti lsp2fourd lsp2yokogawa lsp2eeglab lsp2neuromag lsp2itab lsp2acpc lsp2spm lsp2mni lsp2fsaverage lsp2tal lsp2scanras lsp2scanlps lsp2dicom lsp2als lsp2ali lsp2ars lsp2ari lsp2pls lsp2pli lsp2prs lsp2pri lsp2las lsp2lai lsp2ras lsp2rai lsp2lps lsp2lpi lsp2rps lsp2rpi lsp2asl lsp2ail lsp2asr lsp2air lsp2psl lsp2pil lsp2psr lsp2pir lsp2sal lsp2ial lsp2sar lsp2iar lsp2spl lsp2ipl lsp2spr lsp2ipr lsp2sla lsp2ila lsp2sra lsp2ira lsp2slp lsp2ilp lsp2srp lsp2irp lsp2lsa lsp2lia lsp2rsa lsp2ria lsp2lsp lsp2lip lsp2rsp lsp2rip
persistent lip2ctf lip2bti lip2fourd lip2yokogawa lip2eeglab lip2neuromag lip2itab lip2acpc lip2spm lip2mni lip2fsaverage lip2tal lip2scanras lip2scanlps lip2dicom lip2als lip2ali lip2ars lip2ari lip2pls lip2pli lip2prs lip2pri lip2las lip2lai lip2ras lip2rai lip2lps lip2lpi lip2rps lip2rpi lip2asl lip2ail lip2asr lip2air lip2psl lip2pil lip2psr lip2pir lip2sal lip2ial lip2sar lip2iar lip2spl lip2ipl lip2spr lip2ipr lip2sla lip2ila lip2sra lip2ira lip2slp lip2ilp lip2srp lip2irp lip2lsa lip2lia lip2rsa lip2ria lip2lsp lip2lip lip2rsp lip2rip
persistent rsp2ctf rsp2bti rsp2fourd rsp2yokogawa rsp2eeglab rsp2neuromag rsp2itab rsp2acpc rsp2spm rsp2mni rsp2fsaverage rsp2tal rsp2scanras rsp2scanlps rsp2dicom rsp2als rsp2ali rsp2ars rsp2ari rsp2pls rsp2pli rsp2prs rsp2pri rsp2las rsp2lai rsp2ras rsp2rai rsp2lps rsp2lpi rsp2rps rsp2rpi rsp2asl rsp2ail rsp2asr rsp2air rsp2psl rsp2pil rsp2psr rsp2pir rsp2sal rsp2ial rsp2sar rsp2iar rsp2spl rsp2ipl rsp2spr rsp2ipr rsp2sla rsp2ila rsp2sra rsp2ira rsp2slp rsp2ilp rsp2srp rsp2irp rsp2lsa rsp2lia rsp2rsa rsp2ria rsp2lsp rsp2lip rsp2rsp rsp2rip
persistent rip2ctf rip2bti rip2fourd rip2yokogawa rip2eeglab rip2neuromag rip2itab rip2acpc rip2spm rip2mni rip2fsaverage rip2tal rip2scanras rip2scanlps rip2dicom rip2als rip2ali rip2ars rip2ari rip2pls rip2pli rip2prs rip2pri rip2las rip2lai rip2ras rip2rai rip2lps rip2lpi rip2rps rip2rpi rip2asl rip2ail rip2asr rip2air rip2psl rip2pil rip2psr rip2pir rip2sal rip2ial rip2sar rip2iar rip2spl rip2ipl rip2spr rip2ipr rip2sla rip2ila rip2sra rip2ira rip2slp rip2ilp rip2srp rip2irp rip2lsa rip2lia rip2rsa rip2ria rip2lsp rip2lip rip2rsp rip2rip

% get the options
requireorigin = ft_getopt(varargin, 'requireorigin', true);

% ensure these are in lower case
original = lower(original);
target = lower(target);

if isempty(initialized)
  % keep all variables persistent, these are only constructed once
  initialized = true;

  generic = {
    'als'; 'ali'; 'ars'; 'ari';...
    'pls'; 'pli'; 'prs'; 'pri';...
    'las'; 'lai'; 'ras'; 'rai';...
    'lps'; 'lpi'; 'rps'; 'rpi';...
    'asl'; 'ail'; 'asr'; 'air';...
    'psl'; 'pil'; 'psr'; 'pir';...
    'sal'; 'ial'; 'sar'; 'iar';...
    'spl'; 'ipl'; 'spr'; 'ipr';...
    'sla'; 'ila'; 'sra'; 'ira';...
    'slp'; 'ilp'; 'srp'; 'irp';...
    'lsa'; 'lia'; 'rsa'; 'ria';...
    'lsp'; 'lip'; 'rsp'; 'rip'}';

  % the specific ones specify the axes and also the origin
  specific = {'ctf', 'bti', 'fourd', 'yokogawa', 'eeglab', 'neuromag', 'itab', 'acpc', 'spm', 'mni', 'fsaverage', 'tal', 'scanras', 'scanlps', 'dicom'};

  if false
    % this section can be used to print all persistent variables that should be retained between calls
    for xxx=[specific generic]
      fprintf('persistent');
      for yyy=[specific generic]
        fprintf(' %s2%s', xxx{1}, yyy{1});
      end
      fprintf('\n');
    end
  end

  %--------------------------------------------------------------------------
  % These approximate alignments are based on transformation matrices that were
  % determined by clicking on the CTF fiducial locations in the canonical T1 template
  % MRI.
  %
  % All of the transformation matrices here are expressed in millimeter.

  % this is based on the ear canals, see ALIGN_CTF2ACPC
  acpc2ctf = [
    0.0000  0.9987  0.0517  34.7467
    -1.0000  0.0000  0.0000   0.0000
    0.0000 -0.0517  0.9987  52.2749
    0.0000  0.0000  0.0000   1.0000
    ];

  % this is based on the ear canals, see ALIGN_NEUROMAG2ACPC
  acpc2neuromag = [
    1.0000  0.0000  0.0000   0.0000
    0.0000  0.9987  0.0517  34.7467
    0.0000 -0.0517  0.9987  52.2749
    0.0000  0.0000  0.0000   1.0000
    ];

  % see http://freesurfer.net/fswiki/CoordinateSystems
  fsaverage2mni = [
    0.9975   -0.0073    0.0176   -0.0429
    0.0146    1.0009   -0.0024    1.5496
    -0.0130   -0.0093    0.9971    1.1840
    0.0000    0.0000    0.0000    1.0000
    ];

  % this is a 90 degree rotation around the z-axis
  ctf2neuromag = [
    0.0000   -1.0000    0.0000    0.0000
    1.0000    0.0000    0.0000    0.0000
    0.0000    0.0000    1.0000    0.0000
    0.0000    0.0000    0.0000    1.0000
    ];

  % these are all combinations of 90 degree rotations and/or flips along one of the axes
  for i=1:length(generic)
    for j=1:length(generic)
      xxx = generic{i};
      yyy = generic{j};
      eval(sprintf('%s2%s = transform_generic(''%s'', ''%s'');', xxx, yyy, xxx, yyy));
    end
  end

  % affine transformation from MNI to Talairach, see http://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach
  % the non-linear (i.e. piecewise linear) transform between MNI and Talairach are implemented elsewhere, see the functions MNI2TAL and TAL2MNI
  mni2tal = [
    0.8800    0.0000    0.0000   -0.8000
    0.0000    0.9700    0.0000   -3.3200
    0.0000    0.0500    0.8800   -0.4400
    0.0000    0.0000    0.0000    1.0000
    ];

  % the CTF and BTI coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
  ctf2bti = eye(4);

  % the CTF and EEGLAB coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
  ctf2eeglab = eye(4);

  % the NEUROMAG and ITAB coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
  neuromag2itab = eye(4);

  % BTI and 4D are different names for exactly the same system, see http://www.fieldtriptoolbox.org/faq/coordsys/
  bti2fourd = eye(4);

  % the YOKOGAWA system expresses positions relative to the dewar, not relative to the head
  % see https://www.fieldtriptoolbox.org/faq/coordsys/#details-of-the-yokogawa-coordinate-system
  yokogawa2als = eye(4);

  % the SPM and MNI coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/acpc/
  spm2mni = eye(4);

  % the SPM (aka MNI) and ACPC coordinate system are not the same, but similar enough, see http://www.fieldtriptoolbox.org/faq/acpc/
  spm2acpc = eye(4);
  mni2acpc = eye(4);

  % the CTF, BTI and 4D coordinate systems are all ALS coordinate systems
  % but the origin is poorly defined in ALS, hence converting from ALS to another is problematic
  ctf2als   = eye(4);
  bti2als   = eye(4);
  fourd2als = eye(4);

  % the NEUROMAG, ITAB, ACPC, MNI, SPM and FSAVERAGE coordinate systems are all RAS coordinate systems
  % but the origin is poorly defined in RAS, hence converting from RAS to another is problematic
  neuromag2ras  = eye(4);
  itab2ras      = eye(4);
  acpc2ras      = eye(4);
  mni2ras       = eye(4);
  spm2ras       = eye(4);
  fsaverage2ras = eye(4);
  tal2ras       = eye(4);

  % the SCANRAS coordinate system is RAS with the origin at the center opf the gradient coil
  scanras2ras     = eye(4);

  % the DICOM and SCANLPS coordinate system are the same, and rotated 180 degrees from SCANRAS
  dicom2scanlps   = eye(4);
  dicom2lps       = eye(4);
  scanlps2lps     = eye(4);
  scanlps2scanras = lps2ras; % this is a 180 degree rotation around the z-axis

  % make the combined and the inverse transformations where possible
  coordsys = [specific generic];
  implemented = zeros(length(coordsys)); % this is only for debugging
  for i=1:numel(coordsys)
    for j=1:numel(coordsys)
      xxx = coordsys{i};
      yyy = coordsys{j};

      if isequal(xxx, yyy)
        % construct the transformations on the diagonal
        eval(sprintf('%s2%s = eye(4);', xxx, yyy));
        implemented(i,j) = 1;
      elseif ~isempty(eval(sprintf('%s2%s', xxx, yyy)))
        % construct the inverse transformations
        eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
        implemented(i,j) = 2;
        implemented(j,i) = 2;
      elseif ismember(xxx, specific) && ismember(yyy, generic)
        % try to make the transformation (and inverse) with a two-step approach
        % since we go from specific to generic and thereby loose the origin information anyway, it is fine to use any intermediate step
        for k=1:numel(coordsys)
          zzz = coordsys{k};
          if ~isempty(eval(sprintf('%s2%s', xxx, zzz))) && ~isempty(eval(sprintf('%s2%s', zzz, yyy)))
            eval(sprintf('%s2%s = %s2%s * %s2%s;', xxx, yyy, zzz, yyy, xxx, zzz));
            eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
            implemented(i,j) = 3;
            implemented(j,i) = 3;
            break
          end
        end % for k
      elseif ismember(xxx, specific) && ismember(yyy, specific)
        % try to make the transformation (and inverse) with a two-step approach
        % do not use the generic orientation triplets (like RAS and ALS) as intermediate steps between two specific coordinate systems
        for k=1:numel(specific)
          zzz = specific{k};
          if ~isempty(eval(sprintf('%s2%s', zzz, yyy))) && ~isempty(eval(sprintf('%s2%s', xxx, zzz)))
            eval(sprintf('%s2%s = %s2%s * %s2%s;', xxx, yyy, zzz, yyy, xxx, zzz));
            eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
            implemented(i,j) = 3;
            implemented(j,i) = 3;
            break
          end
        end % for k
      end

    end % for j
  end % for i

  % these conversions should be done using FT_VOLUMENORMALISE, as they imply scaling
  clear acpc2spm acpc2mni acpc2fsaverage acpc2tal

  % converting to/from TAL is only possible for some specific template coordinate systems
  clear bti2tal ctf2tal fourd2tal itab2tal neuromag2tal
  clear tal2bti tal2ctf tal2fourd tal2itab tal2neuromag

  if requireorigin
    % the origin is poorly defined in generic orientation triplets (like RAS and ALS)
    % hence converting between them is problematic with regards to translations
    for i=1:length(generic)
      for j=1:length(specific)
        xxx = generic{i};
        yyy = specific{j};
        eval(sprintf('clear %s2%s', xxx, yyy));
      end
    end
  end

  if false
    % this is only for checking the coverage of all conversions
    for i=1:length(coordsys)
      for j=1:length(coordsys)
        xxx = coordsys{i};
        yyy = coordsys{j};
        if ~exist(sprintf('%s2%s', xxx, yyy), 'var') || isempty(eval(sprintf('%s2%s', xxx, yyy)))
          % update the previous list of implemented transformations, since some have been cleared
          implemented(i,j) = 0;
        end
      end
    end
    figure; imagesc(implemented); caxis([0 3]);
    xticklabels(coordsys); xticks(1:numel(coordsys));
    yticklabels(coordsys); yticks(1:numel(coordsys));
  end

end % if not initialized

if strcmp(original, '4d')
  xxx = 'fourd'; % '4d' is not a valid variable name
else
  xxx = original;
end

if strcmp(target, '4d')
  yyy = 'fourd'; % '4d' is not a valid variable name
else
  yyy = target;
end

if exist(sprintf('%s2%s', xxx, yyy), 'var') && ~isempty(eval(sprintf('%s2%s', xxx, yyy)))
  transform = eval(sprintf('%s2%s', xxx, yyy));
else
  ft_error('converting from %s to %s is not supported', original, target);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to construct generic transformations such as RAS2ALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = transform_generic(from, to)

ap_in  = find(to=='a'   | to=='p');
ap_out = find(from=='a' | from=='p');
lr_in  = find(to=='l'   | to=='r');
lr_out = find(from=='l' | from=='r');
si_in  = find(to=='s'   | to=='i');
si_out = find(from=='s' | from=='i');

% index the axis according to ap,lr,si
order_in  = [ap_in  lr_in  si_in];
order_out = [ap_out lr_out si_out];

% check whether one of the axis needs flipping
flip = 2.*(0.5-double(to(order_in)~=from(order_out)));

T = zeros(4);
for k = 1:3
  T(order_in(k),order_out(k)) = flip(k);
end
T(4,4) = 1;
