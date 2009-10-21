function elec = read_asa_elc(fn);

% READ_ASA_ELC reads electrodes from an ASA electrode file
% converting the units to mm

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_asa_elc.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.4  2008/11/14 07:21:45  roboos
% newer ASA versions write the labels in front of the positions, like "FPz: 10.4 1.3 -2"
% added support for this, thanks to Thomas Hartmann
% use strcmpi instead of strcmp(lower())
%
% Revision 1.3  2003/12/16 10:23:55  roberto
% attemt to make it slightly more robust for the electrode labels
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

Npnt = read_asa(fn, 'NumberPositions=', '%d');
Ndhk = read_asa(fn, 'NumberPolygons=', '%d');
Unit = read_asa(fn, 'UnitPosition', '%s', 1);
pnt  = read_asa(fn, 'Positions', '%f', Npnt, ':');
dhk  = read_asa(fn, 'Polygons', '%d', Ndhk);
lab  = read_asa(fn, 'Labels', '%s', Npnt);

if strcmpi(Unit,'mm')
  pnt = 1*pnt;
elseif strcmpi(Unit,'cm')
  pnt = 100*pnt;
elseif strcmpi(Unit,'m')
  pnt = 1000*pnt;
elseif ~isempty(Unit)
  error(sprintf('Unknown unit of distance for electrodes (%s)', Unit));
end

if length(lab)==1 && iscell(lab)
  % the electrode labels were on a single line
  % reformat the electrode labels into an appropriately sized cell array
  remainder = lab{1};
  lab = {};
  for i=1:size(pnt,1)
    [lab{i}, remainder] = strtok(remainder);
  end
end

elec.pnt = pnt;
elec.dhk = dhk+1;
elec.label = lab;
