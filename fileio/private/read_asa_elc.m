function elec = read_asa_elc(fn);

% READ_ASA_ELC reads electrodes from an ASA electrode file
% converting the units to mm

% Copyright (C) 2002, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
