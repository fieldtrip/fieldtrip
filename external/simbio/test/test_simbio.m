function test_simbio
% test file for SimBio

vol = [];
vol.labels   = [101 102 103 104];
vol.cond     = [1 0.022 0.33 0.33];
vol.headmesh = ['cubic_spheremodel_2mm.v'];

% read dipoles
fid   = fopen('source_space_N10048.dip', 'r');
s     = fgetl(fid);
[dum] = sscanf(s, '%c%d');
Ndip  = dum(2);
tline = fgetl(fid);

while isempty(findstr(tline, 'PositionsFixed'))
  tline = fgetl(fid);
end

dip = zeros(Ndip,3);
for i = 1:Ndip
  tline = fgetl(fid);
  dip(i,:) = sscanf(tline, '%f%f%f')';
end
fclose(fid);

% read elc
sensname = ['EEG_sensors.elc'];
elc      = ft_read_sens(sensname);
tmp      = str2num(elc.label{1});

for i = 1:size(elc.elecpos, 1)
  lbl{i} = num2str(tmp(i));
end
elc.label = lbl';

% call simbio
here = cd('../../../forward/private');     % Store current PWD, restore later.
try
  [lf] = leadfield_simbio(elc, vol, dip);  % This currently (SVN r6715) fails.
catch err
  cd(here);                                % Go back to current dir.
  rethrow(err);  
end
cd(here);                                  % Go back to current dir.
