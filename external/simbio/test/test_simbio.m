% test file for SimBio
% CAUTION: this test script makes use of private functions, so you should
% first locate your FieldTrip private folder, in the module FORWARD:
%   FTprivate = '$myFTfolder/forward/private';
% and the testing branch, which has been checkout previously
%   FTtesting = 'myFTtestingbranch';
% ..this means for me:
FTprivate = '/home/coherence/crimic/fieldtrip-dev/forward/private/';
FTtesting = '/home/coherence/crimic/fieldtrip-dev/testing/';

vol = [];
vol.labels = [101 102 103 104];
vol.cond = [1 0.022 0.33 0.33];
vol.headmesh = [FTtesting 'external/simbio/cubic_spheremodel_2mm.v'];

% read dipoles
filename = [FTtesting 'external/simbio/source_space_N10048.dip'];
fid = fopen(filename,'r');
s  = fgetl(fid); [dum] = sscanf(s,'%c%d'); Ndip = dum(2);
tline = fgetl(fid);
while isempty(findstr(tline,'PositionsFixed'))
  tline=fgetl(fid);
end
dip = zeros(Ndip,3);
for i=1:Ndip
  tline=fgetl(fid);
  dip(i,:)= sscanf(tline,'%f%f%f')';
end
fclose(fid);

% read elc
sensname = [FTtesting 'external/simbio/EEG_sensors.elc'];
elc = ft_read_sens(sensname);
tmp = str2num(elc.label{1});
for i=1:size(elc.pnt,1)
  lbl{i} = num2str(tmp(i));
end
elc.label = lbl';

% call simbio
cd(FTprivate)
[lf] = leadfield_simbio(elc, vol, dip);
