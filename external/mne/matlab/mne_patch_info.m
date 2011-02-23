function [pinfo] = mne_patch_info(nearest)
%
% [pinfo] = mne_patch_info(nearest)
%
% Generate the patch information from the 'nearest' vector in a source space
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%

me='MNE:mne_patch_info';

if nargin ~= 1
   error(me,'Incorrect number of arguments');
end

if isempty(nearest)
   pinfo = [];
   return;
end

[ sorted, indn ] = sort(nearest);

[uniq,firsti,dum] = unique(sorted,'first');
[uniq,lasti,dum] = unique(sorted,'last');

for k = 1:length(uniq)
   pinfo{k} = indn(firsti(k):lasti(k));
end

return;

