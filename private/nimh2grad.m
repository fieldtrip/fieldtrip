function [grad] = nimh2grad(hdr);

% NIMH2GRAD constructs a gradiometer definition from the res4 header whish
% is read using the NIMH implementation of ctf_read_res4. The grad
% structure is compatible with FieldTrip and Robert Oostenveld's low-level
% forward and inverse routines.
%
% Use as
%   hdr  = ctf_read_res4(dataset);
%   grad = nimh2grad(hdr;
%
% See also CTF2GRAD, FIF2GRAD

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: nimh2grad.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.3  2007/03/07 08:37:55  roboos
% fixed typo
%
% Revision 1.2  2007/03/06 09:37:35  roboos
% small change in determining the MEG channels, thanks to Nicolas Robitaille
%
% Revision 1.1  2006/08/31 13:32:11  roboos
% moved from fieldtrip to fileio module
%
% Revision 1.1  2005/05/26 09:55:55  roboos
% new implementation to complement the NIMH ctf reading routines
%

% only work on the MEG channels
if isfield(hdr.sensor.index, 'meg')
  sel = hdr.sensor.index.meg;
else
  sel = hdr.sensor.index.meg_sens;
end

% start with an empty structure
grad.pnt = [];
grad.ori = [];
grad.tra = [];
grad.label = {};

for i=1:length(sel)
  pnt = hdr.sensor.info(sel(i)).location';
  ori = hdr.sensor.info(sel(i)).orientation';
  numcoils(i) = size(pnt,1);
  if size(ori,1)==1 && size(pnt,1)==1
    % one coil position with one orientation: magnetometer
    ori = ori;
  elseif size(ori,1)==1 && size(pnt,1)==2
    % two coil positions with one orientation: first order gradiometer
    % assume that the orientation of the upper coil is opposite to the lower coil
    ori = [ori; -ori];
  else
    error('do not know how to deal with higher order gradiometer hardware')
  end

  % add this channels coil positions and orientations
  grad.pnt = [grad.pnt; pnt];
  grad.ori = [grad.ori; ori];
  grad.label{i} = hdr.sensor.info(sel(i)).label;

  % determine the contribution of each coil to each channel's output signal
  if size(pnt,1)==1
    % one coil, assume that the orientation is correct, i.e. the weight is +1
    grad.tra(i,end+1) = 1;
  elseif size(pnt,1)==2
    % two coils, assume that the orientation for each coil is correct, i.e. the weights are +1 and +1
    grad.tra(i,end+1) = 1;
    grad.tra(i,end+1) = 1;
  else
    error('do not know how to deal with higher order gradiometer hardware')
  end
end

% prefer to have the labels in a column vector
grad.label = grad.label(:);

% reorder the coils, such that the bottom coils are at the first N
% locations and the top coils at the last N positions. This makes it
% easier to use a selection of the coils for topographic plotting
if all(numcoils==2)
  bot = 1:2:sum(numcoils);
  top = 2:2:sum(numcoils);
  grad.pnt = grad.pnt([bot top], :);
  grad.ori = grad.ori([bot top], :);
  grad.tra = grad.tra(:, [bot top]);
end
