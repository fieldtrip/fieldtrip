function [type] = sensortype(grad)

% SENSORTYPE returns a string that describes the type of sensors (EEG or MEG)
% and the manufacturer of the MEG system. The heuristic approach is to test
% the input sensor definition on a few features (like number of channels,
% number of coils, etc.) and score points for each of them. The system that
% is most "similar" wins.
% 
% Use as
%   [str] = sensortype(sens)
% where the input is a electrode or gradiometer structure, or
%   [str] = sensortype(data)
% where the input data structure should contain either a data.grad
% or a data.elec field.
%
% The output will be a string
%   'electrode'
%   'ctf151'
%   'ctf275'
%   'ctf151_planar'
%   'ctf275_planar'
%   'neuromag122'
%   'neuromag306'
%   'magnetometer'
%   'bti148'
%   'yokogawa160'

% Copyright (C) 2004-2006, Robert Oostenveld
%
% $Log: sensortype.m,v $
% Revision 1.6  2007/12/12 10:39:45  roboos
% added try-end in case no pnt present, return 'unknown' as string instead of []
%
% Revision 1.5  2007/06/11 09:38:01  roboos
% added some code for yokogawa160, not yet complete, since yokogawa2grad is not yet implemented fully
%
% Revision 1.4  2007/06/11 09:17:34  roboos
% imporved detection for bti, better use of labels (start with 'A'), thanks to Nathan
%
% Revision 1.3  2007/05/06 09:08:20  roboos
% return warning instead of error when type is not detected
%
% Revision 1.2  2006/10/04 12:08:31  jansch
% changed variable name sens into grad consistently
%
% Revision 1.1  2006/10/04 08:00:29  roboos
% renamed megsystem to sensortype, added support for bti148 and EEG electrodes, renamed 'simulated magnetometer' into 'magnetometer'
%
% Revision 1.7  2006/08/31 08:03:24  roboos
% add explicit support for data as input (use data.grad), thanks to Floris
%
% Revision 1.6  2006/08/29 20:47:01  roboos
% added support for a simple simulated magnetometer system
%
% Revision 1.5  2006/01/30 14:06:04  roboos
% added square brackets around output variable in function definition
% added copyrights and log
% cleaned up help documentation
%

% the input may be a data structure which then contains a grad/elec structure
if isfield(grad, 'grad')
  grad = grad.grad;
elseif isfield(grad, 'elec')
  grad = grad.elec;
end

if isfield(grad, 'label') && isfield(grad, 'pnt') && ~isfield(grad, 'ori')
  % the input looks like an electrode system
  type = 'electrode';
  return;
end

% detect the different MEG sensor types
description = {
  'ctf151'                     % 1	    
  'ctf275'                     % 2
  'ctf151_planar'              % 3
  'ctf275_planar'              % 4
  'neuromag122'                % 5
  'neuromag306'                % 6
  'magnetometer'               % 7
  'bti148'                     % 8
  'yokogawa160'                % 9
};

% start with an empty counter for each system
similar = zeros(size(description));

% look at the number of channels
Nchan = length(grad.label);
similar(1) = similar(1) + (abs(Nchan-151)   <  20);
similar(2) = similar(2) + (abs(Nchan-275)   <  20);
similar(3) = similar(3) + (abs(Nchan-151*2) <  20);
similar(4) = similar(4) + (abs(Nchan-275*2) <  20);
similar(5) = similar(5) + (abs(Nchan-122)   <  20);
similar(6) = similar(6) + (abs(Nchan-306)   <  20);
similar(8) = similar(8) + (abs(Nchan-148)   < 20);
similar(9) = similar(9) + (abs(Nchan-160)   < 20);

similar(1) = similar(1) + (abs(Nchan-151)   <= abs(Nchan-275));
similar(2) = similar(2) + (abs(Nchan-275)   <= abs(Nchan-151));
similar(3) = similar(3) + (abs(Nchan-151*2) <= abs(Nchan-2*275));
similar(4) = similar(4) + (abs(Nchan-275*2) <= abs(Nchan-2*151));
similar(5) = similar(5) + (abs(Nchan-122)   <= abs(Nchan-306));
similar(6) = similar(6) + (abs(Nchan-306)   <= abs(Nchan-122));

% look at the beginning of BTi channel names
similar(8) = similar(8) + length(find(strmatch('A', grad.label)));
similar(8) = similar(8) - length(find(~strmatch('A', grad.label))); % penalty

% look at the beginning of Neuromag channel names
similar(5) = similar(5) + length(find(strmatch('MEG', grad.label)));
similar(6) = similar(6) + length(find(strmatch('MEG', grad.label)));
similar(7) = similar(7) + length(find(strmatch('MEG', grad.label)));
similar(5) = similar(5) - length(find(~strmatch('MEG', grad.label))); % penalty
similar(6) = similar(6) - length(find(~strmatch('MEG', grad.label))); % penalty
similar(7) = similar(7) - length(find(~strmatch('MEG', grad.label))); % penalty

% look at the beginning of CTF channel names
similar(1) = similar(1) + length(find(strmatch('MZ', grad.label)));
similar(2) = similar(2) + length(find(strmatch('MZ', grad.label)));
similar(3) = similar(3) + length(find(strmatch('MZ', grad.label)));
similar(4) = similar(4) + length(find(strmatch('MZ', grad.label)));
similar(5) = similar(5) - length(find(strmatch('MZ', grad.label))); % penalty
similar(6) = similar(6) - length(find(strmatch('MZ', grad.label))); % penalty
similar(7) = similar(7) - length(find(strmatch('MZ', grad.label)));  % penalty
similar(8) = similar(8) - length(find(strmatch('MZ', grad.label)));  % penalty

% look at whether the names contain _dH and _dV
similar(3) = similar(3) + length(find(a_strmatch('_dH', grad.label)));
similar(3) = similar(3) + length(find(a_strmatch('_dV', grad.label)));
similar(4) = similar(4) + length(find(a_strmatch('_dH', grad.label)));
similar(4) = similar(4) + length(find(a_strmatch('_dV', grad.label)));

% if they do not contain any _dH or _dV, the planar CTF systems is less likely that their axial counterpart
similar(3) = similar(3) - (length(find(a_strmatch('_dH', grad.label)))==0);
similar(3) = similar(3) - (length(find(a_strmatch('_dV', grad.label)))==0);
similar(4) = similar(4) - (length(find(a_strmatch('_dH', grad.label)))==0);
similar(4) = similar(4) - (length(find(a_strmatch('_dV', grad.label)))==0);

try
  % a magnetometer or a BTi system must have the same number of coils as labels
  similar(7) = similar(7) * (length(grad.label)==length(grad.pnt));
  similar(8) = similar(8) * (length(grad.label)==length(grad.pnt));
end

% determine to which MEG ssytem the input data is the most similar
[m, i] = max(similar);
if m==0 || length(find(similar==m))>1
  warning('could not determine the type of MEG system');
  type = 'unknown';
else
  type = description{i};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION search for substring in each element of a cell-array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = a_strmatch(str, strs)
a = zeros(length(strs),1);
for i=1:length(strs)
  a(i) = ~isempty(strfind(strs{i}, str));
end
a = find(a);

