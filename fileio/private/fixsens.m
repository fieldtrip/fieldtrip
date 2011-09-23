function [sens] = fixsens(sens)

% FIXSENS is a helper function to convert the old-style sensor description into 
% a new style sensor description.
%
% Old style (MEG): sens.pnt -> coil positions
%		   sens.ori -> coil orientations
% 	   sens.tra -> balancing matrix from coils to channels
%		   sens.label -> channel labels
%		   sens.balance -> additional structure containing info about the balancing
%          sens.unit
%
% New style (MEG): sens.coilpos -> coil positions
%		   sens.coilori -> coil orientations
% 	   sens.tra -> balancing matrix from coils to channels
%		   sens.label -> channel labels
%		   sens.balance -> additional structure containing info about the balancing
%		   sens.chanpos -> channel positions
%      sens.chanori -> 'orientation' of the channel: this is needed for
%          synthetic planar gradient computation
%      sens.unit
%
% Old style (EEG/ECoG): sens.pnt -> electrode positions
% 		 sens.tra -> balancing matrix from electrodes to channels
%		   sens.label -> channel labels
%          sens.unit
%
% New style (EEG/ECoG): sens.elecpos -> electrode positions
% 	   sens.tra -> balancing matrix from electrodes to channels
%		   sens.label -> channel labels
%		   sens.chanpos -> channel positions
%          sens.unit
%
% Use as [sens] = fixsens(sens)

if isfield(sens, 'ori')
  isgrad    = 1;
  doconvert = 1;
elseif (isfield(sens, 'coilori') && isfield(sens, 'coilpos') && isfield(sens, 'chanpos'))
  isgrad    = 1;
  doconvert = 0;
elseif (isfield(sens, 'elecpos') && isfield(sens, 'chanpos'))
  isgrad    = 0;
  doconvert = 0;
else
  isgrad    = 0;
  doconvert = 1;
end

if isgrad && doconvert
  % sensor description is a MEG sensor-array, containing oriented coils
  
  [chanpos, chanori, tmp] = channelposition(sens, 'channel', 'all');
  sens.coilori = sens.ori; sens = rmfield(sens, 'ori');
  sens.coilpos = sens.pnt; sens = rmfield(sens, 'pnt');
  sens.chanpos = chanpos;
  sens.chanori = chanori;
elseif doconvert
  % sensor description is something else, EEG/ECoG etc
  
  chanpos      = channelposition(sens, 'channel', 'all');
  sens.elecpos = chanpos; sens = rmfield(sens, 'pnt');
  sens.chanpos = chanpos;
end
