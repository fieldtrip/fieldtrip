function [sens] = fixsens(sens)

% FIXSENS is a helper function to convert the old-style sensor description into 
% a new style sensor description.
%
% Old style (MEG): sens.pnt -> coil positions
%		   sens.ori -> coil orientations
% 		   sens.tra -> balancing matrix from coils to channels
%		   sens.label -> channel labels
%		   sens.balance -> additional structure containing info about the balancing
%
% New style (MEG): sens.coilpos -> coil positions
%		   sens.coilori -> coil orientations
% 		   sens.tra -> balancing matrix from coils to channels
%		   sens.label -> channel labels
%		   sens.balance -> additional structure containing info about the balancing
%		   sens.chanpos -> channel positions
%
% Old style (EEG/ECoG): sens.pnt -> electrode positions
% 		   sens.tra -> balancing matrix from electrodes to channels
%		   sens.label -> channel labels
%
% New style (EEG/ECoG): sens.elecpos -> electrode positions
% 		   sens.tra -> balancing matrix from electrodes to channels
%		   sens.label -> channel labels
%		   sens.chanpos -> channel positions
%
% Use as [sens] = fixsens(sens)

if isfield(sens, 'ori')
  % sensor description is a MEG sensor-array, containing oriented coils
  
  chanpos      = channelposition(sens);
  sens.coilori = sens.ori; sens = rmfield(sens, 'ori');
  sens.coilpos = sens.pnt; sens = rmfield(sens, 'pnt');
  sens.chanpos = chanpos;  
else
  % sensor description is something else, EEG/ECoG etc
  
  chanpos      = channelposition(sens);
  sens.elecpos = sens.pnt; sens = rmfield(sens, 'pnt');
end
