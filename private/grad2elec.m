function elec = grad2elec(grad)

% GRAD2ELEC converts a grad structure to an elec structure
%
% See also OPTO2ELEC

elec = keepfields(grad, {'unit', 'coordsys'});
elec.type = 'eeg';
elec.label = grad.label;
elec.elecpos = grad.coilpos;
