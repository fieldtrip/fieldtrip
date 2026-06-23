function elec = opto2elec(opto)

% OPTTO2ELEC converts an opto structure to an elec structure
%
% See also GRAD2ELEC

elec = keepfields(opto, {'unit', 'coordsys'});
elec.type = 'eeg';
elec.label = opto.label;
elec.elecpos = opto.optopos;
