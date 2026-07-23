function elec = opto2elec(opto)

% OPTTO2ELEC converts an opto structure to an elec structure
%
% See also GRAD2ELEC

elec = keepfields(opto, {'label', 'unit', 'coordsys'});
elec.type = 'eeg';
if isfield(grad, 'tra')
  weight = abs(opto.tra);
  for i=1:numel(opto.label)
    weight(i,:) = weight(i,:) / norm(weight(i,:));
  end
  elec.elecpos = weight * opto.optopos;
else
  elec.elecpos = opto.optopos;
end