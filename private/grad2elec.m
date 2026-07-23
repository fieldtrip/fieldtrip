function elec = grad2elec(grad)

% GRAD2ELEC converts a grad structure to an elec structure
%
% See also OPTO2ELEC

elec = keepfields(grad, {'label', 'unit', 'coordsys'});
elec.type = 'eeg';
if isfield(grad, 'tra')
  weight = abs(grad.tra);
  for i=1:numel(grad.label)
    weight(i,:) = weight(i,:) / norm(weight(i,:));
  end
  elec.elecpos = weight * grad.coilpos;
else
  elec.elecpos = grad.coilpos;
end