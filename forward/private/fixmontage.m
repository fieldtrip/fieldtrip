function montage = fixmontage(montage)

% FIXMONTAGE use "old/new" instead of "org/new"

if isfield(montage, 'labelorg')
  montage.labelold = montage.labelorg;
  montage = rmfield(montage, 'labelorg');
end
if isfield(montage, 'chantypeorg')
  montage.chantypeold = montage.chantypeorg;
  montage = rmfield(montage, 'chantypeorg');
end
if isfield(montage, 'chanunitorg')
  montage.chanunitold = montage.chanunitorg;
  montage = rmfield(montage, 'chanunitorg');
end
