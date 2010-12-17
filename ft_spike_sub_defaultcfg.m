function [cfg] = ft_spikestation_sub_defaultcfg(cfg,defaults)

% check if the cfg is empty or not, if empty, we can skip the lengthy error checking
if isempty(cfg), % set all to defaults
  cfg           = defaults; 
  fields        = fieldnames(cfg);
  for iField = 1:length(fields)
    val = cfg.(fields{iField});
    cfg.(fields{iField}) =  val{1};
  end
else

  % check which fields are not present and get the default fields at once.
  fields        = fieldnames(cfg);
  defaultFields = fieldnames(defaults);

  % get the fields of defaultFields that are not in fields
  setDefaults = find(~ismember(defaultFields,fields));
  for k = setDefaults(:)'
    val = defaults.(defaultFields{k});
    cfg.(defaultFields{k}) = val{1}; % default is always first of the cell array
  end

  % issue errors for obsolete fields, to avoid typo mistakes
  unusedFields  = find(~ismember(fields,defaultFields));
  if ~isempty(unusedFields),
    for k = unusedFields(:)'
        sprintf('%s%s %s\n%s',...
        'Input field cfg.', fields{k}, 'is not used by this function.', ...
        'Check for typos and remove unused field')
    end
    error('MATLAB:spikestation:defaults:unusedFields',...
        'Make sure there are no unused fields')      
  end

  % issue errors when the input is not allowed
  for iField = 1:length(fields)
    vals  = defaults.(fields{iField});  
    nVals = length(vals);
    if nVals>1 && ~any(strcmp(cfg.(fields{iField}),vals)) 
      valLabel = cell(1:length(vals));
      for k = 1:length(vals)-1
        valLabel{k} = ['"' vals{k} '"' ', or '];
      end
      valLabel{end} = ['"' vals{end} '".'];
      error(['MATLAB:spikestation:defaults:cfg:' fields{iField} 'optionDoesNotExist'],...
      '%s%s %s %s', 'cfg.', fields{iField}, 'should be', valLabel{:})
    end
  end
end 