function failed_bug2377b

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_senslabel ft_senstype ft_chantype ft_chanunit ft_datatype_sens ft_apply_transform

[pnt, tri] = icosahedron162;

pnt = pnt .* 10;                % convert to cm
sel = find(pnt(:,3)>0);         % take the upper hemisphere
nchan = length(sel);            % there are 71 channels remaining

sens = [];
sens.elecpos = pnt(sel,:);
sens.unit = 'cm';

lab = ft_senslabel('eeg1010');  % take the channel names from this set

% perform the test sequence for a sensor array with 10-20 channel labels
sens.label = lab(1:nchan);
perform_actual_test(sens)

for i=1:nchan
  lab{i} = sprintf('ch%02d', i);
end

% perform the test sequence for a sensor array with unkown labels
sens.label = lab(1:nchan);
perform_actual_test(sens)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function perform_actual_test(sens)

nchan = length(sens.label);

% without most new fields
sens0 = sens;

% with manual chantype
sens.chantype = repmat({'eeg'}, nchan, 1);
sens1 = sens;

% with automatic chantype
sens.chantype = ft_chantype(sens);
sens2 = sens;

% with manual chanunit
sens.chanunit = repmat({'V'}, nchan, 1);
sens3 = sens;

% with automatic chanunit
sens.chanunit = ft_chanunit(sens);
sens4 = sens;

% with a tra
sens.tra = eye(numel(sens.label));
sens5 = sens;

%%

sens0a = ft_datatype_sens(sens0, 'version', 'upcoming');
sens1a = ft_datatype_sens(sens1, 'version', 'upcoming');
sens2a = ft_datatype_sens(sens2, 'version', 'upcoming');
sens3a = ft_datatype_sens(sens3, 'version', 'upcoming');
sens4a = ft_datatype_sens(sens4, 'version', 'upcoming');
sens5a = ft_datatype_sens(sens5, 'version', 'upcoming');

% not all of them have the tra field
assert(isequal(tryrmfield(sens0a, 'tra'), tryrmfield(sens1a, 'tra')));
assert(isequal(tryrmfield(sens0a, 'tra'), tryrmfield(sens2a, 'tra')));
assert(isequal(tryrmfield(sens0a, 'tra'), tryrmfield(sens3a, 'tra')));
assert(isequal(tryrmfield(sens0a, 'tra'), tryrmfield(sens4a, 'tra')));
assert(isequal(tryrmfield(sens0a, 'tra'), tryrmfield(sens5a, 'tra')));

%%

montage = [];
montage.labelold = sens.label;
montage.labelnew = montage.labelold;
montage.tra = detrend(eye(nchan), 'constant');

sens0b = ft_apply_montage(sens0, montage);
sens1b = ft_apply_montage(sens1, montage);
sens2b = ft_apply_montage(sens2, montage);
sens3b = ft_apply_montage(sens3, montage);
sens4b = ft_apply_montage(sens4, montage);
sens5b = ft_apply_montage(sens5, montage);

%%

sens0c = ft_datatype_sens(sens0b, 'version', 'upcoming', 'amplitude', 'uV', 'distance', 'mm');
sens1c = ft_datatype_sens(sens1b, 'version', 'upcoming', 'amplitude', 'uV', 'distance', 'mm');
sens2c = ft_datatype_sens(sens2b, 'version', 'upcoming', 'amplitude', 'uV', 'distance', 'mm');
sens3c = ft_datatype_sens(sens3b, 'version', 'upcoming', 'amplitude', 'uV', 'distance', 'mm');
sens4c = ft_datatype_sens(sens4b, 'version', 'upcoming', 'amplitude', 'uV', 'distance', 'mm');
sens5c = ft_datatype_sens(sens5b, 'version', 'upcoming', 'amplitude', 'uV', 'distance', 'mm');

% they should have all fields by now
assert(isequal(sens0c, sens1c));
assert(isequal(sens0c, sens2c));
assert(isequal(sens0c, sens3c));
assert(isequal(sens0c, sens4c));
assert(isequal(sens0c, sens5c));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = tryrmfield(s, f)
if isfield(s, f)
  s = rmfield(s, f);
end

