function test_ft_appendsens

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY


% regular
elec1.elecpos = [1 1 1; 2 2 2; 3 3 3];
elec1.chanpos = elec1.elecpos;
elec1.label = {'1';'2';'3'};
elec1.unit  = 'mm';
elec1.coordsys  = 'acpc';

elec2.elecpos = [4 4 4; 5 5 5];
elec2.chanpos = elec2.elecpos;
elec2.label = {'4';'5'};
elec2.unit  = 'mm';
elec2.coordsys  = 'acpc';

elec = ft_appendsens([], elec1, elec2);

% duplicate chanpos, elecpos, and channel
elec1.elecpos = [1 1 1; 2 2 2; 3 3 3];
elec1.chanpos = elec1.elecpos;
elec1.label = {'1';'2';'3'};
elec1.unit  = 'mm';
elec1.coordsys  = 'acpc';
elec1.tra  = [1 0 0; 0 1 0; 0 0 1];
elec2.elecpos = [2 2 2; 0 0 0]; % 2 2 2 is a duplicate
elec2.chanpos = elec2.elecpos;
elec2.label = {'2';'7'}; % 2 is a duplicate
elec2.unit  = 'mm';
elec2.coordsys  = 'acpc';
elec2.tra  = [1 0; 0 1];
elec = ft_appendsens([], elec1, elec2);

% duplicate elecpos (eg same elec used for two bipolar derivations)
elec1.elecpos = [1 1 1; 2 2 2];
elec1.chanpos = [1.5 1.5 1.5];
elec1.label = {'1-2'};
elec1.unit  = 'mm';
elec1.coordsys  = 'acpc';
elec1.tra  = [1 -1];
elec2.elecpos = [2 2 2; 3 3 3]; % 2 2 2 is a duplicate
elec2.chanpos = [2.5 2.5 2.5];
elec2.label = {'2-3'}; % 2 is a duplicate
elec2.unit  = 'mm';
elec2.coordsys  = 'acpc';
elec2.tra  = [1 -1];

elec = ft_appendsens([], elec1, elec2);

% duplicate chanpos and elecpos, but no duplicate label (eg two electrodes
% localized to the same location)
elec1.elecpos = [1 1 1; 2 2 2; 3 3 3];
elec1.chanpos = elec1.elecpos;
elec1.label = {'1';'2';'3'};
elec1.unit  = 'mm';
elec1.coordsys  = 'acpc';
elec1.tra  = [1 0 0; 0 1 0; 0 0 1];
elec2.elecpos = [2 2 2; 0 0 0]; % 2 2 2 is a duplicate
elec2.chanpos = elec2.elecpos;
elec2.label = {'8';'7'}; % no duplicate label
elec2.unit  = 'mm';
elec2.coordsys  = 'acpc';
elec2.tra  = [1 0; 0 1];

try
  elec = ft_appendsens([], elec1, elec2);
  % this SHOULD throw an error because two labels for the same chanpos are not allowed
  error('an error was expected but not thrown')
catch
  fprintf('ft_appendsens threw an error as expected\n')
end

