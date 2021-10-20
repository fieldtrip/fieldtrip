function test_ft_senslabel

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_senslabel

if nargout
  % assume that this is called by RUNTESTS
  tests = functiontests(localfunctions);
else
  % assume that this is called from the command line
  func = localfunctions;
  for i=1:numel(func)
    fprintf('evaluating %s\n', func2str(func{i}));
    feval(func{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testSenslabel(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = {'ant128', 'biosemi128', 'biosemi256', 'biosemi64', 'bti148', 'bti148_planar', 'bti248', 'bti248_planar', 'btiref', 'ctf151', 'ctf151_planar', 'ctf275', 'ctf275_planar', 'ctf64', 'ctfheadloc', 'ctfref', 'eeg1005', 'eeg1010', 'eeg1020', 'egi128', 'egi256', 'egi32', 'egi64', 'ext1020', 'itab153', 'itab153_planar', 'itab28', 'yokogawa160', 'yokogawa160_planar', 'yokogawa440', 'yokogawa440_planar', 'yokogawa64', 'yokogawa64_planar', 'yokogawa9', 'neuromag122', 'neuromag122alt', 'neuromag122_combined','neuromag122alt_combined', 'neuromag306', 'neuromag306alt', 'neuromag306_combined', 'neuromag306alt_combined'};

for i=1:numel(type)
  disp(type{i})
  assert(~isempty(ft_senslabel(type{i})));
end