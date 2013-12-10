function test_bug2418

% WALLTIME 00:05:00
% MEM 1500mb

% TEST test_bug2418
% TEST ft_senstype ft_senslabel

type = {
  'ant128'
  'biosemi64'
  'biosemi128'
  'biosemi256'
  'bti148'
  'bti148_planar'
  'bti248'
  'bti248_planar'
%  'btiref'         % returns a subset
  'ctf151'
  'ctf151_planar'
  'ctf275'
  'ctf275_planar'
%  'ctfheadloc'     % returns a subset
%  'ctfref'         % returns a subset
  'eeg1005'
  'eeg1010'
  'eeg1020'
%  'ext1020'        % returns lower and upper case
  'egi32'
  'egi64'
  'egi128'
  'egi256'
  'neuromag122'
  'neuromag122alt'
  'neuromag306'
  'neuromag306alt'
  'itab28'
  'itab153'
  'itab153_planar'
  'yokogawa9'
  'yokogawa64'
  'yokogawa64_planar'
  'yokogawa160'
  'yokogawa160_planar'
  'yokogawa440'
  'yokogawa440_planar'
  };


for i=1:length(type)
  lab = ft_senslabel(type{i});
  lab = lab(:);
  
  disp(type{i});
  assert(strcmp(type{i}, ft_senstype(lab)));
end