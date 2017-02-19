function failed_bug1919

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1919
% TEST ft_senslabel

system = {
  'ant128'
  'btiref'
  'bti148'
  'bti148_planar'
  'bti248'
  'bti248_planar'
  'ctfref'
  'ctfheadloc'
  'ctf64'
  'ctf151'
  'ctf151_planar'
  'ctf275'
  'ctf275_planar'
  'neuromag122'
  'neuromag122alt'
  'neuromag306'
  'neuromag306alt'
  'eeg1020'
  'eeg1010'
  'eeg1005'
  'ext1020'
  'biosemi64'
  'biosemi128'
  'biosemi256'
  'egi32'
  'egi64'
  'egi128'
  'egi256'
  'itab28'
  'itab153'
  'itab153_planar'
  'yokogawa9'
  'yokogawa64'
  'yokogawa64_planar'
  'yokogawa160'
  'yokogawa160_planar'
  'yokogawa440'
  'eeg'
  'electrode'
  };

for i=1:length(system)
  disp(system{i});
  % clear ft_senslabel
  
  label = ft_senslabel(system{i});
end

%% test the combination of planar channels
system = {
  'bti148_planar'
  'bti248_planar'
  'ctf151_planar'
  'ctf275_planar'
  'neuromag122'
  'neuromag122alt'
  'neuromag306'
  'neuromag306alt'
  'itab153_planar'
  'yokogawa64_planar'
  'yokogawa160_planar'
  };

for i=1:length(system)
  disp(system{i});
  % clear ft_senslabel
  
  label1 = ft_senslabel(system{i}, 'output', 'normal');         % all physical channels
  label2 = ft_senslabel(system{i}, 'output', 'planarcombined'); % the virtual planar channels plus the combination thereof
  
  if i~=7 && i~=8
    % all systems except neuromag306
    assert(size(label1,2)==2); % one column
    assert(size(label2,2)==3); % three columns
    assert(size(label2,1)==size(label1,1));
  else
    assert(size(label1,2)==3); % one column
    assert(size(label2,2)==3); % three columns
    assert(size(label2,1)==size(label1,1));
  end
end

