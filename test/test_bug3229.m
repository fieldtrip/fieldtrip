function test_bug3229

% WALLTIME 00:10:00
% MEM 1gb

elec = [];
elec.label   = {'1';'2';'3';'4'};
elec.elecpos = [1 1 1; 2 2 2; 3 3 3; 4 4 4];
elec.chanpos = elec.elecpos;
elec.tra     = eye(4);

%%

bipolar.labelold  = {'1', '2', '3', '4'};
bipolar.labelnew  = {'1', '2', '3'};
bipolar.tra       = [
  +1 -1  0  0
  0 +1 -1  0
  0  0 +1 -1
  ];

elec_bi = ft_apply_montage(elec, bipolar);

% channel names are the same, so keep the same position
assert(isequal(elec_bi.chanpos(1:3,:), elec.chanpos(1:3,:)));

assert(isequal(elec_bi.chanpos(1:3,:), elec.chanpos(1:3,:)));

%%

bipolar.labelold  = {'1',   '2',   '3',   '4'};
bipolar.labelnew  = {'1-2', '2-3', '3-4'};
bipolar.tra       = [
  +1 -1  0  0
   0 +1 -1  0
   0  0 +1 -1
  ];

elec_bi = ft_apply_montage(elec, bipolar);

% channel names are not same, so do not keep the same position
assert(~isequal(elec_bi.chanpos(1:3,:), elec.chanpos(1:3,:)));

% channel positions should be half-way
assert(isequal(elec_bi.chanpos(1:3,:), (elec.chanpos(1:3,:)+elec.chanpos(2:4,:))/2));


