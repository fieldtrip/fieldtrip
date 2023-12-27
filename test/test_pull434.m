function test_pull434

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_appendsens
% DATA no

% this script tests the appending of sensor structures
% and specifically the handling of overlapping channels and/or electrodes
% see https://github.com/fieldtrip/fieldtrip/pull/434
% and the preceding https://github.com/fieldtrip/fieldtrip/pull/428

%% construct some electrode structures

elec1.label = {'1'; '2'; '3'};
elec1.elecpos = [1 1 1; 2 2 2; 3 3 3];
elec1.chanpos = [1 1 1; 2 2 2; 3 3 3];
elec1.tra = eye(3);

elec2.label = {'3'; '4'};
elec2.elecpos = [3 3 3; 4 4 4];
elec2.chanpos = [3 3 3; 4 4 4];
elec2.tra = eye(2);

elec3.label = {'4'; '5'; '6'};
elec3.elecpos = [4 4 4; 5 5 5; 6 6 6];
elec3.chanpos = [4 4 4; 5 5 5; 6 6 6];
elec3.tra = eye(3);

%% append them in different ways

cfg = [];

append1   = ft_appendsens(cfg, elec1);
append11  = ft_appendsens(cfg, elec1, elec1);
append111 = ft_appendsens(cfg, elec1, elec1, elec1);

assert(isequal(append1.label,     elec1.label));
assert(isequal(append1.chanpos,   elec1.chanpos));
assert(isequal(append1.tra,       elec1.tra));
assert(isequal(append11.label,    elec1.label));
assert(isequal(append11.chanpos,  elec1.chanpos));
assert(isequal(append11.tra,      elec1.tra));
assert(isequal(append111.label,   elec1.label));
assert(isequal(append111.chanpos, elec1.chanpos));
assert(isequal(append111.tra,     elec1.tra));

%%

append13  = ft_appendsens(cfg, elec1, elec3);
append12  = ft_appendsens(cfg, elec1, elec2);
append23  = ft_appendsens(cfg, elec2, elec3);

assert(numel(append13.label)==6);
assert(numel(append12.label)==4);
assert(numel(append23.label)==4);

assert(size(append13.chanpos,1)==6);
assert(size(append12.chanpos,1)==4);
assert(size(append23.chanpos,1)==4);

assert(size(append13.elecpos,1)==6);
assert(size(append12.elecpos,1)==4);
assert(size(append23.elecpos,1)==4);

assert(isequal(size(append13.tra), [6 6]));
assert(isequal(size(append12.tra), [4 4]));
assert(isequal(size(append23.tra), [4 4]));

%%

append123 = ft_appendsens(cfg, elec1, elec2, elec3);

assert(numel(append123.label)==6);
assert(size(append123.chanpos,1)==6);
assert(size(append123.elecpos,1)==6);
assert(isequal(size(append123.tra), [6 6]));

%%

elec1b.label = {'1b', '2b', '3b'};
elec1b.elecpos = [1 1 1; 2 2 2; 3 3 3];      % same elecpos
elec1b.chanpos = [1 1 1; 2 2 2; 3 3 3] + 10; % different chanpos
elec1b.tra = eye(3);

elec1c.label = {'1c', '2c', '3c'};
elec1c.elecpos = [1 1 1; 2 2 2; 3 3 3] + 10; % different elecpos
elec1c.chanpos = [1 1 1; 2 2 2; 3 3 3];      % same chanpos
elec1c.tra = eye(3);

%%

append11b = ft_appendsens(cfg, elec1, elec1b); % same elecpos, different label

assert(numel(append11b.label)==6);
assert(size(append11b.chanpos,1)==6);
assert(size(append11b.elecpos,1)==3);
assert(isequal(size(append11b.tra), [6 3]));

%%

append11c = ft_appendsens(cfg, elec1, elec1c); % different elecpos

assert(numel(append11c.label)==6);
assert(size(append11c.chanpos,1)==6);
assert(size(append11c.elecpos,1)==6);
assert(isequal(size(append11c.tra), [6 6]));

%%

elec1bip.label = {'1-2'; '2-3'};
elec1bip.elecpos = [1 1 1; 2 2 2; 3 3 3];
elec1bip.chanpos = [1.5 1.5 1.5; 2.5 2.5 2.5];
elec1bip.tra = [
  1 -1  0
  0  1 -1
  ];

append11bip = ft_appendsens(cfg, elec1, elec1bip); % bipolar version of same electrodes

assert(numel(append11bip.label)==5);
assert(size(append11bip.chanpos,1)==5);
assert(size(append11bip.elecpos,1)==3);
assert(isequal(size(append11bip.tra), [5 3]));

%%

grad1x.label = {'1_bx', '2_bx', '3_bx'}';
grad1x.coilpos = [1 1 1; 2 2 2; 3 3 3];
grad1x.coilori = [1 0 0; 1 0 0; 1 0 0];
grad1x.tra = eye(3);

grad1y.label = {'1_by', '2_by', '3_by'}';
grad1y.coilpos = [1 1 1; 2 2 2; 3 3 3];
grad1y.coilori = [0 1 0; 0 1 0; 0 1 0];
grad1y.tra = eye(3);

grad1z.label = {'1_bz', '2_bz', '3_bz'}';
grad1z.coilpos = [1 1 1; 2 2 2; 3 3 3];
grad1z.coilori = [0 0 1; 0 0 1; 0 0 1];
grad1z.tra = eye(3);

%%

append1xx = ft_appendsens(cfg, grad1x, grad1x);
assert(numel(append1xx.label)==3);
assert(size(append1xx.coilpos,1)==3);
assert(size(append1xx.coilori,1)==3);
assert(isequal(size(append1xx.tra), [3 3]));

%%

append1xy = ft_appendsens(cfg, grad1x, grad1y);
assert(numel(append1xy.label)==6);
assert(size(append1xy.coilpos,1)==6);
assert(size(append1xy.coilori,1)==6);
assert(isequal(size(append1xy.tra), [6 6]));

%%

append1xyz = ft_appendsens(cfg, grad1x, grad1y, grad1z);
assert(numel(append1xyz.label)==9);
assert(size(append1xyz.coilpos,1)==9);
assert(size(append1xyz.coilori,1)==9);
assert(isequal(size(append1xyz.tra), [9 9]));
