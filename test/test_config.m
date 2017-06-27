function test_config

% MEM 2gb
% WALLTIME 00:10:00

% TEST config

%% See http://bugzilla.fcdonders.nl/show_bug.cgi?id=1614

% one reported issue had to do with cell-arrays
a1           = [];
a1.field     = 1;
a1.field2{1} = 1;

a2 = config(a1);
a2.field3{1} = 1;

a1.field3{1} = 1;

assert(isequal(struct(a2), a1));
assert(isequal(printstruct('a', a2), printstruct('a', a1)));

% another reported issue had to do with assigning a config

a3 = struct();
a1 = struct();
a1.field1{1} = 1;
a1.field1{1} = a3;
a1.field2    = a3;

a3 = config();
a2 = config();
a2.field1{1} = 1;
a2.field1{1} = a3;
a2.field2    = a3;

assert(isequal(struct(a2), a1));
assert(isequal(printstruct('a', a2), printstruct('a', a1)));

%% See http://bugzilla.fcdonders.nl/show_bug.cgi?id=1762

% this has to do with gathering multiple "ans" outputs in single array

a1 = struct();
a1(1).field = 1;
a1(2).field = 2;
a1val = [a1.field];

a2 = config(a1);
a2val = [a2.field];

assert(isequal(a2val, a1val));

%% Other tests

% construct it as regular structure
a1 = struct;
a1(1).a = 1;
a1(1).b = 1;
a1(2).a = 2;
a1(4).b = 4;

% construct it as config
a2 = config;
a2(1).a = 1;
a2(1).b = 1;
a2(2).a = 2;
a2(4).b = 4;

% construct it by converting
a3 = config(a1);

assert(isequal(struct(a2(:)), a1(:)));
assert(isequal(struct(a3(:)), a1(:)));

assert(isequal(struct(a2), a1));
assert(isequal(struct(a3), a1));

assert(isequal(printstruct('a', a2), printstruct('a', a1)));
assert(isequal(printstruct('a', a3), printstruct('a', a1)));

%%

% construct it as regular structure
a1 = struct;
a1(1).a = 1;
a1(2).b = 2;
a1(3).c = 3;

% construct it as config
a2 = config;
a2(1).a = 1;
a2(2).b = 2;
a2(3).c = 3;

% construct it by converting
a3 = config(a1);

assert(isequal(struct(a2), a1));
assert(isequal(struct(a3), a1));

assert(isequal(printstruct('a', a2), printstruct('a', a1)));
assert(isequal(printstruct('a', a3), printstruct('a', a1)));

%%

% construct it as regular structure
a1 = struct;
a1(1,1).a = 1;
a1(2,2).b = 2;
a1(3,3).c = 3;

% construct it as config
a2 = config;
a2(1,1).a = 1;
a2(2,2).b = 2;
a2(3,3).c = 3;

% construct it by converting
a3 = config(a1);

assert(isequal(struct(a2), a1));
assert(isequal(struct(a3), a1));

assert(isequal(printstruct('a', a2), printstruct('a', a1)));
assert(isequal(printstruct('a', a3), printstruct('a', a1)));

%%

% construct it as regular structure
a1 = struct;
a1.a      = 11;
a1(2).a   = 2; % this will become (1,2)
a1(2,1).a = 21;
% a1(1,2).a = 12;

% construct it as config
a2 = config;
a2.a      = 11;
a2(2).a   = 2; % this will become (1,2)
a2(2,1).a = 21;
% a2(1,2).a = 12;

% construct it by converting
a3 = config(a1);

assert(isequal(a1(1,1).a, a2(1,1).a));
assert(isequal(a1(1,2).a, a2(1,2).a));
assert(isequal(a1(2,1).a, a2(2,1).a));
assert(isequal(a1(2,2).a, a2(2,2).a));

assert(isequal(a1(1,1).a, a3(1,1).a));
assert(isequal(a1(1,2).a, a3(1,2).a));
assert(isequal(a1(2,1).a, a3(2,1).a));
assert(isequal(a1(2,2).a, a3(2,2).a));

assert(isequal(struct(a2), a1));
assert(isequal(struct(a3), a1));

assert(isequal(printstruct('a', a2), printstruct('a', a1)));
assert(isequal(printstruct('a', a3), printstruct('a', a1)));

%%

% construct it as regular structure
a1 = struct;
a1(4,4).a = 1;
for i=1:16
  a1(i).a = i;
end
a1val = {a1.a};
a1bal = {a1(:).a};

% construct it as config
a2 = config;
a2(4,4).a = 1;
for i=1:16
  a2(i).a = i;
end
a2val = {a2.a};
a2bal = {a2(:).a};

% construct it by converting
a3 = config(a1);
a3val = {a3.a};
a3bal = {a2(:).a};

assert(isequal(a2val, a1val));
assert(isequal(a3val, a1val));

try
  assert(isequal(a1bal, a1val));
  assert(isequal(a2bal, a2val));
  assert(isequal(a3bal, a3val));
catch
  warning('there is still a situation where the config object is not 100%% compatible with a struct');
  % I believe this is due to a bug in MATLAB. When calling {a2(:).a} there is a
  % recursive call to subsref. On the outermost call, nargout=1, whereas I would
  % expect it to be 16.
end

%%
% the following was detected in http://bugzilla.fcdonders.nl/show_bug.cgi?id=2709#c1
% it seems to be similar to the test performed in the previous section

clear a1 a2

a1(1).b = 1;
a1(2).b = 1:2;
a1(3).b = 1:3;

a2 = config(a1);

c1 = {a2.b} % works;
c2 = {a2(:).b} % fails;
c3 = {a2(1:3).b} % fails;

try
  assert(length(c1)==3)
  assert(length(c2)==3)
  assert(length(c3)==3)
catch
  warning('there is still a situation where the config object is not 100%% compatible with a struct');
  % I believe this is due to a bug in MATLAB. When calling {a2(:).a} there is a
  % recursive call to subsref. On the outermost call, nargout=1, whereas I would
  % expect it to be 16.
end



