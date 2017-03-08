function test_randomseed

% MEM 1500mb
% WALLTIME 00:10:00


% this is related to bug 1205
% call this for every possible MATLAB version
ft_defaults;

state1=randomseed([]);
x=rand;
randomseed(state1);
y=rand;
if x~=y, error('1'), end

state2=randomseed(state1);
z=rand;
randomseed(state2);
a=rand;
if z~=y, error('2'), end
if z~=a, error('3'), end

state1=randomseed(4);
x=rand;
state2=randomseed(state1);
y=rand;
if x~=y, error('4'), end

state1=randomseed(state2);
z=rand;
state2=randomseed(state1);
a=rand;
if z~=y, error('5'), end
if z~=a, error('6'), end


try
    state3=randomseed([5 6]);
    error('this should not work')
catch 
end

try
    state3=randomseed('no');
    error('this should not work')
catch 
end




