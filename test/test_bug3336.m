function test_bug3336

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

% there are two versions: one with error() and one with ft_error()
% these should give the same result

% compare the two ways of giving an error

try
  subfunction1err
catch
  err1 = lasterror;
end

try
  subfunction2err
catch
  err2 = lasterror;
end

assert(~isequal(err1.message, err2.message))
assert(isequal(err1.identifier, err2.identifier))

% idem for warnings
% but here I cannot get access to the stack

subfunction1warn
[msg1, id1] = lastwarn;

subfunction2warn
[msg2, id2] = lastwarn;

assert(~isequal(msg1, msg2))
assert(isequal(id1, id2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERSION 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subfunction1err
subsubfunction1err

function subsubfunction1err
error('SOME:ERROR', 'something went wrong 1');

function subfunction1warn
subsubfunction1warn

function subsubfunction1warn
warning('SOME:WARNING', 'something went wrong 1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERSION 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subfunction2err
subsubfunction2err

function subsubfunction2err
ft_error('SOME:ERROR', 'something went wrong 2');

function subfunction2warn
subsubfunction2warn

function subsubfunction2warn
warning('SOME:WARNING', 'something went wrong 2');
