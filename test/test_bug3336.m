function test_bug3336

% there are two versions: one with error() and one with ft_error()
% these should give the same result

try
  subfunction1
catch
  err1 = lasterror;
end

try
  subfunction2
catch
  err2 = lasterror;
end

assert(~isequal(err1.message, err2.message))
assert(isequal(err1.identifier, err2.identifier))
assert(numel(err1.stack)==3)
assert(numel(err2.stack)==3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERSION 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subfunction1
subsubfunction1

function subsubfunction1
error('SOME:ERROR', 'something went wrong 1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERSION 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subfunction2
subsubfunction2

function subsubfunction2
ft_error('SOME:ERROR', 'something went wrong 2');
