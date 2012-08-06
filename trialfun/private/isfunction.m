function b = isfunction(funcname)
%ISFUNCTION tests whether the function of the specified name is a callable
% function on the current Matlab path.
%
% Note that this is *not* equivalent to calling exist(funcname, 'file'),
% since that will return 7 in case funcname exists as a folder.

b = ~isempty(which(funcname));
end