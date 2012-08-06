function b = iscompatwrapper(funcname)
%ISCOMPATWRAPPER Checks whether the specified function name will invoke a
% compatibility wrapper or not.

ftPath = which('ft_defaults');
ftPath = ftPath(1:end-numel('ft_defaults.m'));

x = which(funcname);
b = strcmp(x(1:end-numel(funcname)-2), [ftPath 'compat/']);

end

