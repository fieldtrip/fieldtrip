function out = full(fa)
% Convert to numeric form
% FORMAT full(fa)
% fa - a file_array
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


[vo{1:ndims(fa)}] = deal(':');
out = subsref(fa,struct('type','()','subs',{vo}));
