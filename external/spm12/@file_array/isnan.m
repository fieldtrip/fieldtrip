function out = isnan(fa)
% Logical array containing true where the elements of file_array are NaN's
% FORMAT isnan(fa)
% fa - a file_array
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


bs  = 10240;
m   = size(fa);
fa  = reshape(fa,prod(m),1);
n   = prod(m);
out = false(m);
for i=1:ceil(n/bs)
    ii      = ((((i-1)*bs)+1):min((i*bs),n))';
    tmp     = subsref(fa,struct('type','()','subs',{{ii}}));
    out(ii) = isnan(tmp);
end
