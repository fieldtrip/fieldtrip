function en = end(a,k,n)
% Overloaded end function for file_array objects.
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


dim = size(a);
if k>length(dim)
    en = 1;
else
    if n<length(dim)
    dim = [dim(1:(n-1)) prod(dim(n:end))];
    end
    en = dim(k);
end
