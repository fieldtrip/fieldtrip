function a = reshape(b,varargin)
% Overloaded reshape function for file_array objects
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


if length(struct(b))~=1, error('Can only reshape simple file_array objects.'); end

args = [];
for i=1:length(varargin)
    args = [args varargin{i}(:)'];
end
if prod(args)~=prod(b.dim)
    error('To RESHAPE the number of elements must not change.');
end
a = b;
a.dim = args;
