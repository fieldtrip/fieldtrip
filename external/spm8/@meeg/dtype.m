function res = dtype(obj)
% returns datatype of embedded file_array object
% FORMAT dtype(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: dtype.m 1533 2008-05-01 14:29:03Z spm $

res = obj.data.y.dtype;
