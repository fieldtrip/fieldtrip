function [no, fc] = gridsurf(x, y, z, varargin)
%
% [no,fc]=gridsurf(x,y,z)
%    or
% [no,fc]=gridsurf(x,y,z, opt)
% [no,fc]=gridsurf(x,y,z, 'param1', value1, 'param2', value2, ...)
%
% convert a grid-shaped surface (used as input for surf) to a quad or triangular mesh
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%   x,y,z: parameter used as the input for surf()/mesh()
%   : a surface mesh triangle list (ne x 3)
%   opt: a list of optional parameters, currently surfacenorm supports:
%        'Type': [3|4|int] if set to 3 (default), output triangular mesh; if set to 4,
%                output a quad mesh where fc is an Nx4 array; otherwise,
%                output a quad mesh where fc is a cell array
%        'Nodup': [0|1] if set to 0 (default), no duplicated nodes in x/y/z are removed;
%                if set to 1, duplicated nodes are removed
%
% output:
%   no: output surface node coordinates
%   fc: output surface connections - if 'type' set to 3 or 4, fc is a numerical array,
%       otherwise, fc is a cell array with each element containing 4 integers
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

opt = varargin2struct(varargin{:});

s = size(x);
if (~all(s == size(y) & s == size(z)))
    error('x/y/z must be a 2D array of the same size');
end

no = [x(:), y(:), z(:)];
nodelen = size(no, 1);
fc = zeros((s(1) - 1) * s(2), 4);
row = [(1:s(1) - 1)', (2:s(1))', (s(1) + 2:2 * s(1))', (s(1) + 1:2 * s(1) - 1)'];

for i = 0:s(2) - 1
    fc(i * (s(1) - 1) + 1:(i + 1) * (s(1) - 1), :) = row + i * s(1);
end

fc(fc > nodelen) = fc(fc > nodelen) - nodelen;

if (jsonopt('nodup', 0, opt))
    [no, fc] = removedupnodes(no, fc);
end

outputtype = jsonopt('type', 3, opt);
if (outputtype == 3)
    fc = [fc(:, [1 2 3]); fc(:, [1 3 4])];
elseif (outputtype == 4)
    fc = num2cell(fc, 2);
end
