% normalizeObjectArray() Normalize JSON array-of-objects to cell array of structs.
%
%   Input:
%       obj  - result of jsondecode for some JSON array
%       path - string used for error messages, e.g. 'EEGModality.Electrodes'
%
%   Output:
%       objs   - cell array of structs (each element one object)
%       errors - cell array of error strings (empty if ok)
%
% Author: Ramon Martinez-Cancino, Brain Products GmbH, 2025
%
% Copyright (C) 2025 Brain Products GmbH
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [objs, errors] = normalizeObjectArray(obj, path)
errors = {};
objs   = {};

if isstruct(obj)
    % jsondecode returned a struct array: one struct per element
    if isempty(obj)
        errors{end+1} = sprintf('"%s" must be a non-empty array of objects.', path);
        return;
    end
    objs = num2cell(obj); % 1xN cell, each a scalar struct
    return;
end

if iscell(obj)
    if isempty(obj)
        errors{end+1} = sprintf('"%s" must be a non-empty array of objects.', path);
        return;
    end
    % Check all elements are structs
    for i = 1:numel(obj)
        if ~isstruct(obj{i})
            errors{end+1} = sprintf('"%s{%d}" must be an object.', path, i);
        end
    end
    objs = obj;
    return;
end

% Anything else is wrong (numeric, string, etc.)
errors{end+1} = sprintf('"%s" must be an array of objects (struct array or cell array).', path);
end
