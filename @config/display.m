function display(x)

% DISPLAY Display function for a config object.

% Copyright (C) 2012-2015, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% this flag detemines how much feedback is given
% which is usefull for debugging
fb = false;

if ~isempty(inputname(1))
  fprintf('\n');
  fprintf('%s =\n\n', inputname(1));
end

if fb
  fprintf('----------- value      -----------\n');
end

if numel(x)==0
  disp('1x1 config array with no fields.');
elseif numel(x)==1
  disp(x.value);
else
  siz = sprintf('%dx', size(x)); % construct a string like 1x2x3x, the last 'x' has to be removed
  siz = siz(1:(end-1));
  fprintf('%s config array with fields:\n', siz);
  key = fieldnames(x);
  for i=1:length(key)
    fprintf('    %s\n', key{i});
  end
end

if fb
  fprintf('----------- assignment -----------\n');
  disp(x.assign);
  fprintf('----------- reference  -----------\n');
  disp(x.reference);
end

