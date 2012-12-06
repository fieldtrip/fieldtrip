function y = subsref(x, index, inc)

% SUBSREF Return the value of a specified field in a config objects and increment its reference counter.

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

if nargin<3
  inc = true;
end

if length(index)==1
  switch index.type
    case '.'
      y = get(x, index.subs, inc);
    case '{}'
      error('Cell contents reference from a non-cell array object.');
    case '()'
        y = x(index.subs{1});
    otherwise
      error('Incorrect contents reference');
  end
else
  % use recursion to find the subfield that is being indexed
  y = subsref(subsref(x, index(1)), index(2:end));
end


% TEST CODE
% function varargout = subsref(x, index, inc)
% % SUBSREF Return the value of a specified field in a config objects and increment its reference counter.
% 
% y = [];
% if nargin<3
%   inc = true;
% end
% 
% if length(index)==1
%   switch index.type
%     case '.'
%       if numel(x)>1
%         y = cell(1,numel(x));
%         % fields from multiple (sub)structures are requested, loop over each
%         for iobj = 1:numel(x)
%           y{iobj} = get(x(iobj), index.subs, inc);
%         end
%       else
%         y{1} = get(x, index.subs, inc);
%       end
%     case '{}'
%       error('Cell contents reference from a non-cell array object.');
%     case '()'
%         y{1} = x(index.subs{1});
%     otherwise
%       error('Incorrect contents reference');
%   end
% else
%   % use recursion to find the subfield that is being indexed
%     y{1} = subsref(subsref(x, index(1)), index(2:end));
% end
% varargout = y;
