function y = subsref(x, index, inc)

% SUBSREF Return the value of a specified field in a config objects and increment its reference counter.

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



% function varargout = subsref(x, index, inc)
% % SUBSREF Return the value of a specified field in a config objects and increment its reference counter.
% 
% y = [];
% nargout
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
% nargout
