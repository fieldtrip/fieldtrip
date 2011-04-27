function out = horzcat(varargin)
%HORZCAT Horizontal concatenation of data acquisition objects.


% Determine if the objects are valid and have the same parent.
valid_parent = 1;
for i = 1:length(varargin)
    valid_parent = valid_parent&isa(varargin{i},'chanhandle');
    if ~valid_parent,break;end;
end;
if valid_parent == 0
   error('Invalid object.');
end

%Concatenate the handles of each input into one object.
c=[];
for i = 1:nargin
   if ~isempty(varargin{i}),
      if isempty(c),
         c=varargin{i};
      else
         try
            c.handle = [c.handle varargin{i}.handle];
         catch
            error(lasterr);
         end
      end      
   end
end

% Determine if a matrix of channels was constructed if so error
% since only vectors are allowed.
if length(c.handle) ~= prod(size(c.handle))
   error('Only a row or column vector of channels can be created.')
end

% Assign the new channel vector to the output.  
out = c;
