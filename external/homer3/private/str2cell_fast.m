function [C, k] = str2cell_fast(str, delimiters, options)

% Option tells weather to keep leading whitespaces. 
% (Trailing whitespaces are always removed)
if ~exist('options','var')
    options = '';
end

if ~strcmpi(options, 'keepblanks')
    str = strtrim(str);
end
str = deblank(str);

if ~exist('delimiters','var') || isempty(delimiters)
    delimiters{1} = sprintf('\n');
elseif ~iscell(delimiters)
    foo{1} = delimiters;
    delimiters = foo;
end

% Get indices of all the delimiters
k = 0;
for kk=1:length(delimiters)
    k = [k, find(str==delimiters{kk})];
end
if isempty(str)
    nstr = 0;
elseif ismember(str(end), delimiters)
    nstr = length(k)-1;
else
    nstr = length(k);
end
C = cell(1,nstr);
for ii = 1:nstr
    iS = k(ii)+1;
    if ii<length(k)
        iE = k(ii+1)-1;
    else
        iE = length(str);        
    end
    C{ii} = str(iS:iE);
end

% Verification code

% [C2, k2] = str2cell(str, delimiters, options);
% if length(C) ~= length(C2)
%     d3=1;
% end
% for ii = 1:length(C)
%     try 
%         if ~strcmp(C{ii}, C2{ii})
%             d1=1;
%         end
%     catch
%         d2=1;
%     end
% end

