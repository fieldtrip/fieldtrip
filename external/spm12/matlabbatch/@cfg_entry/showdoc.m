function str = showdoc(item, indent)

% function str = showdoc(item, indent)
% Display help text for a cfg_entry item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showdoc.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

str = showdoc(item.cfg_item, indent);
switch item.strtype
    case {'e'},
        str{end+1} = 'Evaluated statements are entered.';
        str{end+1} = shownum(item.num);
    case {'n'},
        str{end+1} = 'Natural numbers are entered.';
        str{end+1} = shownum(item.num);
    case {'r'},
        str{end+1} = 'Real numbers are entered.';
        str{end+1} = shownum(item.num);
    case {'w'},
        str{end+1} = 'Whole numbers are entered.';
        str{end+1} = shownum(item.num);
    case {'s'},
        str{end+1} = 'A String is entered.';
        if isempty(item.num)
            str{end+1} = 'The character array may have arbitrary size.';
        elseif isfinite(item.num(2))
            str{end+1} = sprintf(['The string must have between %d and %d ' ...
                                'characters.'], item.num(1), ...
                                 item.num(2));
        else
            str{end+1} = sprintf(['The string must have at least %d ' ...
                                'characters.'], item.num(1));
        end;
end;

function numstr = shownum(num)
if isempty(num)
    numstr = ['The entered data may have an arbitrary number of dimensions ' ...
              'and elements.'];
else
    for k=1:numel(num)
        if isfinite(num(k))
            numstr1{k} = sprintf('%d',num(k));
        else
            numstr1{k} = 'X';
        end;
    end;
    numstr = sprintf('%s-by-', numstr1{:});
    numstr = sprintf('An %s array must be entered.', numstr(1:end-4));
end;