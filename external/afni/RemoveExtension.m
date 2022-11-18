function [S, xtr] = RemoveExtension (S, xt)
%
%   [Sx, xtr] = RemoveExtension (S, [xt])
%
%Purpose:
%   removes the extension xt from the end of S
%
%
%Input Parameters:
%   S : string
%   xt: string of characters to be removed from the end of S
%      xt can be | delimited strings. Trainling blanks will be removed
%   	if xt is empty (default) then the characters following and including
%       the first . will be removed.
%
%Output Parameters:
%   Sx : string, S without the extension
%   xtr: The extension that was removed
%
%Key Terms:
%
%More Info :
%   S = 'ajh_d.BRIK';
%    [St, xtr] = RemoveExtension (S,'.HEAD|.BRIK')
%
%   S = 'ajh_d';
%   [St, xtr] = RemoveExtension (S,'.HEAD|.BRIK')
%
%     Author : Ziad Saad
%     Date : Mon Oct 9 15:08:08 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'RemoveExtension';

if nargin < 2,
	xt = '';
end

xtr='';
if isempty(xt)
    first_dot_pos=find(S=='.',1,'first');
    if ~isempty(first_dot_pos)
        [S,xtr]=cut_string(S,first_dot_pos);
    end
else
    %find the number of words in xt
    end_strings=get_end_strings(xt);
    n = numel(end_strings);

    for i=1:n
        end_string=end_strings{i};
        if string_ends_with(end_string,S)
            cut_pos=numel(S)-numel(end_string)+1;
            [S,xtr]=cut_string(S,cut_pos);
        end
    end

end

function [first_part,last_part]=cut_string(s, start_last_part)
    first_part=s(1:(start_last_part-1));
    last_part=s(start_last_part:end);


function tf=string_ends_with(end_str,s)
    tf=numel(s)>=numel(end_str) && ...
            strcmp(s(end+((1-numel(end_str)):0)),end_str);

function end_strings=get_end_strings(xt)
    m=regexp(xt,'\s*(\S[^|]*\S)\s*([|]|$)\s*','tokens');

    % get first element of each match
    end_strings=cellfun(@(x)x{1},m,'UniformOutput',false);

