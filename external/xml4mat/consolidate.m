function y=consolidate(x)

%CONSOLIDATE field names in nested cell arrays produced by xml2mat of mbmling results
%
%Syntax: y=consolidate(x)
%
%Description:
% the nested cell arrays produced by XML2CELL emcapsulate the individual
% data structures. CONSOLIDATE will remove the cell encapsulation,
% returning the nested structure. CONSOLIDATEALL will apply CONSOLIDATE to
% all the cells in the array.
%
% See also CONSOLIDATEALL
%
% Jonas Almeida, almeidaj@musc.edu, 20 May 2003, MAT4NAT Tbox

if strcmp(class(x),'cell')
    if strcmp(class(x{1}),'struct')
        n=length(x);
        for i=1:n
            f(i)=fieldnames(x{i});
            I=strmatch(f(i),f(1:i-1));
            if ~isempty(I);
                j=I(1);
                %eval(['y.',f{j},'=consolidate(x{i}.',f{j},');'])
                if isfield(x{i},f{j})
                    eval(['y.',f{j},'{end+1}=consolidate(x{i}.',f{j},');'])
                    %warning(['field does not exist: ',f{j}])
                end
            else
                j=i;
                %eval(['y.',f{j},'=consolidate(x{i}.',f{j},');'])
                eval(['y.',f{j},'{1}=consolidate(x{i}.',f{j},');'])
            end
        end
    else
        y=x;
    end
else
    y=x;
end