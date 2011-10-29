function y=consolidateall(x)

% CONSOLIDATEALL applies consolidate to all the cells of a nested cell arrays
%
% Syntax: function y=consolidateall(x)
%
% Description:
% See description of CONSOLIDATE for details. By removing cell
% encapsulation of all the cells in the cell array, CONSOLIDATEALL will
% produce a dimensional structure. Therefore this fucntion will convert the
% product of XML2CELL into the output of XML2STRUCT
%
% See also: CONSOLIDATE, XML2STRUCT
%
% Jonas Almeida, almeidaj@musc.edu, 30 June 2003, MAT4NAT Tbox

z=consolidate(x);
if strcmp(class(z),'struct')
    f=fieldnames(z);
    for i=1:length(f)
        %disp(['(...).',f{i}])
        eval(['k=z.',f{i},';'])
        for j=1:length(k)
            eval(['y(',num2str(j),').',f{i},'=consolidateall(k{',num2str(j),'});'])
        end
    end
else
    y=z;
end