function [err, indx, SubLabel] = WhichSubBricks (InfoStat, Type, Label)
%
%   [err, indx, SubLabel] = WhichSubBricks (InfoStat, Type, Label)
%
%Purpose:
%
%
%
%Input Parameters:
%
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%   indx: a vector containing sub-brick indices (+1 from AFNI indexing)
%         with labels containing both of Type and Label strings
%   SubLabel: a vector of structures containing the label strings.
%         For example if indx = [1 4] then the labels you were searching
%         for are SubLabel(1).str and SubLabel(4).str which correspond to
%         AFNI sub-bricks 0 and 3 in AFNI's indexing world.
%
%
%Key Terms:
%
%More Info :
%
%
%
%
%     Author : Ziad Saad
%     Date : Fri Jul 18 16:56:58 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'WhichSubBricks';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

indx = zeros(1,InfoStat.DATASET_RANK(2));

%change the format of BRICK_LABS
cend = 0;
for (i=1:1:InfoStat.DATASET_RANK(2)),
   [err, SubLabel(i).str, cend] = NextString(InfoStat.BRICK_LABS,'~',cend+1);
end

cnt = 0;
for (i=1:1:InfoStat.DATASET_RANK(2)),
   if (~isempty(findstr(SubLabel(i).str, Label)) & ~isempty(findstr(SubLabel(i).str, Type))),
      cnt = cnt + 1;
      indx(cnt) = i;
   end
end

if (cnt == 0), indx = [];
else indx = indx(1:cnt);
end


err = 0;
return;

