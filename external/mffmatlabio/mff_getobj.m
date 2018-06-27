% mff_getobj - import information from MFF object
%
% Usage:
%   values = mff_setobj(mffObj, structure);
%
% Inputs:
%  mffObj    - MFF object of any type
%  structure - description of variables to assign
%
% Outputs:
%  values    - values for the structure
% 
% Note: this function may take any MFF object structure

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function matStruct = mff_getobj(mffObj, variables)

for iVar = 1:size(variables,1)
    varName  = variables{iVar,1};
    varType  = variables{iVar,2};
    varArray = variables{iVar,3};
    
    if ~isequal(varType, 'array')
        eval( [ 'matStruct.(lower(varName)) = ' varType '(mffObj.get' varName '());' ] );
    else
        if isempty(varName)
            for iArray = 1:mffObj.size
                matStruct.value(iArray) = mffObj.get(iArray-1);
            end
        else
            eval( [ 'tmparray = mffObj.get' varName '();' ] );
            if ~isempty(tmparray)
                for iArray = 1:tmparray.size
                    tmpval = tmparray.get(iArray-1);
                    if ~isempty(varArray)
                        matStruct.(lower(varName))(iArray) = mff_getobj(tmpval, varArray);
                    else
                        matStruct.(lower(varName))(iArray) = tmpval;
                    end
                end
            else
                matStruct.(lower(varName)) = [];
            end
        end
    end
end
