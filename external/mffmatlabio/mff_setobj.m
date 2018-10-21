% mff_setobj - export information into MFF object
%
% Usage:
%   mffObj = mff_setobj(mffObj, structure, values);
%
% Inputs:
%  mffObj    - MFF object of any type
%  structure - description of variables to assign
%  values    - values for the structure
% 
% Output:
%  mffObj    - MFF object of any type
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

function mffObj = mff_setobj(mffObj, matStruct, variables)

for iVar = 1:size(variables,1)
    varName  = variables{iVar,1};
    varType  = variables{iVar,2};
    varArray = variables{iVar,3};
    
    if ~isequal(varType, 'array')
        eval( [ 'mffObj.set' varName '(matStruct.(lower(varName)));' ] );
    else
        if isempty(varName)
            % special case for TilingSet
            for iArray = 1:length(matStruct.value)
                mffObj.add(javaObject('java.lang.Integer', matStruct.value(iArray)));
            end
        else
            if ~isempty(matStruct.(lower(varName)))
                jList = javaObject('java.util.ArrayList');
                for iArray = 1:length(matStruct.(lower(varName)))
                    if isempty(varArray) || isempty(varArray{1})
                        mffObj2 = javaObject('java.util.ArrayList');
                    else
                        mffObj2 = javaObject( [ 'com.egi.services.mff.api.' varName(1:end-1) ]);
                    end
                    if isreal(matStruct.(lower(varName))(iArray))
                        jList.add(javaObject('java.lang.Integer', matStruct.(lower(varName))(iArray)));
                    else
                        jList.add(mff_setobj(mffObj2, matStruct.(lower(varName))(iArray), varArray));
                    end
                end
                eval( [ 'mffObj.set' varName '(jList);' ]);
            end
        end
    end
end
