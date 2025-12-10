function ft_test_update_dependency(depmat, inlist, outlist)

% FT_TEST_UPDATE_DEPENDENCY documentation is included inside ft_test
% documentation.
% 
% See also FT_TEST, READLINES, WRITELINES

% Copyright (C) 2023, Konstantinos Tsilimparis
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% Replace inlist with the full filename, including path
for i = 1:length(inlist) 
  fullname = which(inlist{i});

  if ~isempty(fullname)
    inlist{i} = fullname;
  else
    inlist{i} = fullfile(p,[f x]);
  end

end

% Replace outlist with the filename, excluding path
for i = 1:length(outlist) 
  [p, f, x]  = fileparts(outlist{i});
  outlist{i} = f;
end

% Update test scripts based on their dependencies
for  i = 1:length(inlist)

    % Read the lines from the file
    lines  = readlines(inlist{i});  

    % Loop through each line to check if any line contains "% DEPENDENCY"
    containsDependency = false;
    for k = 1:numel(lines)
        if contains(lines{k}, '% DEPENDENCY')
            containsDependency = true;
            break; 
        end
    end
    
    if containsDependency
        line = find(startsWith(lines, '% DEPENDENCY'));
    
        if length(line) == 1
            lines{line}(length('% DEPENDENCY')+1:end) = []; % Delete the old content after '% DEPENDENCY '
        end
    
        for j = 1:length(outlist)
            if depmat(i,j) == 2
                lines{line} = [lines{line}, ' ', outlist{j}];  % Append the new content from outlist
            end
        end
    
        % Save the modified content back to the file
        writelines(lines, inlist{i}); % Note: writelines exists since MATLAB 2022a

    else 
        
        line  = find(startsWith(lines, '% WALLTIME'));        
        lines = [lines(1:line);  '% DEPENDENCY'; lines(line+1:end)]; % List the dependencies a line after the line that contains the "WALLTIME"

         if length(line) == 1
             for j = 1:length(outlist)
                  if depmat(i,j) == 2
                      lines{line+1} = [lines{line+1}, ' ', outlist{j}];  % Append the new content from outlist
                  end
              end
         end

        % Save the modified content back to the file
        writelines(lines, inlist{i}); % Note: writelines exists since MATLAB 2022a

    end
end 
