% mff_importcategories - import information from MFF 'categories.xml' file
%
% Usage:
%   cat = mff_exportcategories(mffFile, version);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file
%  vesion  - file version (optional - default is 3)
%
% Output:
%  cat - Matlab structure containing informations contained in the MFF file.

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

function cat = mff_importcategories(mffFile, version)

cat = [];
p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

if nargin < 2
    version = 3;
end
if version == 0
    divider = 1000;
else
    divider = 1;
end

% Create a factory.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

%% create Segment to load time
categoriesRType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Categories'));
catURI = [ mffFile filesep 'categories.xml' ];
catsResource = mfffactorydelegate.openResourceAtURI(catURI, categoriesRType);
if ~isempty(catsResource)
    if catsResource.loadResource()
        categories = catsResource.getCategories();
        fprintf('Importing categories.xml ressource: %d categories\n', categories.size);
        
        for iCat = 1:categories.size
            category = categories.get(iCat-1);
            cat(iCat).name = char(category.getName());
            
            % Get the list of segments for this category.
            segments = category.getSegments();
            fprintf('Category %s, %d trials\n', char(category.getName()), segments.size);
            
            if ~isempty(segments)
                
                for iSeg = 1:segments.size
                    
                    segment = segments.get(iSeg-1);
                    cat(iCat).trials(iSeg).name = char(segment.getName());
                    cat(iCat).trials(iSeg).status = char(segment.getStatus());
                    cat(iCat).trials(iSeg).begintime = segment.getBeginTime()/divider;
                    cat(iCat).trials(iSeg).endtime = segment.getEndTime()/divider;
                    cat(iCat).trials(iSeg).eventbegin = segment.getEventBegin()/divider;
                    cat(iCat).trials(iSeg).eventend = segment.getEventEnd()/divider;
                    
                    if segment.getClockStartTimePresent()
                        cat(iCat).trials(iSeg).clockstarttime = char(segment.getClockStartTime());
                    end
                end
            end
        end
    end
end
