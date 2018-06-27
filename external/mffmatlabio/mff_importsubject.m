% mff_importsubject - import subject information from MFF ressources
%
% Usage:
%   subject = mff_importsubject(mffFile);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file
%
% Output:
%  subject - Matlab structure containing informations contained in the MFF file

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
% GNU General Public License for more detailsmff_setobj.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function subject = mff_importsubject(mffFile)

layout = [];
rVal = true;

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
import com.egi.services.mff.api.MFFFactory;
import com.egi.services.mff.api.MFFResourceType;
import com.egi.services.mff.api.LocalMFFFactoryDelegate;
import com.egi.services.mff.utility.ResourceUnmarshalException;
import com.egi.services.mff.api.SensorLayout;
import com.egi.services.mff.api.Sensor;
import com.egi.services.mff.api.Key;
import com.egi.services.mff.api.Neighbor;
warning('on', 'MATLAB:Java:DuplicateClass');

% Create an MFFFactory object.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
sURI = fullfile(mffFile, 'subject.xml');
objectType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Subject'));

sObject = mfffactory.openResourceAtURI(sURI, objectType);

variables = { 'Fields' 'array' { 'Name' 'char' {}; 'Data' 'string' {}; 'DataType' 'string' {} } };

subject = [];
if ~isempty(sObject)
    try
        if sObject.loadResource()
            subject = mff_getobj(sObject, variables);
        end
    catch
        disp('Failed to load subject ressource');
    end
end
