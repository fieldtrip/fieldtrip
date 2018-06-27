% mff_createmff - create empty MFF file/folder
%
% Usage:
%   mff_createmff(mffFile);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file (MFF file/folder must
%            already exist)

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

function events = mff_createmff(mffURI)

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
import com.egi.services.mff.api.MFFFactory;
import com.egi.services.mff.api.MFFResourceType;
import com.egi.services.mff.api.LocalMFFFactoryDelegate;
import com.egi.services.mff.utility.ResourceUnmarshalException;
import com.egi.services.mff.api.Signal;
import com.egi.services.mff.api.SignalBlock;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
warning('on', 'MATLAB:Java:DuplicateClass');

% Test script for matlab binding of Java based MFF file access library.
% The MFFFactory class is used to access all resources contained in a
% MFF package file. This example loads an event track from the example file
% and writes out information about events in that file.

% create a factory
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate); %MFFResourceType.MFFResourceTypes.
mffFileResourceType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_MFFFile'));

if mfffactory.resourceExistsAtURI(mffURI, mffFileResourceType)
    disp('The resource already exists');
else
    disp('The resource does not exist');
    
    if mfffactory.createResourceAtURI(mffURI, mffFileResourceType)
        disp('Success creating MFF file.');
    else
        disp('An error occured creating the MFF file.');
    end
end
