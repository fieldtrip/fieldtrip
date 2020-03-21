%% mff_valid.m
%  Matlab File
%  author Colin Davey
%  date 12/3/2013
%  Copyright 2013 EGI. All rights reserved.
%  Support routine for MFF Matlab code. Not intended to be called directly.
%
%  Tests whether a file is valid, and throws an exception if not. Gives an
%  informational warning if filename doesn't end in '.mff'. 
%%
function mff_valid(filePath)
valid = false;
filePathExist = exist(filePath);
% Check that filePath exists. 
if filePathExist == 0
    theException = MException('EGI_MFF:MFF_NO_EXIST', 'MFF does not exist.');
% Check that filePath isn't a regular file. 
elseif filePathExist == 2
    theException = MException('EGI_MFF:MFF_NOT_DIR', 'MFF is not a folder or package.');
% Check that filePath is a directory. 
elseif filePathExist ~= 7
    theException = MException('EGI_MFF:MFF_NOT_DIR_OR_FILE', 'MFF is not a folder, package or file.');
else
    % Check that filePath/info.xml exists and is a file. 
    if exist([filePath filesep 'info.xml']) ~= 2
        theException = MException('EGI_MFF:INVALID_NO_INFO', 'MFF is not valid. There is no info.xml file.');
    % Check that filePath/info1.xml exists and is a file. 
    elseif exist([filePath filesep 'info1.xml']) ~= 2
        theException = MException('EGI_MFF:INVALID_NO_INFO1', 'MFF is not valid. There is no info1.xml file.');
    % Check that filePath/signal1.bin exists and is a file. 
    elseif exist([filePath filesep 'signal1.bin']) ~= 2
        theException = MException('EGI_MFF:INVALID_NO_SIGNAL1', 'MFF is not valid. There is no signal1.bin file.');
    % Check the version. 
    else
        infoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info, 'info.xml', filePath);
        ver = infoObj.getMFFVersion;
        if ver ~= 3;
            theException = MException('EGI_MFF:WRONG_VER', 'MFF is not version 3. Please convert using EGI''s MFF File Converter. Contact supportteam@egi.com form more information.');
        else
            valid = true;
        end
    end
end
if ~valid
    throw(theException);
end
% Delete DS_Store file that gets generated when moving to a PC. 
fullDS_StorePath = [filePath filesep '._.DSStore'];
if exist(fullDS_StorePath, 'file') == 2
    delete(fullDS_StorePath);
end

% Give a warning if name doesn't end in '.mff'. 
mff_warning = '*** Filename does not end in ''.mff''. This can cause a problem for some versions of NetStation. ***';
if size(filePath,2) < size('.mff',2) + 1
%     warning('EGI_MFF:EXT_WARNING', mff_warning);
    fprintf('%s\n', mff_warning);
elseif ~(strcmp(lower(filePath(end-3:end)), '.mff'))
%     warning('EGI_MFF:EXT_WARNING', mff_warning);
    fprintf('%s\n', mff_warning);
end
