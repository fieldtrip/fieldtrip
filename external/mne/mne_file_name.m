function [res] = mne_file_name(dir,name)
%
%   [name] = mne_file_name(dir,name)
%
%   Compose a file name under MNE_ROOT
%
%   dir     - Name of the directory containing the file name
%   name    - Name of the file under that directory
%   

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%

me='MNE:mne_file_name';

if ~ispref('MNE','MNE_ROOT')
    mne_root = getenv('MNE_ROOT');
    if isempty(mne_root)
        error(me,'MNE_ROOT not defined');
    end
else
    mne_root=getpref('MNE','MNE_ROOT');
end

if nargin == 2
    res = sprintf('%s/%s/%s',mne_root,dir,name);
elseif nargin == 1
    res = sprintf('%s/%s',mne_root,dir);
else
    error(me,'incorrect number of arguments');
end 

return;

end
        
