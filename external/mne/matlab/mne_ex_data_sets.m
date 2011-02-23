function [res] = mne_ex_data_sets(fname)
%   
%   Find all evoked response data from a given file
%
%   [res] = mne_ex_data_sets(fname)
%
%   fname   - Name of the file to look at
%   res     - Structure containing the result
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/21 17:31:07  msh
%   Added the examples.
%   Modified the formatting of some informative output.
%
%
me='MNE:mne_ex_data_sets';

if nargin ~= 1
    error(me,'Wrong number of arguments');
end

try
    res = fiff_find_evoked(fname);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end

if isempty(res)
    fprintf(1,'No evoked response data sets in %s\n',fname);
else
    fprintf(1,'Found %d evoked response data sets in %s :\n',length(res),fname);
    for k = 1:length(res)
        fprintf(1,'\t%s (%s)\n',res(k).comment,res(k).aspect_name);
    end
end


return;

end

