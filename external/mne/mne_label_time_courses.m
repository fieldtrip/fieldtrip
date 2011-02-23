function [ values, times, vertices ] = mne_label_time_courses(labelfile,stcfile)
%
% function [ values, times ] = mne_label_time_courses(labelfile,stcfile)
%
% Extract the time courses corresponding to a label file from an stc file
%
% labelfile - The name of the label file
% stcfile   - The name of the stc file (must be on the same subject and
% hemisphere as the stc file
%
% values    - The time courses
% times     - The time points
% vertices  - The vertices corresponding to the time points
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%
me='MNE:mne_label_time_courses';

if nargin ~= 2
    error(me,'Incorrect number of arguments');
end
try
    stc = mne_read_stc_file(stcfile);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end

try
    lab = mne_read_label_file(labelfile);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end

[vertices,ia,ib] = intersect(double(stc.vertices),double(lab.vertices));
if length(vertices) == 0 
    error(me,'No vertices match the label in the stc file');
end

values = stc.data(ia,:);
times  = zeros(1,size(stc.data,2));
for k = 0:length(times)-1
    times(k+1) = stc.tmin + k*stc.tstep;
end

end

