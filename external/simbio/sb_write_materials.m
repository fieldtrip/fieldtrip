function sb_write_materials(filename,values,labels,resolution)
%
% SB_WRITE_MATERIALS writes a .mat file (read by the simbio-vgrid mesher)
% in which each tissue compartment is assigned to an integer value 
%
% Use as
%   sb_write_materials(values,labels,weights)
% 
% where:
%   values are the integer values assigned to the different compartments
%   labels is a cell with the labels descibing the comartments (default: {'scalp' 'skull' 'brain'})
%   weights are the weights assigned to the different compartments (default: [1.0 1.0 1.0])
%
% $Id$

% Copyright (C) 2011, Cristiano Micheli

if nargin<3 & length(values)==3
  weights = [1. 1. 1.];
  labels  = {'skin', 'skull', 'brain'}; % it assumes this order
elseif nargin<3 & length(values)~=3
  error('You should specify the labels for models with more (or less) than 3 compartments')
end

if isempty(resolution)
  resolution = 1;
end

% open the file and write the header
try
  fid = fopen(filename, 'w');
  fprintf(fid,'%s\n','material bg 1 0 0   2   1.0');
  for i=1:numel(values)
    fprintf(fid,'%s%s %d %d %d %d %0.1f\n','material ',labels{i},resolution,values(i),values(i),1,1);
  end
  fclose(fid);
catch
  disp('Error in writing the file')
  rethrow(ME)
end


