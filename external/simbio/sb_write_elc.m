function sb_write_elc(pnt,labels,file_name,flag)
%
% sb_write_elc writes a set of electrodes on disk
%
% Use as
%   sb_write_elc(pnt,labels,file_name)
% 
% where:
%   pnt       are the positions of the electrodes in cartesian coordinates and millimiters
%   labels    is a cellstr of labels containing the names of the electrodes
%   flag      is 1 for DEEP electrodes

% Copyright (C) 2011, Felix Lucka

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information


% pnt units have to be in millimiters, cartesian coordinates
if nargin==3
  flag = 0;
end
N_sens = size(pnt,1);

% open the file and write the header
try
  fid = fopen(file_name, 'w');
  header{1} = 'ReferenceLabel	avg';
  header{2} = ['NumberPositions=    ' num2str(N_sens)];
  header{3} = 'UnitPosition	mm';
  header{4} = 'Positions';
  fprintf(fid,'%s\n%s\n%s\n%s\n',header{1},header{2},header{3},header{4});
  
  % write the points
  fprintf(fid, '%4.8e\t%4.8e\t%4.8e\n',pnt');
  
  % write additional fields
  % FIXME: is this optional?
  %  why the number of polygons is 0?
    mid{1} = 'NumberPolygons=    0';
    mid{2} = 'TypePolygons=   3';
    mid{3} = 'Polygons';
    fprintf(fid,'%s\n%s\n%s\n',mid{1},mid{2},mid{3});
  
  % write the labels
  fprintf(fid,'%s\n','Labels');
  for i=1:N_sens
    % true for depth electrodes
    if flag
      str = ['DEPTH ' labels{i}];
    else
      str = labels{i};
    end
    fprintf(fid,'%s\t',str);
  end
  
  % close the file
  fclose(fid);
catch
  disp('Error in writing the electrodes file')
  rethrow(ME)
end
