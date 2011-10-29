function sb_write_dip(pos,file_name)

% sb_write_dip writes a set of dipoles on disk
%
% Use as
%   sb_write_dip(pos,file_name)
% 
% where pos are the positions of the dipoles in cartesian coordinates and
% millimiters

% Copyright (C) 2011, Felix Lucka

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information


% number of dipoles (locations only)
N_sl = size(pos,1);

% write header
try
  fid = fopen(file_name,'w');
  fprintf(fid,['# ' num2str(3 * N_sl) ' dipole(s)\n']);
  fprintf(fid,'# Explanation with regard to the magnitude matrix:\n');
  fprintf(fid,'# Each column of the magnitude matrix below results in one forward simulated\n');
  fprintf(fid,'# sample. Each row corresponds to the respective dipole. If, like in this source file,\n');
  fprintf(fid,'# each column has only one nonzero 1 entry, then only the corresponding single dipole contributes\n');
  fprintf(fid,'# to the forward simulated field sample. Two nonzero 1 entries in a single column would result\n');
  fprintf(fid,'# in a field sample produced by the two corresponding simultaneously active dipoles etc..\n');
  fprintf(fid,'UnitPosition	mm\n');
  fprintf(fid,'UnitMoment	nAm\n');
  fprintf(fid,'UnitTime	ms\n');
  fprintf(fid,'NumberPositions=\t%d\n',3 * N_sl);
  fprintf(fid,'NumberTimeSteps=\t%d\n',1);
  fprintf(fid,'TimeSteps	0(1)%d\n',1);
  fprintf(fid,'FirstTimeStep	0\n');
  fprintf(fid,'LastTimeStep\t%d\n',1);
  fprintf(fid,'PositionsFixed\n');
  
  % write the points' locations (cartesian coordinates)
  for i=1:N_sl
    fprintf(fid,'%.6f\t%.6f\t%.6f\n',pos(i,1),pos(i,2),pos(i,3));
    fprintf(fid,'%.6f\t%.6f\t%.6f\n',pos(i,1),pos(i,2),pos(i,3));
    fprintf(fid,'%.6f\t%.6f\t%.6f\n',pos(i,1),pos(i,2),pos(i,3));   
  end
  
  % write the orientations components
  fprintf(fid,'MomentsFixed\n');
  for i=1:N_sl
    fprintf(fid,'%.6f\t%.6f\t%.6f\n',1,0,0);
    fprintf(fid,'%.6f\t%.6f\t%.6f\n',0,1,0);
    fprintf(fid,'%.6f\t%.6f\t%.6f\n',0,0,1);
  end
  
  % write the magnitudes
  fprintf(fid,'Magnitudes\n');
  for i=1:3*N_sl
    fprintf(fid,'%d\n',1);
  end
  
  % additional field
  fprintf(fid,'NoLabels');
  fclose(fid);
  
catch
  disp('Error in writing the dipoles file')
  rethrow(ME)
end
