function dip = sb_read_dip(filename)  
% reads dipole file for simbio toolbox
fid = fopen(filename,'r');
  s  = fgetl(fid); 
  [dum] = sscanf(s,'%c%d'); 
  Ndip  = dum(2);
  tline = fgetl(fid);
  while isempty(findstr(tline,'PositionsFixed'))
    tline=fgetl(fid);
  end
  dip = zeros(Ndip,3);
  for i=1:Ndip
    tline=fgetl(fid);
    dip(i,:)= sscanf(tline,'%f%f%f')';
  end
  fclose(fid);