d = dir('test_*.m');
filename = {d.name}';

medianmem = nan(size(filename));
walltime = nan(size(filename));

for i=1:numel(filename)
  disp(i/numel(filename));
  
  fid = fopen(filename{i});
  
  while ~feof(fid)
    l = fgetl(fid);
    
    if contains(l, 'WALLTIME')
      val = split(l, ' ');
      val = split(val{end}, ':');
      h = str2double(val{1});
      m = str2double(val{2});
      s = str2double(val{3});
      walltime(i) = (h*3600 + m*60 + s)/3600; % in fractional hours
    end
    
    if contains(l, 'MEM')
      val = split(l, ' ');
      val = val{3};
      r = str2double(val(1:end-2));
      if contains(lower(val), 'mb')
        r = r/1024;
      end
      mem(i) = r; % in GB
    end
    
  end % while
  fclose(fid);
  
end % for
