% UPDATE generates html help files for all classes and updates funct.html

d = pwd;

% read implemented methods

D = dir('../+dml/*.m');
if isempty(D)
  error('update should be called from within the toolbox html directory');
end

% write method help

s = which('help2html');

f = strfind(s,'/');
if isempty(f)
  cd(d);
  error('path not found');
end

ds = s(1:(f(end)-1));
cd(ds);

for i=1:length(D)
  nme = D(i).name(1:(end-2));
  s = help2html(['dml.' nme]);  
  fname = fullfile(d,[nme '.html']);
  if ~exist(fname)
    fprintf('added %s.html\n',nme);
  end
  fid = fopen(fname,'w');
  fprintf(fid,'%s',s);
  fclose(fid);
end

% write funct.m

fid = fopen(fullfile(d,'funct.m'),'w');
fprintf(fid,'%%%% DML Toolbox Functions\n%%\n');
for i=1:length(D)
  nme = D(i).name(1:(end-2));
  x = help(sprintf('dml.%s',nme));
  b = strfind(lower(x),nme)+length(nme)+1; b = b(1);
  e = min(b+50,(strfind(x,'.')-1)); e = e(1);
  x = x(b:e);
  fprintf(fid,['%% * <' nme '.html |dml.' nme '|> - ' x '\n']);
end
fclose(fid);

cd(d);

publish('funct.m','outputDir',d);
