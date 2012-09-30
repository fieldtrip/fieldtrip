% UPDATE generates html help files for all classes and updates funct.html

d = '/Users/marcelge/Code/DMLT/html';
cd(d);

% read implemented methods

D = dir('../+dml/*.m');

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
  bf = exist(fname);
  fid = fopen(fname,'w');
  fprintf(fid,'%s',s);
  fclose(fid);
  if ~bf
    fprintf('added %s.html\n',nme);
    cd(d);
    system(sprintf('git add %s.html',nme));
    cd(ds);
  end
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
