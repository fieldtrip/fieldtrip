function [pnt, dhk] = loadtri(name);

% LOADTRI reads vertices and triangles from a MBFYS triangulation file
%	[pnt, tri] = loadtri(filename)
%

% Thom Oostendorp, 2010

f = fopen(name);
if (f==-1)
  fprintf('\nCannot open %s\n\n', name);
  pnt=[];
  dhk=[];
  return;
end

[npnt,nr]=fscanf(f,'%d',1);
if (nr~=1)
  fclose(f);
  f=fopen(name);
  [magic ,nr]=fread(f,8,'char');
  if (char(magic')==';;mbftri') | (char(magic')==';;mbfpnt')
    fread(f,1,'char');
    hs=fread(f,1,'long');
    fread(f,1,'char');
    fread(f,1,'char');
    fread(f,1,'char');
    npnt=fread(f,1,'long');
    nc=fread(f,1,'long');
    fread(f,1,'char');
    ndhk=fread(f,1,'long');
    nv=fread(f,1,'long');
    if (npnt==0)
      pnt=[];
    else
      pnt=fread(f,[nc,npnt],'double');
      pnt=pnt';
    end;
    if (ndhk==0) | (char(magic')==';;mbfpnt')
      ndhk=0;
      dhk=[];
    else
      dhk=fread(f,[nv,ndhk],'long');
      dhk=dhk+1;
      dhk=dhk';
    end
  else
    fclose(f);
    fprintf('\n%s is of unknown binary type\n\n', name);
    pnt=[];
    dhk=[];
    return;
  end
else
  pnt  = fscanf(f, '%f', [4, npnt]);
  pnt  = pnt(2:4,:)';
  % if present, read the triangles
  if (~(feof(f)))
    [ndhk, count] = fscanf(f, '%d', 1);
    if (count ~= 0)
      dhk = fscanf(f, '%d', [4, ndhk]);
      dhk = dhk(2:4,:)';
    end
  else
    dhk = [];
  end
end
fclose(f);
return;

