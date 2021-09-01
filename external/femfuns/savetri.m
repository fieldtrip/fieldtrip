function savetri(name,pnt,dhk);

% SAVETRI saves vertices and triangles to a MBFYS triangulation file
%	savetri(filename, pnt, tri)
%

[nr nc]=size(pnt);
if (nc~=3)
  if (nr==3)
    pnt=pnt';
  else
    disp('Savetri error: Vertex parameter should have 3 rows or columns');
    return;
  end
end
[nr nc]=size(dhk);
if (nc~=3)
  if (nr==3)
    dhk=dhk';
  else
    disp('Savetri error: Triangle parameter should have 3 rows or columns');
    return;
  end
end

[npnt nc]=size(pnt);
[ndhk nv]=size(dhk);
f = fopen(name, 'w');

% see triio.cc for description of header

fwrite(f,';;mbftri','char');
fwrite(f,0,'char');
fwrite(f,20,'long');
fwrite(f,1,'char');
fwrite(f,1,'char');
fwrite(f,8,'char');
fwrite(f,npnt,'long');
fwrite(f,nc,'long');
fwrite(f,4,'char');
fwrite(f,ndhk,'long');
fwrite(f,nv,'long');
fwrite(f,pnt','double');
fwrite(f,dhk'-1,'long');
fclose(f);


% This was the ascii  version
%fprintf(f,'%d\n ',nver);
%for i=1:nver;
%      fprintf(f,'%d %8.4f %8.4f %8.4f\n',i,VER(i,1:3));
%end
%
%[ntri dim]=size(ITRI);
%fprintf(f,'%d\n ',ntri);
%for i=1:ntri;
%      fprintf(f,'%d %d %d %d\n',i,ITRI(i,1:3));
%end
%fprintf('\ntriangle specs written to file: %s\n',fn); 
%fclose(f);
