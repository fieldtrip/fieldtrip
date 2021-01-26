function [vertex, faces, texture, textureIdx] = read_obj_new(filename)

% Faster .obj reader for headshape information acquired using a structure
% sensor

% from http://boffinblogger.blogspot.com/2015/05/faster-obj-file-reading-in-matlab.html
% modified output structure and added scan for texture and textureIdx to handle vt
% Modifications by Robert Oostenveld and Robert Seymour

fid = fopen(filename);
if fid<0
  error(['Cannot open ' filename '.']);
end
[str, count] = fread(fid, [1,inf], 'uint8=>char');
fprintf('Reading %d characters from %s\n', count, filename);
fclose(fid);

tic;
vertex_lines = regexp(str,'v [^\n]*\n', 'match');
vertex = zeros(length(vertex_lines), 3);
for i = 1: length(vertex_lines)
  v = sscanf(vertex_lines{i}, 'v %f %f %f');
  vertex(i, :) = v';
end

texture_lines = regexp(str,'vt [^\n]*\n', 'match');
texture = zeros(length(texture_lines), 2);
for i = 1: length(texture_lines)
  vt = sscanf(texture_lines{i}, 'vt %f %f');
  texture(i, :) = vt';
end

face_lines = regexp(str,'f [^\n]*\n', 'match');
faces = zeros(length(face_lines), 3);
textureIdx = zeros(length(face_lines), 3);

for i = 1: length(face_lines)
  %     f = sscanf(face_lines{i}, 'f %d//%d %d//%d %d//%d');
  %     if (length(f) == 6) % face
  %         faces(i, 1) = f(1);
  %         faces(i, 2) = f(3);
  %         faces(i, 3) = f(5);
  %         continue
  %     end
  %     f = sscanf(face_lines{i}, 'f %d %d %d');
  %     if (length(f) == 3) % face
  %         faces(i, :) = f';
  %         continue
  %     end
  f = sscanf(face_lines{i}, 'f %d/%d %d/%d %d/%d');
  if (length(f) == 6) % face
    faces(i, 1) = f(1);
    faces(i, 2) = f(3);
    faces(i, 3) = f(5);
    textureIdx(i,1) = f(2);
    textureIdx(i,2) = f(4);
    textureIdx(i,3) = f(6);
    continue
  end
  f = sscanf(face_lines{i}, 'f %d/%d/%d %d/%d/%d %d/%d/%d');
  if (length(f) == 9) % face
    faces(i, 1) = f(1);
    faces(i, 2) = f(4);
    faces(i, 3) = f(7);
    textureIdx(i,1) = f(2);
    textureIdx(i,2) = f(5);
    textureIdx(i,3) = f(8);
    continue
  end
end
fprintf('.obj file took %f seconds to load\n',toc)
return