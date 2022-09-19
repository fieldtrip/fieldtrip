function outpoints = tal2mni(inpoints)
% Converts coordinates to MNI brain best guess
% from Talairach coordinates
% FORMAT outpoints = tal2mni(inpoints)
% Where inpoints is N by 3 or 3 by N matrix of coordinates
%  (N being the number of points)
% outpoints is the coordinate matrix with MNI points
% Matthew Brett 2/2/01

dimdim = find(size(inpoints) == 3);
if isempty(dimdim)
  ft_error('input must be a N by 3 or 3 by N matrix')
end
if dimdim == 2
  inpoints = inpoints';
end

% Transformation matrices, different zooms above/below AC
rotn  = [
    1.0000    0.0000    0.0000    0.0000
    0.0000    0.9988    0.0500    0.0000
    0.0000   -0.0500    0.9988    0.0000
    0.0000    0.0000    0.0000    1.0000
    ];

upz = [
    0.9900    0.0000    0.0000    0.0000
    0.0000    0.9700    0.0000    0.0000
    0.0000    0.0000    0.9200    0.0000
    0.0000    0.0000    0.0000    1.0000
     ];

downz = [
    0.9900    0.0000    0.0000    0.0000
    0.0000    0.9700    0.0000    0.0000
    0.0000    0.0000    0.8400    0.0000
    0.0000    0.0000    0.0000    1.0000
    ];

inpoints = [inpoints; ones(1, size(inpoints, 2))];
% Apply inverse translation
inpoints = inv(rotn)*inpoints;

tmp = inpoints(3,:)<0;  % 1 if below AC
inpoints(:, tmp) = inv(downz) * inpoints(:, tmp);
inpoints(:, ~tmp) = inv(upz) * inpoints(:, ~tmp);
outpoints = inpoints(1:3, :);
if dimdim == 2
  outpoints = outpoints';
end

