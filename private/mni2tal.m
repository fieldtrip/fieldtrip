function outpoints = mni2tal(inpoints)
% Converts coordinates from MNI brain to best guess
% for equivalent Talairach coordinates
% FORMAT outpoints = mni2tal(inpoints)
% Where inpoints is N by 3 or 3 by N matrix of coordinates
%  (N being the number of points)
% outpoints is the coordinate matrix with Talairach points
% Matthew Brett 10/8/99

dimdim = find(size(inpoints) == 3);
if isempty(dimdim)
  ft_error('input must be a N by 3 or 3 by N matrix')
end
if dimdim == 2
  inpoints = inpoints';
end

% Transformation matrices, different zooms above/below AC
upT = [
    0.9900    0.0000    0.0000    0.0000
    0.0000    0.9688    0.0460    0.0000
    0.0000   -0.0485    0.9189    0.0000
    0.0000    0.0000    0.0000    1.0000
      ];

downT = [
    0.9900    0.0000    0.0000    0.0000
    0.0000    0.9688    0.0420    0.0000
    0.0000   -0.0485    0.8390    0.0000
    0.0000    0.0000    0.0000    1.0000
    ];

tmp = inpoints(3,:)<0;  % 1 if below AC
inpoints = [inpoints; ones(1, size(inpoints, 2))];
inpoints(:, tmp) = downT * inpoints(:, tmp);
inpoints(:, ~tmp) = upT * inpoints(:, ~tmp);
outpoints = inpoints(1:3, :);
if dimdim == 2
  outpoints = outpoints';
end

