function rgb = char2rgb(c)

% CHAR2RGB converts the line color character into the corresponding RGB triplet
%
% see https://nl.mathworks.com/help/matlab/ref/colorspec.html
% and https://nl.mathworks.com/matlabcentral/fileexchange/48155-convert-between-rgb-and-color-names

rgb = zeros(numel(c), 3);
for i=1:numel(c)
  switch c(i)
    case 'y'
      rgb(i,:) = [1 1 0];
    case 'm'
      rgb(i,:) = [1 0 1];
    case 'c'
      rgb(i,:) = [0 1 1];
    case 'r'
      rgb(i,:) = [1 0 0];
    case 'g'
      rgb(i,:) = [0 1 0];
    case 'b'
      rgb(i,:) = [0 0 1];
    case 'w'
      rgb(i,:) = [1 1 1];
    case 'k'
      rgb(i,:) = [0 0 0];
    otherwise
      ft_error('unknown color specification');
  end
end