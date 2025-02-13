function rgb = colorspec2rgb(c, n)

% COLORSPEC2RGB converts the string color specification into the corresponding RGB
% triplet(s), unless the string equals 'none', or is a hexadecimal (starting with #).
% The optional second input argument determines the number of rows in the output Nx3 matrix, 
% when applicable.
%
% If the first input argument equals 'none', or starts with a '#', the output will
% be the same as the input argument, and the assumption is that the downstream function
% that uses the colorspec knows how to deal with this. Otherwise, a Nx3 matrix or
% 1x3 vector will be returned.
%
% If the input string contains only characters from the
% following sequence 'ymcrgbwk', an Mx3 matrix will be returned, where M is the number
% of characters in the input string. If a second input argument N is defined (N>M), the 
% output will be expanded to have N number of rows.
% Otherwise, colorspec2rgb checks whether the input is a valid name for one of the colors
% htmlcolors or standardcolors. If this also fails, colorspec2rgb checks whether the colorspec
% defines a supported FieldTrip colormap. If this also fails, an error is thrown.
%
% If this also fails, an error is thrown
%
% see https://nl.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html
% and https://nl.mathworks.com/matlabcentral/fileexchange/48155-convert-between-rgb-and-color-names
%
% see also: HTMLCOLORS, STANDARDCOLORS, FT_COLORMAP

if nargin==1
  n = 1;
end

if ischar(c) && (isequal(c, 'none') || startsWith(c, '#'))
  rgb = c;
elseif isnumeric(c) && size(c,2)~=3
  ft_error('a numeric input should have 3 columns');
elseif all(ismember(c, 'ymcrgbwk'))
  % translate the single character color specifications
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
elseif ismember(c, htmlcolors)
  rgb = htmlcolors(c);
elseif ismember(c, standardcolors)
  rgb = standardcolors(c);
elseif ~isnumeric(c)
  try
    rgb = ft_colormap(c, n);
  catch
    ft_error('unkown color specification');
  end
end

if isnumeric(rgb) && n>size(rgb,1)
  rgb = repmat(rgb, [ceil(n/size(rgb,1)) 1]);
  rgb = rgb(1:n,:);
end

