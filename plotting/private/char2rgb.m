function rgb = char2rgb(c)

% CHAR2RGB converts the line color character or string into the corresponding RGB
% triplet(s), unless the string equals 'none'. In that case the function returns
% 'none', under assumption that the downstream function knows how to handle this.
% Subsequently, if the input string contains only characters from the
% following sequence 'ymcrgbwk', an Nx3 matrix will be returned. Otherwise,
% char2rgb checks whether the input is a valid name of one of the colors in
% htmlcolors. If this also fails, char2rgb checks whether the specified
% color is from the following list:
%
%   red               = [255   0   0]/255;
%   green             = [  0 192   0]/255;
%   blue              = [  0   0 255]/255;
%   magenta           = [255 255   0]/255;
%   cyan              = [  0 255 255]/255;
%   yellow            = [255 255   0]/255;
%   white             = [255 255 255]/255;
%   black             = [  0   0   0]/255;
%
%   skull             = [140  85  85]/255
%   cortex            = [255 213 119]/255;
%   cortex_light      = [199 194 169]/255;
%   cortex_dark       = [100  97  85]/255;
%   skin              = [249 223 192]/255;
%   skin_light        = [249 223 192]/255;
%   skin_medium_light = [225 194 158]/255;
%   skin_medium       = [188 142 106]/255;
%   skin_medium_dark  = [155 102  65]/255;
%   skin_dark         = [ 91  71  61]/255;
%
% If this also fails, an error is thrown
%
% see https://nl.mathworks.com/help/matlab/ref/colorspec.html
% and https://nl.mathworks.com/matlabcentral/fileexchange/48155-convert-between-rgb-and-color-names
%
% see also: HTMLCOLORS

list = {'red', 'green', 'blue', 'magenta', 'cyan', 'yellow', 'white', 'black', ...
  'skull', 'cortex', 'cortex_light', 'cortex_dark', 'skin', 'skin_light', ...
  'skin_medium_light', 'skin_medium', 'skin_medium_dark', 'skin_dark'};

if isequal(c, 'none')
  rgb = 'none';
elseif all(ismember(c, 'ymcrgbwk'))
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
elseif ismember(c, list)
  rgb = feval(str2func(c));
else
  ft_error('unkown color specification');
end
