function [output] = volumefillholes(input, along)

% VOLUMEFILLHOLES is a helper function for segmentations
%
% See also VOLUMETHRESHOLD, VOLUMESMOOTH, VOLUMEPAD, VOLUMESELECTLARGEST

% ensure that SPM is available, needed for spm_bwlabel
ft_hastoolbox('spm8up', 3) || ft_hastoolbox('spm2', 1);

if nargin<2

  inflate = volumepad(input, 1);                    % pad the volume
  [lab, num] = spm_bwlabel(double(~inflate), 18);   % note that 18 is consistent with imfill, 26 is not

  if num>1
    inflate(lab~=lab(1)) = true;
    output  = inflate(2:end-1, 2:end-1, 2:end-1);   % trim the edges
  else
    output = input;
  end
  
else
  output = input;
  dim    = size(input);
  switch along
    case 1
      for i=1:dim(1)
        slice = reshape(input(i,:,:),dim([2 3]));
        im = imfill(slice,8,'holes');
        output(i,:,:) = im;
      end
      
    case 2
      for i=1:dim(2)
        slice = reshape(input(:,i,:),dim([1 3]));
        im = imfill(slice,8,'holes');
        output(:,i,:) = im;
      end
      
    case 3
      for i=1:dim(3)
        slice = reshape(input(:,:,i),dim([1 2]));
        im = imfill(slice,8,'holes');
        output(:,:,i) = im;
      end
      
    otherwise
      ft_error('invalid dimension %d to slice the volume', along);
  end % switch along
end % if nargin
