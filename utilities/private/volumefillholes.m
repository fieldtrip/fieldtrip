function [output] = volumefillholes(input, along)

% VOLUMEFILLHOLES is a helper function for segmentations
%
% See also VOLUMETHRESHOLD, VOLUMESMOOTH

% check for any version of SPM
if ~ft_hastoolbox('spm')
  % add SPM8 to the path
  ft_hastoolbox('spm8', 1);
end

if nargin<2
  inflate = false(size(input)+2);                   % grow the edges along each dimension
  inflate(2:end-1, 2:end-1, 2:end-1) = (input~=0);  % insert the original volume
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
        slice=squeeze(input(i,:,:));
        im = imfill(slice,8,'holes');
        output(i,:,:) = im;
      end
      
    case 2
      for i=1:dim(2)
        slice=squeeze(input(:,i,:));
        im = imfill(slice,8,'holes');
        output(:,i,:) = im;
      end
      
    case 3
      for i=1:dim(3)
        slice=squeeze(input(:,:,i));
        im = imfill(slice,8,'holes');
        output(:,:,3) = im;
      end
      
    otherwise
      error('invalid dimension along which to slice the volume');
  end % switch
end % if nargin
