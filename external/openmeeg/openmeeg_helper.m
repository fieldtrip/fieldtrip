function lf=openmeeg_helper(pos,vol,dimdip)
%
% Splits the calculation of leadfields in smaller chunks
% pos: positions of dipoles
% vol: volume conductor
% dimdip: number of dipoles to calculate at a time

lf   = [];
ndip = size(pos,1);
splitdim = dimdip;
numsol   = fix(ndip./splitdim);
nleft    = rem(ndip,splitdim);
lf = zeros(size(vol.mat,1),ndip*3);
if isfield(vol,'mat')
  if ft_hastoolbox('openmeeg', 1);
    for i = 1:numsol
      dsm = openmeeg_dsm(pos((splitdim*(i - 1) + 1):splitdim*i,:),vol);
      tmp = vol.mat*dsm;
      lf(:,(splitdim*(i - 1)*3 + 1):splitdim*i*3) = tmp;
    end
    dsm = openmeeg_dsm(pos(ndip-nleft+1:ndip,:),vol);
    tmp = vol.mat*dsm;
    lf(:,ndip-nleft+1:ndip,:) = tmp;
  else
    error('Openmeeg toolbox not installed')
  end
else
  error('No system matrix is present, BEM head model not calculated yet')
end
