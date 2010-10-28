function STM = ft_omri_slice_time_init(X0, TR, deltaT);
% function STM = ft_omri_slice_time_init(X0, TR, deltaT);
%
% Initialize simple slice time correction structure.
% The algorithm will use plain linear interpolation.

% 2010 S.Klanke

STM.dims = size(X0)
if nargin == 3
  if STM.dims(3)~=length(deltaT) | any(deltaT>TR) | any(deltaT<0) 
    error 'You must pass in a correct timing description vector';
  end
else
  deltaT = (0:(STM.dims(3)-1))*TR/STM.dims(3);
end

STM.deltaT = deltaT(:)';
STM.TR = TR;

xy = STM.dims(1)*STM.dims(2);

STM.weightOld = reshape(ones(xy,1)*(deltaT/TR), STM.dims);
STM.weightNew = reshape(ones(xy,1)*(1-deltaT/TR), STM.dims);

STM.lastScan = X0;