
function [freq] = specest_tfr_small(data, M, fsample, foi)
%convolves data in chan by time matrix (from one trial) with 1 wavelet for
%each specified frequency of interest
freq = zeros(size(data,1),length(foi), size(data,2));
  for k=1:size(data,1) %nchans
    for j=1:length(foi)
      cTmp = conv(data(k,:),M{j});
      cTmp = 2*(abs(cTmp).^2)/fsample;
      cTmp = cTmp(ceil(length(M{j})/2):length(cTmp)-floor(length(M{j})/2));
      freq(k,j,:) = cTmp;
    end
    
  end
  
  %output should be chan by freq by time