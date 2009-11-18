function m = nanmean(data,dim)
%NANMEAN take mean ignoring nans along any dimension

    if nargin<2, dim=1; end

    tmpdata = data;
    tmpdata(isnan(tmpdata(:))) = 0;
            
    m = mean(tmpdata,dim);

end