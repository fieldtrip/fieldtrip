function [dmeg, dspc] = block_delay(bmeg, bspc, d)

Nblock = length(bmeg);
dmeg = cell(Nblock,1);
dspc = cell(Nblock,1);
for bi=1:Nblock
    common = min(size(bmeg{bi},1),length(bspc{bi}));
    dmeg{bi} = bmeg{bi}((1+d):common,:);
    dspc{bi} = bspc{bi}(1:(common-d));
end
dmeg = cell2mat(dmeg);
dspc = cell2mat(dspc);