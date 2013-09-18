function [cmbindx, n, blocklabel] = blockindx2cmbindx(labelcmb, blockindx, block)

% this is a helper function which JM at some point wrote, but it is not
% clear to him anymore, what the function actually does. as soon as he
% finds out, the documentation will be updated.

y = cell(size(labelcmb));
z = cell(size(labelcmb));
okorig = zeros(size(labelcmb));
cmbindx = cell(size(block));
n       = cell(size(block));
for k = 1:numel(y)
  x1   = strfind(labelcmb{k},'[');
  x2   = strfind(labelcmb{k},']');
  y{k} = labelcmb{k}(1:x1-1);
  z{k} = labelcmb{k}(x1+1:x2-1);
  okorig(k) = numel(z{k});
end

for k = 1:numel(block)
  tmpsel = [];
  tmp    = {};
  for j = 1:numel(block{k})
    dummy  = find(ismember(blockindx{2},block{k}(j)));
    tmpsel = [tmpsel;dummy];
    n{k}(j,1) = numel(dummy);
  end
  tmp    = blockindx{1}(tmpsel);
  blocklabel{k} = cat(2,y{tmpsel});

  ok = okorig;
  for j = 1:numel(labelcmb)
    for jj = 1:numel(tmp)
      if isempty(strfind(z{j}, tmp{jj})),
        ok(j) = nan;
      else
        ok(j) = ok(j) - numel(tmp{jj});
      end 
    end
  end
  
  list = []; 
  if sum(ok(:,1)==0)==numel(tmp).^2,
    list = zeros(numel(tmp)*[1 1]);
    ok   = find(ok(:,1)==0);
    for j = 1:numel(ok)
      x1 = match_str(tmp, y(ok(j),1));
      x2 = match_str(tmp, y(ok(j),2));
      list(x1,x2) = ok(j);
    end
  end
  cmbindx{k} = list(:);
end
