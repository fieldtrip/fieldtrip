s = cell2mat(source.avg.mom(inside_idx));

sflip = s;

ss = s*s';
ssflip = ss;
% 
% for ii=1:length(s)
%     flipflag(ii) = -sign(mean(ssflip(ii,:)));
%     
%     sflip(ii,:) = flipflag(ii)*s(ii,:);
%     
% end

flipflag = diag(-sign(mean(ss)));

sflip = flipflag*sflip;

if(sum(sflip(:)) < 0)
    sflip = -sflip;
end

for ii=1:length(s)
    source.avg.mom{inside_idx(ii)} = sflip(ii,:);
end