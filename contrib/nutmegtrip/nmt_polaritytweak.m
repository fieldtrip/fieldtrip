function source = nmt_polaritytweak(cfg,source);

inside_idx = find(source.inside);

s = cell2mat(source.avg.mom(inside_idx));
ori = cell2mat(source.avg.ori(inside_idx));

tsel = dsearchn(source.time',cfg.toilim');
tsel = tsel(1):tsel(2);


if(1)
    [~,b]=max(nmt_rownorm(s(:,tsel)));
    flipper=sparse(diag(sign(s(b,tsel)*s(:,tsel)'))); 
    s = flipper*s;
    ori = ori*flipper; 
end

% if(0)
%     for ii=1:size(ss,1)
%         ss(ii,ii)=0;
%     end
%     [ssranked,idx] = sort(abs(ss(:)),'descend');
%
% end
% if(0)
% for ii=1:10
% ss = s(:,tsel)*s(:,tsel)';
% ssabs = abs(s(:,tsel))*abs(s(:,tsel)');
% flipper = diag(sign(max(ss./ssabs)));
% s = flipper*s;
% end
% 
% flipflag = diag(sign(mean(s(:,tsel),2)));
% s = flipflag*s;
% 
% switch(2)
%     case 1 % simple method: ensure sign of overall mean is positive:
%         flipflag = diag(sign(mean(s(:,tsel),2)));
%     case 2 % more complex: ensure covariances are more positive (so that voxels move together)
%         ss = s(:,tsel)*s(:,tsel)';
%         ssflip = ss;
%         
%         flipflag = diag(sign(mean(ss)));
%         
% end
% end


for ii=1:length(s)
    source.avg.mom{inside_idx(ii)} = s(ii,:);
    source.avg.ori{inside_idx(ii)} = ori(:,ii);
end