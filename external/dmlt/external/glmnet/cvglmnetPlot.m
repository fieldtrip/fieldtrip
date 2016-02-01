function cvglmnetPlot(CVglmnet)
% Plot MSE w. errorbars for a cross validation fit of glmnet
% Call: cvglmnetPlot(CVglmnet)
%   CVglmnet: output from cvglmnet(...) function
    errorbar(log(CVglmnet.glmnetOptions.lambda),CVglmnet.cvm,CVglmnet.stderr,'linewidth',2)
    hold on
    plot(log(CVglmnet.glmnetOptions.lambda),CVglmnet.cvm,'r-o','linewidth',2)
    [val idx2] = min(fliplr(CVglmnet.cvm));
    vals=flipud(log(CVglmnet.glmnetOptions.lambda));
    plot(log([CVglmnet.lambda_min CVglmnet.lambda_min]),minmax([CVglmnet.cvlo CVglmnet.cvup]),'g--',...
        'linewidth',3)
    if CVglmnet.lambda_min ~=CVglmnet.lambda_1se
        plot(log([CVglmnet.lambda_1se CVglmnet.lambda_1se]),minmax([CVglmnet.cvlo CVglmnet.cvup]),'k--',...
            'linewidth',3)
    end
    xlabel('log(\lambda)','fontsize',12,'fontweight','bold')
    ylabel('Mean squared error +/- std. err','fontsize',12,'fontweight','bold')
    title('Nonzero elements','fontsize',12,'fontweight','bold') %use as x-top-label.
    axis tight
    H1=gca;
    ax=axis;
    H2=axes('position',get(H1,'position'));
    set(H2,'color','none')
    axis(ax);
    idx=floor(linspace(1,length(CVglmnet.glmnetOptions.lambda),7));
    labl=flipud(CVglmnet.glmnet_object.df);
    set(H2,'XAxisLocation','top')
    set(H2,'xtick',vals(idx))
    set(H2,'xticklabel',labl(idx))
    set(H2,'YAxisLocation','right')
    set(H2,'yticklabel','')
    set([H1 H2],'box','off','fontsize',12,'fontweight','bold')
    text(vals(idx2),val,['min, nz = ' num2str(labl(idx2))],'fontweight','bold','fontsize',14);
end