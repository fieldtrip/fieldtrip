% Demonstrate how negative values arise from sampling variability when
% there is no effect

% trials per class
Ntrl_class = 100;
% sampling repetitions
Nrep = 1000;

Inobc = zeros(1,Nrep);
Ibc = zeros(1,Nrep);
for ri=1:Nrep
    % generate data
    % the class has no effect on the distribution the data is drawn from
    x = randn(2*Ntrl_class,1);
    cx = copnorm(x);
    y = [zeros(Ntrl_class,1); ones(Ntrl_class,1)];
    
    % without bias correction
    Inobc(ri) = mi_model_gd(cx,y,2,false,true);
    % with bias correction
    Ibc(ri) = mi_model_gd(cx,y,2,true,true);
end

%%
figure
subplot(1,2,1)
hist(Inobc,20);
vline(mean(Inobc),'r')
title('mi\_model\_gd, no bias correction')
set(gca,'FontSize',12)

subplot(1,2,2)
hist(Ibc,20);
vline(mean(Ibc),'r')
title('mi\_model\_gd, with bias correction')
set(gca,'FontSize',12)

