% $Id: lagextraction_robustness.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

close all; clear; clc;

USE_WOODY = false;
% USE_WOODY = true;

SNRs = 0.5:0.5:4;
% SNRs = 0.5:1:4;
% noise_amps = linspace(3,30,3);
% noise_amps = linspace(3,50,30);
% noise_amps = linspace(3,50,5);

N = 10;
err = zeros(length(SNRs),6,N);

for tt=1:N

    rseed = tt;

    template_id = 1;
    USE_AR = true;
    lagextraction_robustness_aux
    err(:,1,tt) = abs([cref(:) - cref_gt(:)]);

    USE_AR = false;
    lagextraction_robustness_aux
    err(:,2,tt) = abs([cref(:) - cref_gt(:)]);

    refs(:,1) = ref(:);

    template_id = 2;
    USE_AR = true;
    lagextraction_robustness_aux
    err(:,3,tt) = abs([cref(:) - cref_gt(:)]);

    USE_AR = false;
    lagextraction_robustness_aux
    err(:,4,tt) = abs([cref(:) - cref_gt(:)]);

    refs(:,2) = ref(:);

    template_id = 3;
    USE_AR = true;
    lagextraction_robustness_aux
    err(:,5,tt) = abs([cref(:) - cref_gt(:)]);

    USE_AR = false;
    lagextraction_robustness_aux
    err(:,6,tt) = abs([cref(:) - cref_gt(:)]);

    refs(:,3) = ref(:);

    eval(['save err_',num2str(tt),' err']);

end

err_mean = mean(err,3);
err_std = std(err,1,3);

smart_figure('Corr with ref'); clf; hold on
plot(10 * log10(SNRs(end)),err(end,2),'xk',10 * log10(SNRs(end)),err(end,2),'ok','Linewidth',2);
legend({'AR noise','White noise'})
plot(10 * log10(SNRs),err(:,1),'-xb',10 * log10(SNRs),err(:,2),'-ob','Linewidth',2);
plot(10 * log10(SNRs),err(:,3),'--xg',10 * log10(SNRs),err(:,4),'--og','Linewidth',2);
plot(10 * log10(SNRs),err(:,5),'-.xr',10 * log10(SNRs),err(:,6),'-.or','Linewidth',2);
xlabel('SNR (dB)')
ylabel('Error')
axis tight
ylim([0,max(err(:))])
grid on
hold off
savefig([pref,'realign_corr_ref'],22,{'png','pdf'});

smart_figure('Ref');
plot(1:size(refs,1),refs(:,1),'b-',1:size(refs,1),refs(:,2),'g--', ...
     1:size(refs,1),refs(:,3),'r-.','Linewidth',2);
xlabel('Time (ms)')
legend({'Template 1','Template 2','Template 3'},'Location','Southwest');
savefig([pref,'realign_templates'],22,{'png','pdf'});

% archive_name = [pref,'_results_'];
% eval(['save ',archive_name,' ref err']);
