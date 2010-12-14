doZoom = 1;

PQNexp_CRF_plot(1,1,8936,19);
if doZoom
    xlim([1 500]);
    ylim([1.5e4 4.5e4]);
end

PQNexp_GGM_plot(1,.1,50);
if doZoom
    figure(2);
    xlim([1 100]);
    ylim([-1200 -400]);
    figure(3);
    xlim([1 1000]);
    ylim([-650 -450]);
end

PQNexp_GrGGM_plot(1,.1,50);
if doZoom
    figure(4);
    xlim([1 200]);
    ylim([-1200 -200]);
    figure(5);
    xlim([1 1000]);
    ylim([-420 -260]);
end

PQNexp_MRF_plot(1,5400,10,50,'exact');
if doZoom
    figure(6);
    xlim([1 100]);
    ylim([3.6e4 4.5e4]);
    figure(7);
    xlim([1 1000]);
    ylim([3.6e4 4.35e4]);
end

PQNexp_GrGGM2_plot(1,50);
if doZoom
    figure(8);
    xlim([1e-4 1.25e-2]);
    ylim([-543 -530]);
end