PQNexp_CRF(0,1,146,4);
PQNexp_CRF_plot(0,1,146,4);
xlim([1 50]);
fprintf('(paused)');
pause;

PQNexp_GGM(0,.1,10);
PQNexp_GGM_plot(0,.1,10);
xlim([1 50]);
fprintf('(paused)');
pause;

PQNexp_GrGGM(0,.1,10);
PQNexp_GrGGM_plot(0,.1,10);
xlim([1 15]);
fprintf('(paused)');
pause;

PQNexp_MRF(0,1000,5,10,'pseudo');
PQNexp_MRF_plot(0,1000,5,10,'pseudo');
fprintf('(paused)');
pause;

PQNexp_GrGGM2(0,1);
PQNexp_GrGGM2_plot(0,1);