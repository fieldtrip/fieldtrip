function test_ctf2grad

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_read_header ft_read_sens ctf2grad

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf151');
dataset = fullfile(datadir, 'Subject01.ds');
coildeffile = fullfile(datadir, 'ctf151dccn_coil_def.dat');
grad1 = ft_read_sens(dataset, 'senstype', 'meg', 'coordsys', 'dewar');
grad1 = ft_convert_units(grad1,'m');
grad2 = ft_read_sens(dataset, 'senstype', 'meg', 'coordsys', 'dewar', 'coilaccuracy', 0);
grad3 = ft_read_sens(dataset, 'senstype', 'meg', 'coordsys', 'dewar', 'coilaccuracy', 0, 'coildeffile', coildeffile);
% the order of the MEGGRAD-REF channels is different, as is the orientation of one of
% the coils per gradiometer. This should be reflected by the tra



coillist1ref = 303:356;
coillist2ref = 1:54;
coillist1meg = 1:302; coillist1meg = reshape(coillist1meg,151,2)'; coillist1meg = coillist1meg(:);
coillist2meg = 55:356;

chanlist1ref = 152:184;
chanlist2ref = 1:33;
chanlist1meg = 1:151;
chanlist2meg = 34:184;

% compare positions -> this only gets as good is it gets w.r.t. the
% accuracy of the coil_def info with regard to the system under
% consideration
deltameg = sqrt(sum((grad1.coilpos(coillist1meg,:)-grad2.coilpos(coillist2meg,:)).^2,2));
deltameg13 = sqrt(sum((grad1.coilpos(coillist1meg,:)-grad3.coilpos(coillist2meg,:)).^2,2));

deltaref = sqrt(sum((grad1.coilpos(coillist1ref,:)-grad2.coilpos(coillist2ref,:)).^2,2));
deltaref13 = sqrt(sum((grad1.coilpos(coillist1ref,:)-grad3.coilpos(coillist2ref,:)).^2,2));

sumtra = sum(grad2.tra)';
orimeg = grad1.coilori(coillist1meg,:)-grad2.coilori(coillist2meg,:).*sumtra(coillist2meg);
oriref = grad1.coilori(coillist1ref,:)-grad2.coilori(coillist2ref,:).*sumtra(coillist2ref);



cp1=grad1.coilpos(303:356,:); % after 151*2 coils
co1=grad1.coilori(303:356,:);
t1=grad1.tra(152:end,303:356);

cp2=grad2.coilpos(1:54,:); % take the first 356-151*2 coils
co2=grad2.coilori(1:54,:);
t2=grad2.tra(1:33,1:54);

[i1,i2] = match_str(grad1.label,grad2.label);

figure;imagesc(grad1.tra)
figure;imagesc(grad2.tra)
% grad1 has the ref sensors + coils at the end, grad2 at the beginning


figure;ft_plot_sens(grad1, 'chantype', 'refgrad', 'coil', true);
figure;ft_plot_sens(grad2, 'chantype', 'refgrad', 'coil', true);

% in head coordinates
grad1b = ft_read_sens('Subject01.ds', 'senstype', 'meg');
grad1b = ft_convert_units(grad1b,'m');
grad2b = ft_read_sens('Subject01.ds', 'senstype', 'meg', 'coilaccuracy', 0);
grad3b = ft_read_sens('Subject01.ds', 'senstype', 'meg', 'coilaccuracy', 0, 'coildeffile', coildeffile);

figure;ft_plot_sens(grad1b, 'chantype', 'refgrad', 'coil', true);
figure;ft_plot_sens(grad2b, 'chantype', 'refgrad', 'coil', true);

headmodel = ft_read_headmodel(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.hdm'));
headmodel = ft_convert_units(headmodel,'m');
sourcemodel.pos = [0 0 0.1];

cfg = [];
cfg.grad = grad1b;
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel;
lf1 = ft_prepare_leadfield(cfg);
cfg.grad = grad2b;
lf2 = ft_prepare_leadfield(cfg);
cfg.grad = grad3b;
lf3 = ft_prepare_leadfield(cfg);

figure;plot(lf1.leadfield{1}(i1,:));
figure;plot(lf2.leadfield{1}(i2,:));

figure;plot(lf1.leadfield{1}(i1,:)-lf2.leadfield{1}(i2,:));
figure;plot(lf1.leadfield{1}(i1,:)-lf3.leadfield{1}(i2,:));
