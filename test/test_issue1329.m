function test_issue1329

% WALLTIME 00:60:00
% MEM 2gb
% DEPENDENCY ft_convert_coordsys

%%

[ftver, ftpath] = ft_version;

elecfile = fullfile(ftpath, 'template', 'electrode', 'standard_1020.elc');
elec = ft_read_sens(elecfile);
% the electrodes are in MNI coordinates, since these are placed on the template headmodel, which in turn is based on Colin27
% ft_determine_coordsys(elec)
elec.coordsys = 'mni';

mrifile = fullfile(ftpath, 'template', 'headmodel', 'standard_mri.mat');
mri = ft_read_mri(mrifile);
mri.coordsys = 'mni';

%%

method = 0;

% the specific ones are explicit in the orientation of the axes AND in the place of the origin
% the generic ones are explicit in the orientation of the axes BUT NOT in the place of the origin

specific = {'ctf' '4d' 'bti' 'eeglab' 'neuromag' 'itab' 'acpc' 'spm' 'mni' 'fsaverage' 'tal'};
generic  = {'als' 'ali' 'ars' 'ari' 'pls' 'pli' 'prs' 'pri' 'las' 'lai' 'ras' 'rai' 'lps' 'lpi' 'rps' 'rpi' 'asl' 'ail' 'asr' 'air' 'psl' 'pil' 'psr' 'pir' 'sal' 'ial' 'sar' 'iar' 'spl' 'ipl' 'spr' 'ipr' 'sla' 'ila' 'sra' 'ira' 'slp' 'ilp' 'srp' 'irp' 'lsa' 'lia' 'rsa' 'ria' 'lsp' 'lip' 'rsp' 'rip'};
coordsys = [specific generic];

%%
% convert MNI to anything, including TAL

for i=1:length(coordsys)
  disp(coordsys{i})
  dum = ft_convert_coordsys(elec, coordsys{i}, method);
  dum = ft_convert_coordsys(mri, coordsys{i}, method);
end

%%
% convert CTF to anything, note that ACPC is used instead of conversion to MNI, SPM and TAL

elec_ctf = ft_convert_coordsys(elec, 'ctf', method);
mri_ctf  = ft_convert_coordsys(mri, 'ctf', method);

for i=1:length(coordsys)
    disp(coordsys{i})
    dum = ft_convert_coordsys(elec_ctf, coordsys{i}, method);
    dum = ft_convert_coordsys(mri_ctf, coordsys{i}, method);
end

%%
% convert RAS to any other generic triplet

elec_ras = ft_convert_coordsys(elec, 'ras', method);
mri_ras  = ft_convert_coordsys(mri, 'ras', method);

for i=1:length(generic)
  disp(generic{i})
  dum = ft_convert_coordsys(elec, generic{i}, method);
  dum = ft_convert_coordsys(mri, generic{i}, method);
end
