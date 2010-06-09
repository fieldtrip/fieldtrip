% SPM2 (c) 1991,1994-2003
% Statistical Parametric Mapping           -                       SPM2 
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  
% \__ \ )___/ )    (   Statistical Parametric Mapping
% (___/(__)  (_/\/\_)  SPM - http://www.fil.ion.ucl.ac.uk/spm
%_______________________________________________________________________
%
% This Contents.m file holds the version ID for this release of Matlab,
% and contains a manifest of the included functions and their version numbers.
%
% SPM2 is written for Matlab v6.0.0 under UNIX and Windows
% ( Compiled binaries of external MEX functions are provided for:       )
% (                   Solaris2, Linux, and Windows                      )
%
% See spm.man for details of this release.
% See the README for information on installation and getting started.
% See spm_motd.man for last minute release details.
%
% SPM (being the collection of files given in the manifest below) is
% free but copyright software, distributed under the terms of the GNU
% General Public Licence as published by the Free Software Foundation
% (either version 2, as given in file spm_LICENCE.man, or at your option,
% any later version). Further details on "copyleft" can be found at
% http://www.gnu.org/copyleft/.
%
%_______________________________________________________________________
% @(#)Contents.m	2.11 03/05/12
%
% SPM2 - Manifest
%-----------------------------------------------------------------------

help Contents

%=======================================================================
% PROGRAMMERS NOTE:
% This (Contents.m) is the contents file for SPM, used by spm('Ver') to
% recover the version number and copyright information. MatLab's ver
% also uses Contents.m files to identify toolbox versions.
% Line1: Version (first word) & copyright information (rest of line).
% Line2: One line description
%
% This file was produced with the help of:
%      what Makefile *.{m,c,h,man} */*/*.{m,c,h,man} | nawk '{if (NF>1) printf("%c %-30s %s\n", 37, $1, $2)}' >> Contents.m
%=======================================================================
% Makefile                       2.16
% spm_image.m                    2.20
% spm_load.m                     2.1
% spm_P_RF.m                     2.12
% spm_reml.m                     2.22
% spm.man                        2.15
% spm_mip.m                      2.7
% spm_get_ons.m                  2.39
% spm_contrasts.m                2.3
% spm_input.m                    2.8
% spm_input.m                    2.8
% spm_progress_bar.m             2.2
% spm_sections.m                 2.14
% win32mmap.c                    2.1
% connex.h                       1.2
% spm.m                          2.87
% spm_P.m                        2.3
% spm_powell.m                   2.3
% spm_add.c                      2.10
% spm_add.m                      2.7
% spm_resss.m                    2.9
% spm_Ce.m                       2.7
% spm_coreg_ui.m                 2.9
% spm_meanby.m                   2.1
% spm_clusters.c                 2.2
% spm_help.m                     2.47
% spm_max.c                      2.4
% spm_PEB.m                      2.11
% spm_graph.m                    2.40
% spm_list.m                     2.43
% spm_nlsi_GN.m                  2.9
% spm_project.c                  2.4
% spm_Q.m                        2.2
% spm_int.m                      2.4
% spm_list_files.c               2.10
% spm_list_files_expand.m        1.5
% spm_nlsi.m                     2.6
% spm_print.m                    1.3
% spm_bi_reduce.m                2.5
% spm_brainwarp.m                2.2
% spm_close_vol.m                2.4
% spm_en.m                       1.1
% spm_mean_ui.m                  2.5
% spm_unlink.m                   2.1
% spm_platform.m                 2.11
% spm_smooth.m                   2.3
% spm_spm.m                      2.66
% spm_t2z.m                      2.2
% spm_Tpdf.m                     2.2
% spm_vol_utils.c                2.5
% spm_brainwarp.c                2.12
% spm_diff.m                     2.1
% spm_get_space.m                2.13
% spm_regions.m                  2.17
% spm_unlink.c                   1.2
% spm_dicom_headers.m            2.13
% spm_global.c                   2.4
% spm_mask.m                     2.12
% spm_mip_ui.m                   2.12
% spm_coreg.m                    2.4
% spm_defaults.m                 2.23
% spm_dicom_convert.m            2.12
% spm_FcUtil.m                   2.15
% spm_getSPM.m                   2.50
% spm_hrf.m                      2.8
% spm_log.m                      2.4
% spm_mfx.m                      2.10
% spm_Ppdf.m                     2.2
% spm_render.m                   2.19
% spm_slice_vol.c                2.1
% spm_u.m                        2.2
% spm_VOI.m                      2.18
% spm_vol_ecat7.m                2.11
% spm_max.m                      2.2
% spm_dcm_priors.m               2.6
% spm_imcalc_ui.m                2.7
% spm_sample_vol.c               2.1
% spm_str_manip.m                2.10
% spm_vol_utils.h                2.2
% spm_detrend.m                  2.1
% spm_smooth_ui.m                2.11
% spm_atranspa.m                 2.1
% spm_bsplinc.c                  2.4
% spm_dctmtx.m                   2.1
% spm_Pcdf.m                     2.2
% spm_reslice.m                  2.12
% spm_conv_vol.c                 2.4
% spm_fMRI_design_show.m         2.22
% spm_sptop.m                    1.7
% spm_XYZreg.m                   2.6
% spm_bsplins.c                  2.5
% spm_motd.man                   2.8
% spm_vol_ana.m                  2.6
% spm_write_sn.m                 2.17
% spm_dcm_display.m              2.5
% spm_clf.m                      2.1
% spm_get.m                      2.43
% spm_resels_vol.c               2.4
% spm_resels_vol.m               2.2
% spm_adjmean_ui.m               2.8
% spm_barh.m                     2.2
% spm_bias_ui.m                  2.4
% spm_datatypes.h                2.2
% spm_DesMtx.m                   2.10
% spm_Fcdf.m                     2.2
% spm_Fpdf.m                     2.2
% spm_Gcdf.m                     2.2
% spm_Gpdf.m                     2.2
% spm_imatrix.m                  2.1
% spm_invGcdf.m                  2.2
% spm_invTcdf.m                  2.2
% spm_invXcdf.m                  2.3
% spm_Ncdf.m                     2.2
% spm_normalise.m                2.8
% spm_normalise.m                2.8
% spm_normalise.m                2.8
% spm_Npdf.m                     2.2
% spm_sp.m                       2.14
% spm_surf.m                     2.3
% spm_Tcdf.m                     2.2
% spm_type.m                     2.3
% spm_write_filtered.m           2.7
% spm_Xcdf.m                     2.2
% spm_Xpdf.m                     2.2
% spm_clusters.m                 2.2
% spm_conman.m                   2.18
% spm_conv.m                     1.5
% spm_conv_vol.m                 2.1
% spm_DesRep.m                   2.31
% spm_invFcdf.m                  2.2
% spm_matx.m                     2.5
% spm_P_FDR.m                    2.4
% spm_SpUtil.m                   2.18
% spm_invNcdf.m                  2.2
% spm_RandFX.man                 2.2
% spm_global.m                   2.2
% spm_grid.m                     2.2
% spm_resels.m                   2.2
% spm_adjmean_fmri_ui.m          2.9
% spm_bsplinc.m                  2.3
% spm_nCr.m                      2.1
% spm_spm_ui.m                   2.49
% spm_spm_ui.m                   2.49
% spm_sys_deps.h                 2.6
% spm_Volt_W.m                   1.2
% spm_hdm_ui.m                   2.15
% spm_list_files.m               2.1
% spm_hist2.c                    2.9
% spm_matrix.m                   2.1
% spm_project.m                  2.2
% spm_render_vol.m               2.1
% spm_sample_vol.m               2.1
% spm_slice_vol.m                2.1
% spm_spm_Bayes.m                2.7
% spm_uw_apply.m                 1.7
% spm_format.man                 2.2
% spm_filter.m                   2.10
% spm_kernels.m                  2.2
% spm_read_hdr.m                 2.1
% spm_figure.m                   2.37
% spm_render_vol.c               2.2
% spm_atranspa.c                 2.1
% spm_getdata.c                  2.3
% spm_segment.m                  2.24
% spm_make_lookup.c              2.3
% spm_mapping.c                  2.7
% spm_normalise_ui.m             2.9
% spm_Pec_resels.m               2.1
% spm_Bcdf.m                     2.2
% spm_Bpdf.m                     2.2
% spm_Icdf.m                     2.2
% spm_Ipdf.m                     2.2
% spm_bsplins.m                  2.3
% spm_results_ui.m               2.42
% spm_vol_access.c               2.2
% spm_vol_minc.m                 2.16
% spm_write_vol.m                2.9
% spm_check_registration.m       2.5
% spm_chi2_plot.m                1.2
% spm_peb_ppi.m                  2.10
% spm_realign.m                  2.36
% spm_XYZreg_Ex2.m               2.1
% spm_dcm_ui.m                   2.9
% spm_fmri_spm_ui.m              2.51
% spm_imcalc.m                   2.6
% spm_est_smoothness.m           2.6
% spm_hdm_priors.m               2.6
% spm_invBcdf.m                  2.2
% spm_invIcdf.m                  2.2
% spm_invPcdf.m                  2.2
% spm_LICENCE.man                2.4
% spm_uc_FDR.m                   2.6
% spm_hist2.m                    2.3
% spm_mvNpdf.m                   2.1
% spm_orthviews.m                2.38
% spm_krutil.c                   2.2
% spm_write_plane.m              2.17
% spm_bias_apply.m               2.4
% spm_getdata.h                  2.3
% spm_dilate.c                   1.2
% spm_read_netcdf.m              2.1
% spm_read_vols.m                2.5
% spm_make_lookup.h              2.2
% spm_bias_estimate.m            2.5
% spm_bias_mex.c                 2.5
% spm_uw_estimate.m              1.4
% spm_mapping.h                  2.2
% spm_vol_access.h               2.2
% spm_vol_check.m                1.2
% spm_smoothto8bit.m             2.2
% spm_get_bf.m                   2.22
% spm_reml_ancova.m              2.1
% spm_flip_analyze_images.m      2.4
% spm_svd.m                      2.2
% spm_vol.m                      2.15
% spm_win32utils.c               2.3
% spm_normalise_disp.m           2.5
% spm_segment_ui.m               2.4
% spm_win32utils.m               2.1
% spm_lambda_HRF.m               2.3
% spm_lx_dcm.m                   2.3
% spm_fx_dcm.m                   2.2
% spm_fx_HRF.m                   2.2
% spm_bilinear.m                 2.2
% spm_ancova.m                   2.1
% spm_orth.m                     2.1
% spm_slice_timing.m             2.16
% spm_templates.man              2.11
% spm_matfuns.c                  2.1
% spm_Volterra.m                 2.3
% spm_expm.m                     2.3
% spm_uc_RF.m                    2.4
% spm_P_Bonf.m                   2.3
% spm_uc.m                       2.2
% spm_uc_Bonf.m                  2.2
% spm_create_vol.m               2.13
% spm_transverse.m               2.21
% spm_setup_satfig.m             2.1
% spm_defaults_edit.m            2.25
% spm_non_sphericity.m           2.2
% spm_affreg.m                   2.3
% spm_dcm_reshape.m              2.2
% spm_get_def.m                  1.1
% spm_get_image_def.m            1.1
% spm_realign_ui.m               2.13
% spm_uw_show.m                  1.2
% spm_get_data.m                 2.1
% spm_fMRI_design.m              2.34
% spm_invdef_ui.m                1.2
% spm_sn2def.m                   1.3
% spm_applydef_ui.m              1.3
% spm_warp_ui.m                  1.3
% spm_warp.c                     1.1
% spm_combdef_ui.m               1.2
% spm_def2det_ui.m               1.2
% spm_procrustes_ui.m            1.2
% spm_procrustes_ui.m            1.2
% spm_affdef.c                   1.1
% spm_def2det.c                  1.1
% spm_Deformations.m             1.2
% spm_invdef.c                   1.1
% spm_invdef_vox.c               1.1
% spm_load_float.m               1.1
% spm_DICOM.m                    1.1
