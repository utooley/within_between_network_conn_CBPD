% Example code for how to run the spin test
% Medial wall removal is now included
% SMW 07/31/2020

% Step 1: SpinPermuFs.m to obtain 'spins' of the data 
% (or use SpinPermuCIVET.m):

%Set up paths
fshome = "/Applications/freesurfer/";

% left and right surfaces (group-averaged at every vertex):
readleft = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_linear_on_CT_from_FS.csv'; 
readright = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_linear_on_CT_from_FS.csv'; 
permno = 1000; % how many spins
wsname = sprintf('~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/CT_spins_1000x');

%SpinPermuFS(readleft,readright,permno,wsname)

% Step 2: pvalvsNull.m to run the spin test
% left and right hemispheres for the second modality:
readleft1 = readleft;
readright1 = readright;
readleft2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_pos_sig_edges_0.05_rev.csv';
readright2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_pos_sig_edges_0.05_rev.csv';

% indicate (with 0's and 1's) which vertices in the left and right
% hemispheres are part of the medial wall
[vl, left_labels, ctl] = read_annotation(fullfile(fshome,'/subjects/fsaverage/label/lh.aparc.a2009s.annot'));
%v_exclude_left = left_labels==1644825; % label of vertices in the medial wall is 1644825
v_exclude_left = left_labels==0; % label of vertices in the medial wall is 0
[vr,right_labels,ctr] = read_annotation(fullfile(fshome,'/subjects/fsaverage/label/rh.aparc.a2009s.annot'));
%v_exclude_right = right_labels==1644825;
v_exclude_right = right_labels==0;

pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right)

%Pos at 0.01
readleft2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_pos_sig_edges_0.01_rev.csv';
readright2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_pos_sig_edges_0.01_rev.csv';
pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right)

%Pos at 0.001
readleft2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_pos_sig_edges_0.001_rev.csv';
readright2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_pos_sig_edges_0.001_rev.csv';
pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right)

%Pos at 0.0001
readleft2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_pos_sig_edges_0.00001_rev.csv';
readright2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_pos_sig_edges_0.00001_rev.csv';
pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right)

%% SA age-squared from Freesurfer

% left and right surfaces (group-averaged at every vertex):
readleft = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_sq_on_SA_from_FS.csv'; 
readright = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_sq_on_SA_from_FS.csv'; 
permno = 1000; % how many spins
wsname = sprintf('~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/SA_spins_1000x');

SpinPermuFS(readleft,readright,permno,wsname)

% Step 2: pvalvsNull.m to run the spin test
% left and right hemispheres for the second modality:
readleft1 = readleft;
readright1 = readright;
readleft2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_pos_sig_edges_0.05.csv';
readright2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_pos_sig_edges_0.05.csv';

% indicate (with 0's and 1's) which vertices in the left and right
% hemispheres are part of the medial wall
[vl, left_labels, ctl] = read_annotation(fullfile(fshome,'/subjects/fsaverage/label/lh.aparc.a2009s.annot'));
%v_exclude_left = left_labels==1644825; % label of vertices in the medial wall is 1644825
v_exclude_left = left_labels==0; % label of vertices in the medial wall is 0
[vr,right_labels,ctr] = read_annotation(fullfile(fshome,'/subjects/fsaverage/label/rh.aparc.a2009s.annot'));
%v_exclude_right = right_labels==1644825;
v_exclude_right = right_labels==0;

pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right)

%Neg at 0.01
readleft2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_pos_sig_edges_0.01.csv';
readright2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_pos_sig_edges_0.01.csv';
pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right)

%Neg at 0.001
readleft2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_pos_sig_edges_0.001.csv';
readright2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_pos_sig_edges_0.001.csv';
pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right)

%Neg at 0.0001
readleft2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/lh_age_pos_sig_edges_1e-04.csv';
readright2 = '~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/spin_tests/rh_age_pos_sig_edges_1e-04.csv';
pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right)

