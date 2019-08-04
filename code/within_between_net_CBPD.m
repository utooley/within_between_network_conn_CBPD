%Running on the cluster
datadir=fullfile('/data/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_nospkreg/')
listdir='/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
outdir='/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/Schaefer400zNetworks'

%running with the cluster mounted locally
% MODIFY THIS FOR DIFFERENT PIPELINES
datadir=fullfile('~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_spkreg_fd1.25dvars2_drpvls//')
listdir='/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
z_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_spkreg_fd1.25dvars2_drpvls/Schaefer200zNetworks'
noz_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_spkreg_fd1.25dvars2_drpvls/Schaefer200Networks'

%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n76_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_80119.csv'),'Delimiter',',','ReadVariableNames', 1)
%subjlist=readtable('/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n64_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_060419.csv', 'Delimiter',',')
%subjlist=subjlist(:,1:2);
%% Z-score FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.id0(n)) %look at this
    run=char(subjlist.id1(n))
    file=fullfile(datadir,strcat(sub,'/',run,'/fcon/schaefer200/',sub,'_',run,'_schaefer200_network.txt'));
    %file=fullfile(datadir,strcat(sub,'/fcon/schaefer400/',sub,'_schaefer400_network.txt'));
    outfile=fullfile(z_outdir, strcat(num2str(sub),'_',run,'_Schaefer200MNI_znetwork.txt'));
    %outfile=fullfile(z_outdir, strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'));
    if (exist(outfile)==2) %if it's already written don't do it again
       fprintf('Sub %s already exists. \n', sub);
    else
        subfcmat=load(file);
        %make into adjacency matrix and save out
        %size_vec=tril(ones(400,400),-1);
        size_vec=tril(ones(200,200),-1);
        adj_mat=size_vec;
        adj_mat(adj_mat==1)=subfcmat;
        subfcmat=adj_mat+adj_mat';
        outfile=fullfile(noz_outdir,strcat(sub,'_',run,'_schaefer200_network.txt'));
        %outfile=fullfile(noz_outdir,strcat(sub,'_schaefer400_network.txt'));
        csvwrite(outfile, subfcmat);
        %elimate parcel 52 (parcel index 103), delete row 103
        %subfcmat=removerows(subfcmat, 'ind', [103]);
        %remove column 103
       % subfcmat(:,103)=[]; %never checked parcel coverage for this.
        %replace the diagonal of 1's with 0's
        for x=1:200
            subfcmat(x,x)=0;
        end
        %create an empty z-matrx
        zfc=[];
        for i=1:200
            %cycle through each column of the FC matrix and do a fisher r-to-z
            %for each value
            zfc(:,i)=fisherz(subfcmat(:,i));
        end
        outfile=fullfile(z_outdir, strcat(num2str(sub),'_',run,'_Schaefer200MNI_znetwork.txt'));
        %outfile=fullfile(z_outdir, strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'));
        csvwrite(outfile, zfc);
    end
end

%% Within and between network connectivity
%cluster mounted locally
%modify pipeline here
clear modul
clear avgweight
clear num_comms_modul
clear system_segreg
clear mean_within_sys
clear mean_between_sys
clear system_conn
clear part_coef_pos
clear part_coef_neg
datadir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_spkreg_fd1.25dvars2_drpvls/Schaefer200zNetworks'
%yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer200/schaefer200x7CommunityAffiliation.1D')
system_segreg=zeros(height(subjlist),1);
mean_within_sys=zeros(height(subjlist),1);
mean_between_sys=zeros(height(subjlist),1);
part_coef_pos=zeros(height(subjlist),1);
part_coef_neg=zeros(height(subjlist),1);

for n=1:height(subjlist)
    sub=char(subjlist.id0(n)) %look at this
    run=char(subjlist.id1(n))
    file=fullfile(datadir,strcat(num2str(sub),'_',run,'_Schaefer200MNI_znetwork.txt'))
    %file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
    subfcmat = load(file);
    for x=1:200
        subfcmat(x,x)=0;
    end

%average whole system segregation, from Micaela Chan 2018
avgweight(n,1)=mean(subfcmat(subfcmat~=0));
[S, W, B] = segregation(subfcmat,yeo_nodes);
system_segreg(n,1)=S;
mean_within_sys(n,1)=W;
mean_between_sys(n,1)=B;

%Within and between connectivity, adapted from Micaela Chan 2018
Ci=yeo_nodes;
nCi = unique(Ci);
M=subfcmat;

for i = 1:length(nCi) % loop through communities
    for j = 1:length(nCi)
       Wi = Ci == nCi(i); % find index for within communitiy edges
       Bi = Ci == nCi(j); % find index for between communitiy edges to specific community
       
       Wv_temp = M(Wi,Wi); % extract within communitiy edges
       Bv_temp = M(Wi,Bi); % extract between communitiy edges to specific community
       
       %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
       Bv = [Bv_temp(:)'];
       system_connectivity(i,j)=mean(Bv(Bv~=0));
       %system_between(i,1)=mean(Bv);
       %if i==j
       %else
    end
end

%transpose these (or something so they can be saved out on a subject basis.
system_connectivity;
system_conn_vector = reshape(system_connectivity',[],1)';

system_conn(n,:)=system_conn_vector;

%modularity using Gen Louvain classic modularity maximization, negative
%asymmetric weighting
[M Q]=community_louvain(subfcmat, 1, [], 'negative_asym');
modul(n,1)=Q;
num_comms_modul(n,1)=length(unique(M));%also save the number of communities detected.

%average participation coefficient!
[Ppos, Pneg]=participation_coef_sign(subfcmat, yeo_nodes);
part_coef_pos(n,1)=mean(Ppos);
part_coef_neg(n,1)=mean(Pneg);

end
outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_spkreg_fd1.25dvars2_drpvls/'
header={'ID', 'run', 'avgweight', 'modul', 'num_comms_modul','part_coef_pos','part_coef_neg', 'system_segreg', 'mean_within_sys', 'mean_between_sys', 'sys1to1','sys1to2','sys1to3','sys1to4','sys1to5','sys1to6','sys1to7','sys2to1','sys2to2','sys2to3','sys2to4','sys2to5','sys2to6','sys2to7','sys3to1','sys3to2','sys3to3','sys3to4','sys3to5','sys3to6','sys3to7','sys4to1','sys4to2','sys4to3','sys4to4','sys4to5','sys4to6','sys4to7','sys5to1','sys5to2','sys5to3','sys5to4','sys5to5','sys5to6','sys5to7','sys6to1','sys6to2','sys6to3','sys6to4','sys6to5','sys6to6','sys6to7','sys7to1','sys7to2','sys7to3','sys7to4','sys7to5','sys7to6','sys7to7'}

outfile=table(char(subjlist.id0), char(subjlist.id1), avgweight, modul, num_comms_modul, part_coef_pos, part_coef_neg, system_segreg, mean_within_sys, mean_between_sys, system_conn)
outfile2=splitvars(outfile, 'system_conn')
outfile2.Properties.VariableNames=header

save(fullfile(outdir, 'n76_within_between_Yeo7_Schaefer200_withmodulpartcoef.mat'), 'outfile2')
writetable(outfile2,fullfile(outdir,'n76_within_between_Yeo7_Schaefer200_withmodulpartcoef.csv'))

