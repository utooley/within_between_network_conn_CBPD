%Running on the cluster
datadir=fullfile('/data/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_nospkreg/')
listdir='/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
outdir='/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/Schaefer400zNetworks'

%running with the cluster mounted locally
% MODIFY THIS FOR DIFFERENT PIPELINES
datadir=fullfile('~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_nospkreg_nodespike_bestrunonly/')
listdir='/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
z_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_nospkreg_nodespike_bestrunonly/Schaefer400zNetworks'
noz_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_nospkreg_nodespike_bestrunonly/Schaefer400Networks'

%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n64_total_subjs_usable_t1_rest_1mm_outliers_10_2mm_060419'),'Delimiter',',','ReadVariableNames', 0)
%subjlist=readtable('/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n64_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_060419.csv', 'Delimiter',',')
%subjlist=subjlist(:,1:2);
%% Z-score FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.Var1(n)) %look at this
    %run=char(subjlist.id1(n))
    %file=fullfile(datadir,strcat(sub,'/',run,'/fcon/schaefer400/',sub,'_',run,'_schaefer400_network.txt'));
    file=fullfile(datadir,strcat(sub,'/fcon/schaefer400/',sub,'_schaefer400_network.txt'));
    %outfile=fullfile(z_outdir, strcat(num2str(sub),'_',run,'_Schaefer400MNI_znetwork.txt'));
    outfile=fullfile(z_outdir, strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'));
    if (exist(outfile)==2) %if it's already written don't do it again
       fprintf('Sub %s already exists. \n', sub);
    else
        subfcmat=load(file);
        %make into adjacency matrix and save out
        size_vec=tril(ones(400,400),-1);
        adj_mat=size_vec;
        adj_mat(adj_mat==1)=subfcmat;
        subfcmat=adj_mat+adj_mat';
        %outfile=fullfile(noz_outdir,strcat(sub,'_',run,'_schaefer400_network.txt'));
        outfile=fullfile(noz_outdir,strcat(sub,'_schaefer400_network.txt'));
        csvwrite(outfile, subfcmat);
        %elimate parcel 52 (parcel index 103), delete row 103
        %subfcmat=removerows(subfcmat, 'ind', [103]);
        %remove column 103
       % subfcmat(:,103)=[]; %never checked parcel coverage for this.
        %replace the diagonal of 1's with 0's
        for x=1:359
            subfcmat(x,x)=0;
        end
        %create an empty z-matrx
        zfc=[];
        for i=1:400
            %cycle through each column of the FC matrix and do a fisher r-to-z
            %for each value
            zfc(:,i)=fisherz(subfcmat(:,i));
        end
        %outfile=fullfile(z_outdir, strcat(num2str(sub),'_',run,'_Schaefer400MNI_znetwork.txt'));
        outfile=fullfile(z_outdir, strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'));
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
datadir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_nospkreg_nodespike_bestrunonly/Schaefer400zNetworks'
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
system_segreg=zeros(height(subjlist),1);
mean_within_sys=zeros(height(subjlist),1);
mean_between_sys=zeros(height(subjlist),1);

for n=1:height(subjlist)
    sub=char(subjlist.Var1(n)) %look at this
    %run=char(subjlist.id1(n))
    %file=fullfile(datadir,strcat(num2str(sub),'_',run,'_Schaefer400MNI_znetwork.txt'))
    file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
    subfcmat = load(file);
    for x=1:359
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

end
outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_nospkreg_nodespike_bestrunonly/'
outfile=dataset(char(subjlist.Var1), avgweight, modul, num_comms_modul, system_segreg, mean_within_sys, mean_between_sys, system_conn)
%this isn't working at the moment-figure out header.
header={'subjlist', 'system_segreg', 'mean_within_sys', 'mean_between_sys', '1to1','1to2','1to3','1to4','1to5','1to6','1to7','2to1','2to2','2to3','2to4','2to5','2to6','2to7','3to1','3to2','3to3','3to4','3to5','3to6','3to7','4to1','4to2','4to3','4to4','4to5','4to6','4to7','5to1','5to2','5to3','5to4','5to5','5to6','5to7','6to1','6to2','6to3','6to4','6to5','6to6','6to7','7to1','7to2','7to3','7to4','7to5','7to6','7to7'}
 
save(fullfile(outdir, 'n64_within_between_Yeo7_Schaefer400_withmodul.mat'), 'outfile')
export(outfile,'File',fullfile(outdir,'n64_within_between_Yeo7_Schaefer400_withmodul.csv'),'Delimiter',',')

