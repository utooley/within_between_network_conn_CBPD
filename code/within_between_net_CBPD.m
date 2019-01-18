%Running on the cluster
datadir=fullfile('/data/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_nospkreg/')
listdir='/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
outdir='/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/Schaefer400zNetworks'

%running with the cluster mounted locally
datadir=fullfile('~/Desktop/cluster/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_nospkreg/')
listdir='~/Desktop/cluster/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
z_outdir='~/Desktop/cluster/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/Schaefer400zNetworks'
noz_outdir='~/Desktop/cluster/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/Schaefer400Networks'


%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n47_cohort_usable_t1_rest_1mm_outliers_10_2mm_11718.csv'),'Delimiter',',')
subjlist=subjlist(:,1);
%% Z-score FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.id0(n)) %look at this
    file=fullfile(datadir,strcat(sub,'/fcon/schaefer400/',sub,'_schaefer400_network.txt'));
    subfcmat=load(file);
    %make into adjacency matrix and save out
    size_vec=tril(ones(400,400),-1);
    adj_mat=size_vec;
    adj_mat(adj_mat==1)=subfcmat;
    subfcmat=adj_mat+adj_mat';
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
    outfile=fullfile(z_outdir, strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'));
    csvwrite(outfile, zfc);
end

%% Within and between network connectivity
%cluster mounted locally
datadir='~/Desktop/cluster/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/Schaefer400zNetworks'
yeo_nodes=dlmread('~/Desktop/cluster/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
system_segreg=zeros(height(subjlist),1);
mean_within_sys=zeros(height(subjlist),1);
mean_between_sys=zeros(height(subjlist),1);

for n=1:height(subjlist)
    sub=char(subjlist.id0(n)) %look at this
    file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
    subfcmat = load(file);
    for x=1:359
        subfcmat(x,x)=0;
    end

%average whole system segregation, from Micaela Chan 2018
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

end
outdir='~/Desktop/cluster/Ursula/projects/in_progress/within_between_network_conn_CBPD/output/data'
outfile=dataset(char(subjlist.id0), system_segreg, mean_within_sys, mean_between_sys, system_conn)
%this isn't working at the moment-figure out header.
header={'subjlist', 'system_segreg', 'mean_within_sys', 'mean_between_sys', '1to1','1to2','1to3','1to4','1to5','1to6','1to7','2to1','2to2','2to3','2to4','2to5','2to6','2to7','3to1','3to2','3to3','3to4','3to5','3to6','3to7','4to1','4to2','4to3','4to4','4to5','4to6','4to7','5to1','5to2','5to3','5to4','5to5','5to6','5to7','6to1','6to2','6to3','6to4','6to5','6to6','6to7','7to1','7to2','7to3','7to4','7to5','7to6','7to7'}
 

export(outfile,'File',fullfile(outdir,'n47_within_between_Yeo7_Schaefer400.csv'),'Delimiter',',')

