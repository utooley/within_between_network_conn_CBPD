%Running on the cluster
datadir=fullfile('/data/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_nospkreg/')
listdir='/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
outdir='/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/Schaefer400zNetworks'
%subjlist=readtable('/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n64_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_060419.csv', 'Delimiter',',')
%subjlist=subjlist(:,1:2);

%% For each parcellation and each pipeline
parcellations={'schaefer200', 'schaefer400'}
pipelines={'nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls','gsr_censor_5contig_fd1.25dvars2_drpvls','nogsr_spkreg_fd1.25dvars2_drpvls','gsr_censor_5contig_fd0.5dvars1.75_drpvls'}
for p=1:length(parcellations)
    for pl=1:length(pipelines)
        parcellation=parcellations{p}
        pipeline=pipelines{pl}
        dim=str2double(parcellation(end-2:end))%what is the dimensionality of the parcellation
        
%running with the cluster mounted locally
datadir=strcat('~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_',pipeline)
listdir='/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
z_outdir=fullfile('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,strcat(parcellation,'zNetworks'))
noz_outdir=fullfile('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,strcat(parcellation,'Networks'))

%make directories if they don't exist
if ~exist(z_outdir, 'dir')
       mkdir(z_outdir)
end
if ~exist(noz_outdir, 'dir')
       mkdir(noz_outdir)
end
%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n74_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_80119.csv'),'Delimiter',',','ReadVariableNames', 1)
%% Z-score FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.id0(n)) %look at this
    run=char(subjlist.id1(n))
    file=fullfile(datadir,strcat(sub,'/',run,'/fcon/',parcellation,'/',sub,'_',run,'_',parcellation,'_network.txt'));
    %file=fullfile(datadir,strcat(sub,'/fcon/schaefer400/',sub,'_schaefer400_network.txt'));
    outfile=fullfile(z_outdir, strcat(num2str(sub),'_',run,'_',parcellation,'MNI_znetwork.txt'));
    %outfile=fullfile(z_outdir, strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'));
    if (exist(outfile)==2) %if it's already written don't do it again
       fprintf('Sub %s already exists. \n', sub);
    else
        subfcmat=load(file);
        %make into adjacency matrix and save out
        %size_vec=tril(ones(400,400),-1);
        size_vec=tril(ones(dim,dim),-1);
        adj_mat=size_vec;
        adj_mat(adj_mat==1)=subfcmat;
        subfcmat=adj_mat+adj_mat';
        outfile=fullfile(noz_outdir,strcat(sub,'_',run,'_',parcellation,'_network.txt'));
        %outfile=fullfile(noz_outdir,strcat(sub,'_schaefer400_network.txt'));
        csvwrite(outfile, subfcmat);
        %elimate parcel 52 (parcel index 103), delete row 103
        %subfcmat=removerows(subfcmat, 'ind', [103]);
        %remove column 103
       % subfcmat(:,103)=[]; %never checked parcel coverage for this.
        %replace the diagonal of 1's with 0's
        for x=1:dim
            subfcmat(x,x)=0;
        end
        %create an empty z-matrx
        zfc=[];
        for i=1:dim
            %cycle through each column of the FC matrix and do a fisher r-to-z
            %for each value
            zfc(:,i)=fisherz(subfcmat(:,i));
        end
        outfile=fullfile(z_outdir, strcat(num2str(sub),'_',run,'_',parcellation,'MNI_znetwork.txt'));
        %outfile=fullfile(z_outdir, strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'));
        csvwrite(outfile, zfc);
     end
end

%% Average those who have more than one run
datadir=strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,'/',parcellation,'zNetworks')
z_avg_outdir=fullfile('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,strcat(parcellation,'zNetworks_avg'))

if ~exist(z_avg_outdir, 'dir')
       mkdir(z_avg_outdir)
end
% JUST USE THE ALREADY MUTATED IN R DATA TO WEIGHT MATRICES, BUT NEED TO
% FIGURE OUT HOW TO DO THIS IN MATLAB
%qavars2=readtable('~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/mriqc_fd_2_mm/group_bold_with_names.csv','Delimiter',',','ReadVariableNames', 1)
%qavars=readtable(fullfile(datadir,'/XCP_QAVARS_FIXED_n74.csv'),'Delimiter',',','ReadVariableNames', 1)
myqavars=readtable(strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,'/n74_fixed_within_between_Yeo7_', parcellation, '_', pipeline, '_QA.csv'),'Delimiter',',','ReadVariableNames', 1)
dim=str2double(parcellation(end-2:end))%what is the dimensionality of the parcellation
[unique_subs, index]=unique(subjlist.id0)
unique_subjlist= subjlist(index,:)

for n=1:height(unique_subjlist)
    clear aggregmat
        sub=char(unique_subjlist.id0(n)) %look at this
        for r=1:4
            try
                run=strcat('run-0',num2str(r))
                file=fullfile(datadir,strcat(num2str(sub),'_',run,'_',parcellation,'MNI_znetwork.txt'))
                %file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
                subfcmat = load(file);
                for x=1:dim
                    subfcmat(x,x)=0;
                end
                rows = (string(myqavars.ID)==sub & string(myqavars.run)==run);
                vars = {'nVolCensored', 'size_t', 'totalnVolCensored','totalSizet'};
                weights=myqavars(rows, vars);
                %weight by
                %(sizet-nvolscensored)/(totalsizet-totalnvolscensored)
                weighted=subfcmat.*((weights.size_t-weights.nVolCensored)/(weights.totalSizet-weights.totalnVolCensored));
                aggregmat(:,:, r)=weighted;
            catch
                fprintf('There is no run %s for sub %s \n', num2str(r), sub);
            end
        end
        aggregmat(aggregmat == 0) = NaN;    %skip any matrices that are empty.     
meanMatrix = sum(aggregmat,3, 'omitnan'); %doc mean for more info.
%make the diagonal 0's
    for x=1:dim
        meanMatrix(x,x)=0;
    end
outfile=fullfile(z_avg_outdir, strcat(num2str(sub),'_',parcellation,'MNI_zavgnetwork.txt'));
csvwrite(outfile, meanMatrix);
end
    
%% Within and between network connectivity
%cluster mounted locally
%mdo this on the averaged weighted networks for people who have more than
%one.
clear modul
clear avgweight
clear num_comms_modul
clear system_segreg
clear mean_within_sys
clear mean_between_sys
clear system_conn
clear part_coef_pos
clear part_coef_neg
clear avgclustco_both
clear avgclustco_all
clear part_coef_avg_all
datadir=strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,'/',parcellation,'zNetworks_avg')
yeo_nodes=dlmread(fullfile('~/Desktop/cluster/picsl/mackey_group/tools',parcellation,strcat(parcellation,'x7CommunityAffiliation.1D.txt')))
% yeo_nodes=readtable(strcat('~/Desktop/cluster/picsl/mackey_group/tools/',parcellation,'/Schaefer2018_',num2str(dim),'Parcels_7Networks_order_comm.txt'))
% yeo_nodes.Properties.VariableNames={'index','community'}
% yeo_nodes.community=categorical(yeo_nodes.community)
%yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer200/schaefer200x7CommunityAffiliation.1D')
system_segreg=zeros(height(unique_subjlist),1);
mean_within_sys=zeros(height(unique_subjlist),1);
mean_between_sys=zeros(height(unique_subjlist),1);
part_coef_pos=zeros(height(unique_subjlist),1);
part_coef_neg=zeros(height(unique_subjlist),1);
avgclustco_both=zeros(height(unique_subjlist),1);
avgclustco_all=zeros(height(unique_subjlist),dim);
part_coef_avg_all=zeros(height(unique_subjlist),dim);

for n=1:height(unique_subjlist)
    sub=char(unique_subjlist.id0(n)) %look at this
    run=char(unique_subjlist.id1(n))
    file=fullfile(datadir,strcat(num2str(sub),'_',parcellation,'MNI_zavgnetwork.txt'))
    %file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
    subfcmat = load(file);
    for x=1:dim
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
%put this all into a matrix for everyone so that we can see the average
%system connectivity
system_connectivity_all(:,:,n)=system_connectivity;
system_conn_vector = reshape(system_connectivity',[],1)';

system_conn(n,:)=system_conn_vector;

%modularity using Gen Louvain classic modularity maximization, negative
%asymmetric weighting
for c = 1:100 %run it 100x for each subject, get an average number of communities detected.
    [M Q]=community_louvain(subfcmat, 1, [], 'negative_asym');
    modul_temp(c)=Q;
    num_comms_temp(c)=length(unique(M));
%also save the number of communities detected.
end
modul(n,1)=mean(modul_temp(:));
num_comms_modul(n,1)=mean(num_comms_temp(:));
%average participation coefficient!
[Ppos, Pneg]=participation_coef_sign(subfcmat, yeo_nodes);
part_coef_pos(n,1)=mean(Ppos);
part_coef_neg(n,1)=mean(Pneg);
%write out all nodes participation coefficient
part_coef_avg_all(n,:)=((Ppos+Pneg)/2);

%average clustering coefficient
avgclustco_both(n,1)=mean(clustering_coef_wu_sign(subfcmat,3));
%write out all nodes clustering coefficient
avgclustco_all(n,:)=clustering_coef_wu_sign(subfcmat,3);

end
outdir=strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline)
header={'ID', 'avgweight', 'modul_avg', 'avgclustco_both','num_comms_modul_avg','part_coef_pos','part_coef_neg', 'system_segreg', 'mean_within_sys', 'mean_between_sys', 'sys1to1','sys1to2','sys1to3','sys1to4','sys1to5','sys1to6','sys1to7','sys2to1','sys2to2','sys2to3','sys2to4','sys2to5','sys2to6','sys2to7','sys3to1','sys3to2','sys3to3','sys3to4','sys3to5','sys3to6','sys3to7','sys4to1','sys4to2','sys4to3','sys4to4','sys4to5','sys4to6','sys4to7','sys5to1','sys5to2','sys5to3','sys5to4','sys5to5','sys5to6','sys5to7','sys6to1','sys6to2','sys6to3','sys6to4','sys6to5','sys6to6','sys6to7','sys7to1','sys7to2','sys7to3','sys7to4','sys7to5','sys7to6','sys7to7'}

outfile=table(char(unique_subjlist.id0), avgweight, modul, avgclustco_both, num_comms_modul, part_coef_pos, part_coef_neg, system_segreg, mean_within_sys, mean_between_sys, system_conn)
outfile2=splitvars(outfile, 'system_conn')
outfile2.Properties.VariableNames=header

save(fullfile(outdir, strcat('n74_within_between_Yeo7_avgruns_',parcellation,'_withmodulpartcoef.mat')), 'outfile2')
writetable(outfile2,fullfile(outdir,strcat('n74_fixed_within_between_Yeo7_avgruns_',parcellation,'_withmodulpartcoef.csv')))

%also save the mean system connectivity matrix
mean_system_conn_mat=mean(system_connectivity_all,3)
% header={'sys1', 'sys2', 'sys3','sys4','sys5','sys6','sys7'}
% mean_system_conn_mat.Properties.VariableNames=header;
save(fullfile(outdir, strcat('n74_mean_system_conn_avgruns_',parcellation,'.mat')), 'mean_system_conn_mat')

%save the nodal participation coefficient and clustering coefficient
outfile=dataset(char(unique_subjlist.id0), char(unique_subjlist.id1), part_coef_avg_all)
export(outfile,'File',strcat(outdir,'/n74_part_coef_avg_nodewise_avgruns', parcellation,'.csv'),'Delimiter',',')

outfile=dataset(char(unique_subjlist.id0), char(unique_subjlist.id1), avgclustco_all)
export(outfile,'File',strcat(outdir,'/n74_clust_co_avg_nodewise_avgruns', parcellation,'.csv'),'Delimiter',',')

    end
end

%% Make an average connectivity matrix for each pipeline and parcellation
% Only include the weighted connectivity matrix for each subject
[unique_subs, index]=unique(subjlist.id0)
unique_subjlist= subjlist(index,:)
parcellations={'schaefer200', 'schaefer400'}
pipelines={'nogsr_spkreg_fd0.5dvars1.75_drpvls', 'gsr_spkreg_fd0.5dvars1.75_drpvls','gsr_censor_5contig_fd1.25dvars2_drpvls','nogsr_spkreg_fd1.25dvars2_drpvls','gsr_censor_5contig_fd0.5dvars1.75_drpvls'}
for p=1:length(parcellations)
    %for pl=1:length(pipelines)
        parcellation=parcellations{p}
        %pipeline=pipelines{pl}
        dim=str2double(parcellation(end-2:end))%what is the dimensionality of the parcellation
        stackedMatrix = zeros(dim, dim, height(unique_subjlist));
        datadir=strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,'/',parcellation,'zNetworks_avg')
        for n=1:height(unique_subjlist)
        sub=char(unique_subjlist.id0(n)) %look at this
        run=char(unique_subjlist.id1(n))
        file=fullfile(datadir,strcat(num2str(sub),'_',parcellation,'MNI_zavgnetwork.txt'))
        %file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
        subfcmat = load(file);
        for x=1:dim
            subfcmat(x,x)=0;
        end
        stackedMatrix(:,:, n)=subfcmat;
        end
meanMatrix = mean(stackedMatrix,3); %doc mean for more info.
%make the diagonal 0's
    for x=1:dim
        meanMatrix(x,x)=0;
    end
%export mean matrix
csvwrite(fullfile(outdir,strcat('averaged_FC_mat_n74_', parcellation,'.csv')), meanMatrix)
save(fullfile(outdir, strcat('averaged_FC_mat_n74_',parcellation,'.mat')), 'meanMatrix')
    %end
end