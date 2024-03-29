
listdir='/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
subjlist=readtable(fullfile(listdir,'n74_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_80119.csv'),'Delimiter',',','ReadVariableNames', 1)

myqavars=readtable(strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/nogsr_spkreg_fd0.5dvars1.75_drpvls/n74_fixed_within_between_Yeo7_avgruns_schaefer400_nogsr_spkreg_fd0.5dvars1.75_drpvls_withmodulpartcoef_with_QA.csv'))
dim=str2double(parcellation(end-2:end))%what is the dimensionality of the parcellation
[unique_qvars, index]=unique(myqavars.ID)
unique_qvars= myqavars(index,:)

%% For each parcellation and each pipeline
parcellations={'schaefer200', 'schaefer400'}
pipelines={'nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls','gsr_censor_5contig_fd1.25dvars2_drpvls','nogsr_spkreg_fd1.25dvars2_drpvls','gsr_censor_5contig_fd0.5dvars1.75_drpvls'}
for p=1:length(parcellations)
    for pl=1:length(pipelines)
        parcellation=parcellations{p}
        pipeline=pipelines{pl}
        dim=str2double(parcellation(end-2:end))%what is the dimensionality of the parcellation

       datadir=strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,'/',parcellation,'zNetworks_avg')


for n=1:height(unique_subjlist)
    sub=char(unique_subjlist.id0(n)); %look at this
    run=char(unique_subjlist.id1(n));
    file=fullfile(datadir,strcat(num2str(sub),'_',parcellation,'MNI_zavgnetwork.txt'));
    %file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
    subfcmat = load(file);
    for x=1:dim
        subfcmat(x,x)=0;
    end

%average weight
avgweight(n,1)=mean(subfcmat(subfcmat~=0));
end
mean(avgweight)
corrcoef(unique_qvars.fd_mean_avg, avgweight)
    end
end

%% all runs straight from fcon
parcellations={'schaefer200', 'schaefer400'}
pipelines={'nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls','gsr_censor_5contig_fd1.25dvars2_drpvls','nogsr_spkreg_fd1.25dvars2_drpvls','gsr_censor_5contig_fd0.5dvars1.75_drpvls'}
for p=1:length(parcellations)
    for pl=1:length(pipelines)
        parcellation=parcellations{p}
        pipeline=pipelines{pl}
        dim=str2double(parcellation(end-2:end))%what is the dimensionality of the parcellation
datadir=strcat('~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_',pipeline);

for n=1:height(subjlist)
    sub=char(subjlist.id0(n)); %look at this
    run=char(subjlist.id1(n));
    file=fullfile(datadir,strcat(sub,'/',run,'/fcon/',parcellation,'/',sub,'_',run,'_',parcellation,'_network.txt'));
        subfcmat=load(file);
        %make into adjacency matrix and save out
        %size_vec=tril(ones(400,400),-1);
        size_vec=tril(ones(dim,dim),-1);
        adj_mat=size_vec;
        adj_mat(adj_mat==1)=subfcmat;
        subfcmat=adj_mat+adj_mat';
        for x=1:dim
            subfcmat(x,x)=0;
        end
        %make into adjacency matrix and save out
        %size_vec=tril(ones(400,400),-1);
        avgweight(n,1)=mean(subfcmat(subfcmat~=0));
end
fprintf('Average weight is %s \n', mean(avgweight));
fprintf('Correlation coefficient is %s \n', corrcoef(myqavars.fd_mean, avgweight))
end
end
