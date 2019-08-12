%% Extract only edges between DMN-DAN and DMN-VAN and VIS-DAN to look at age and age x SES effects
listdir='/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'

%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n75_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_80119.csv'),'Delimiter',',','ReadVariableNames', 1)
%subjlist=readtable('/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n64_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_060419.csv', 'Delimiter',',')
%subjlist=subjlist(:,1:2);

%% For each parcellation and each pipeline
parcellations={'schaefer200', 'schaefer400'}
pipelines={'gsr_censor_5contig_fd0.5dvars1.75_drpvls', 'gsr_censor_5contig_fd1.25dvars2_drpvls','nogsr_spkreg_fd1.25dvars2_drpvls'}
for p=1:length(parcellations)
    for pl=1:length(pipelines)
        parcellation=parcellations{p}
        pipeline=pipelines{pl}
   
%start here for one pipeline
pipeline='gsr_censor_5contig_fd0.5dvars1.75_drpvls'
parcellation='schaefer400'
dim=str2double(parcellation(end-2:end))
%running with the cluster mounted locally
listdir='/Users/utooley/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
z_outdir=fullfile('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,strcat(parcellation,'zNetworks'))

datadir=strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline,'/',parcellation,'zNetworks')
yeo_nodes=dlmread(fullfile('~/Desktop/cluster/picsl/mackey_group/tools',parcellation,strcat(parcellation,'x7CommunityAffiliation.1D.txt')))

size_vec=triu(ones(400,400),1)
for n=1:height(subjlist)
    %% rewrite edges for each subject into a vector
    sub=char(subjlist.id0(n)) %look at this
    run=char(subjlist.id1(n))
    file=fullfile(datadir,strcat(num2str(sub),'_',run,'_',parcellation,'MNI_znetwork.txt'))
    %file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
    subfcmat = load(file);
    for x=1:dim
        subfcmat(x,x)=0;
    end
    vector_edges(n,:)=subfcmat(size_vec==1);
end
outdir=strcat('~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/',pipeline)
outfile=dataset(subjlist.id0,subjlist.id1, vector_edges)
export(outfile,'File',fullfile(outdir,'zedges_for_each_subj_081219.csv'),'Delimiter',',')


%read in pvalues from R

indexofcols=ismemberf(yeo_nodes,3);
indexofrows=ismemberf(yeo_nodes,4);

%get indices of edges that are DMN-DAN
%what are the significant nodes?
indexofnodes=linspace(1,359,359)
indexofnonsignodes=setdiff( indexofnodes, indexofcols)
%read in the betas for edge weights for age x SES interaction
edgebetas=csvread(fullfile(clustcodir,'edge_betas_agexses_int_scaled.csv'),1, 1 );
%edgebetas(:,1)=[]
%make a matrix with only the 26 nodes of significance and their edgess
%pull out only the 26 nodes that show an age x SES effect colnum-wise.
    for l=1:length(indexofcols);
    i=indexofcols(l);
    signodesonly(:,l)=edgebetas(:,i);
    end
%pull out only the 26 nodes that show an age x SES effect rowise.
    for j=1:length(indexofcols);
    i=indexofcols(j);
    signodesonly2(j,:)=signodesonly(i,:);
    end
    
%look at the avg value within those 26 node
avgbetainside=mean(signodesonly2(signodesonly2~=0))
%look at average value outside the 26 nodes
%make a matrix with only the 26 nodes of significance and their edges
a=edgebetas([indexofcols],:)
avgbetabetween=mean(a(a~=0))

    %make one large matrix with only DMN-DAN edges
    
    %make one large matrix with only DMN-VAN edges
    
    %make one large matrix with only VIS-DAN edges
    
    
    %write these out
    
    
    %do the regression in R for which are significant
    
    
    
    %read the pvals or betas for edges back in
    
    
    %if an edge is significant, mark both nodes with that significant
    %value.
    
    
    %assign those nodes to colors of the yeo annotation and write it out.
    