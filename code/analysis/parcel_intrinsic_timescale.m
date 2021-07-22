%% Setup
addpath(genpath('~/Documents/tools/MATLAB/intrisic-timescales/'))
pipeline='gsr_spkreg_fd0.5dvars1.75_drpvls'
timeseries_dir=strcat("/cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpEngine_", pipeline,"/")

% Set parameters
num_nodes = 400;    % number of time series 
outdir = strcat('~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/', pipeline,"/");    % set directory for saving out images here
mkdir(outdir)
lags = -6:6;    % range of TR shifts; should be sufficient to allow for all autocovariance functions (ACFs) to decay below .5

tr = 2;
motion_thresh = .5;    % important: must match motion criteria used during preproc
min_block_durn = (max(lags)+1)*tr;   % min. block duration (in seconds)

%read in subject list
listdir='/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
subjlist=readtable(fullfile(listdir,'n92_subjlist_for_timescales.csv'),'Delimiter',',','ReadVariableNames', 1)
uniquesubs=unique(subjlist.ID)
subjects=uniquesubs;

%% Loop over subjects

% initialize group matrices
grp_acfs = single(nan(num_nodes,numel(lags),numel(subjects))); % keep all ACFs for all subjects for stats

for s = 1:numel(uniquesubs)
    tic
    sub=char(uniquesubs(s));
    disp(['Processing ' sub]);
    
    % initialize subject matrices
    subj_lags = single(nan(num_nodes)); % peak lags
    subj_ZL = subj_lags;   % zero-lag correlation
    subj_peak = subj_lags; % peak correlation (correlation at optimal lag)
    
   %load timeseries matrix across runs--SHOULD PUT A 'CENSORED' VOLUME BETWEEN
   %FIRST RUN AND SECOND RUN TO DIVIDE UP CONTINUGUOUS-NESS?
    runs = subjlist(string(subjlist.ID)==sub, 'run') %get runs
    clear subts
    clear temp_mask mask
    subts=NaN(1,400);
    mask=[]
    for j=1:height(runs)
        run=runs.run(j)
        tsfile=strcat(timeseries_dir,sub,"/",run,"/fcon/schaefer400/",sub,"_",run,"_schaefer400_uncensored_ts.1D")
        tsfile2=strcat(timeseries_dir,sub,"/",run,"/fcon/schaefer400x7/",sub,"_",run,"_schaefer400x7_uncensored_ts.1D");
        if (exist(tsfile)==2)
             x=load(tsfile); % read in time series matrix
        elseif (exist(tsfile2)==2)
            x=load(tsfile2);
        end
            subts=vertcat(subts,x);
            temp_mask=load(strcat(timeseries_dir,sub,"/",run,"/",sub,"_", run,"-nFlags.1D")); %load nflags vector
            mask=vertcat(mask,temp_mask)
            %add a row of 1s between timeseries runs and at end of
            %temporal motion mask, to make sure an ACF block doesn't cross
            %runs.
            subts=vertcat(subts,repmat(NaN,1,400))
            mask=vertcat(mask,1)
        end
        size(subts)
        size(mask)
        %take off the last "masked frame" from each 
        subts(end,:)=[];
        mask(end,:)=[];
    subts(1,:) = [];%take out the first row with NaNs from timeseries matrix
    BOLD=subts;
    
    %resample to 135 TRs (the minimum amount of data to be included)
%     BOLD=BOLD(1:135,:);
%     mask=mask(1:135,:);

    %code to convert temporal mask/motion time series from 1s and 2s as excluded frames to 1s as included
    %frames, 0s as excluded frames
    [Aval, ~, indAval] = unique(mask);
    Avalnew = [1; 0; 0]; 
    Anew = Avalnew(indAval);
    mask = reshape(Anew, size(mask));
    format = mask ==1 %make format a logical vector
    
    % ignore pre-steady-state frames
    %format(1:2) = false; % ignore first X frames, turned this off.
    
    FORMAT = create_blocks(format,min_block_durn,tr);
    
    %% Construct ACF
    ACFs = single(zeros(num_nodes,numel(lags)));
    nblocks = numel(FORMAT);
    nframes = 0;
    
    % De-mean time series
    run_mean = nanmean(BOLD(format,:),1);
    BOLD = bsxfun(@minus,BOLD,run_mean);
    
    % Loop over blocks of contiguous frames
    for j = 1:numel(FORMAT)
        nframes = nframes + numel(FORMAT{j});
        FHCR = false(1,numel(format));
        FHCR(FORMAT{j}) = true;
        %for i = 1:sum(good)
        for i = 1:num_nodes
            %ACFs(i,:) = ACFs(i,:) + squeeze(lagged_cov(BOLD(FHCR,i),BOLD(FHCR,i),max(lags)))';
            %this line doesn't work transposed....
            ACFs(i,:) = ACFs(i,:) + squeeze(lagged_cov(BOLD(FHCR,i),BOLD(FHCR,i),max(lags)));
        end
    end
    
    % Normalize ACFs based on entire run
    for k = 1:numel(lags)
        ACFs(:,k) = ACFs(:,k)/(nframes - abs(lags(k))*nblocks);
    end
    ACFs = bsxfun(@rdivide,ACFs,ACFs(:,lags==0));
        
    % Store
    grp_acfs(:,:,s) = ACFs;

    toc
end

%Estimate FWHM lag 

acf_mean = nanmean(grp_acfs,3);
acf_mean = bsxfun(@rdivide,acf_mean,acf_mean(:,lags==0));

% fit
hwhm = acf_hwhm(acf_mean',tr); % group-wise

% subject-wise
hwhms = zeros(num_nodes,numel(subjects));
for s = 1:numel(subjects)
    tic
    hwhms(:,s) = acf_hwhm(grp_acfs(:,:,s)',tr);
    toc
end
%save out hwhm (groupwise at parcel level) and hwhms (at the subject level
save(fullfile(outdir, strcat('n92_schaefer400_timescales_censor_block_parcel_group.mat')), 'uniquesubs','hwhm', 'hwhms')

outfile=dataset(char(uniquesubs), hwhms')
export(outfile,'File',strcat(outdir,'/n92_schaefer400_timescales_censor_block_parcel_subjectwise.csv'),'Delimiter',',')

