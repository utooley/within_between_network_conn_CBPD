addpath(genpath('~/Documents/tools/MATLAB/intrisic_timescales'))
pipeline='gsr_spkreg_fd0.5dvars1.75_drpvls'
timeseries_dir=strcat("/cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpEngine_", pipeline,"/")

%read in subject list
listdir='/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/'
subjlist=readtable(fullfile(listdir,'n92_subjlist_for_timescales.csv'),'Delimiter',',','ReadVariableNames', 1)
uniquesubs=unique(subjlist.ID)

for n=1:height(uniquesubs)
    sub=char(uniquesubs(n)) %look at this
    runs = subjlist(string(subjlist.ID)==sub, 'run')
    %empty_ts=NaN(1,400)
    clear subts
    for j=1:height(runs)
        run=runs.run(j)
        tsfile=strcat(timeseries_dir,sub,"/",run,"/fcon/schaefer400/",sub,"_",run,"_schaefer400_uncensored_ts.1D")
        if (exist(tsfile)==2)
             x=load(tsfile);
             subts=horzcat(subts,x)
             %load nflags vector
             
        else
            tsfile=strcat(timeseries_dir,sub,"/",run,"/fcon/schaefer400x7/",sub,"_",run,"_schaefer400x7_uncensored_ts.1D")
            x=load(tsfile);
            subts=horzcat(subts,x)  
        end
        size(subts)
    end
    sub_parcel_ts[[subject]] <- na.omit(subts)
        
    %do intrinsic timescale
    


#     x <- rbind(x,as.matrix(read.table(file =tsfile, header = F)))   
#    else{
#     tsfile=paste0(timeseries_dir,subject,"/",run,"/fcon/schaefer400x7/",subject,"_",run,"_schaefer400x7_ts.1D")
#      x <- rbind(x,as.matrix(read.table(file =tsfile, header = F)))
#    }
#   }
# sub_parcel_ts[[subject]] <- na.omit(x)
# }
    
    
    
    rows = (string(myqavars.ID)==sub & string(myqavars.run)==run);
 vars = {'nVolCensored', 'size_t', 'totalnVolCensored','totalSizet'};
                weights=myqavars(rows, vars);
                
    for (subject in list$ID){
#   print(subject)
#   runs <- list[list$ID==subject,"run"]
#   x <- matrix(NA,nrow=1,ncol = 400)
#   for (run in runs){
#     tsfile=paste0(timeseries_dir,subject,"/",run,"/fcon/schaefer400/",subject,"_",run,"_schaefer400_ts.1D")
#    if(file.exists(tsfile)){
#     x <- rbind(x,as.matrix(read.table(file =tsfile, header = F))) 
#    }
#    else{
#     tsfile=paste0(timeseries_dir,subject,"/",run,"/fcon/schaefer400x7/",subject,"_",run,"_schaefer400x7_ts.1D")
#      x <- rbind(x,as.matrix(read.table(file =tsfile, header = F)))
#    }
#   }
# sub_parcel_ts[[subject]] <- na.omit(x)
# }