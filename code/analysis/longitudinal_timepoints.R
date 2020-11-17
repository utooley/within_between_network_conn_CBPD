library(dplyr)
library(reshape2)
library(ggplot2)
data_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/"
data_maintenance_day_dir="~/Box/Mackey_Lab/CBPD (825656)/Data/Data Maintenance Day/CURRENT_2019-09-Data/"

alldata=read.csv(paste0(data_dir,"CBPD_data_redcapexport_2020.1.31.csv"))
long_t2_time_gaps_only=read.csv(paste0(data_dir,"long_t2_data_gaps_1.31.20.csv"))

# Data Cleaning -----------------------------------------------------------

#make new id to match on, match the longitudinal time gaps from one file to the baseline data in the other
long_t2_time_gaps_only$idlong <- as.character(long_t2_time_gaps_only$record_id)
long_t2_time_gaps_only$record_id <- sub(pattern = "_[23]", "",long_t2_time_gaps_only$idlong)
long_t2_time_gaps_only$record_id

#reshape longitudinal data to wide and then merge with baseline data
long_t2_time_gaps_only <- reshape(long_t2_time_gaps_only, idvar="record_id", timevar="longitudinal_visit_num", direction="wide")

#filter out unneeded cols in redcap data
alldata <- alldata %>% select(., -c(general_information_complete:mri_information_complete))

all <- left_join(long_t2_time_gaps_only, alldata, by= "record_id")
dim(all)

#Look at distribution of baseline ages and age gaps
all$date_behav<- strptime(all$date_behav, format="%m/%d/%y")
all$date_behav.2<- strptime(all$date_behav.2, format="%m/%d/%y")
all$date_behav.3<- strptime(all$date_behav.3, format="%m/%d/%y")

all$gapbetween_T2 <- difftime(all$date_behav.2, all$date_behav) 
all$gapbetween_T2mnths <- as.numeric(difftime(all$date_behav.2, all$date_behav)/30) #assumption of 30-day month

all$gapbetween_T3 <- difftime(all$date_behav.3, all$date_behav.2) 
all$gapbetween_T3mnths <- as.numeric(difftime(all$date_behav.3, all$date_behav.2)/30) #assumption of 30-day month


# Gaps between Timepoints Viz ---------------------------------------------


#age at t1
hist(all$age_behav, main="Age at T1", col="lightblue")
summary(all$age_behav)
#age at longitudinal T2
hist(all$age_behav.2, main="Age at T2", col="blue")
summary(all$age_behav.2)
#age at t3
hist(all$age_behav.3, main="Age at T3", col="blue")
summary(all$age_behav.3)

#gap between T1 and T2
hist(all$gapbetween_T2mnths, main="T1-T2 (months)", col="darkblue")
summary(all$gapbetween_T2mnths)

#gap between T2 and T3
hist(all$gapbetween_T3mnths, main="T2-T3 (months)", col="darkblue")
summary(all$gapbetween_T3mnths)


# Spaghetti Plot ----------------------------------------------------------

small <- all %>% select(., record_id, age_behav, age_behav.2, age_behav.3) %>% melt(.,id_vars="record_id")
small$variable <- "age_behav"
small <- small %>%  arrange(value)
small <- transform(small, record_id=reorder(record_id, -value) ) 

small$record_id<- factor(small$record_id, levels=unique(as.character(small$record_id)))

p <- ggplot(data =small, aes(x = value,y=record_id))
p + geom_line()+geom_point()+xlab(label = "Age (years)")+ylab("")+ ggtitle("Longitudinal data Feb 2020")+ theme(axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks = element_blank())+ theme_classic()+scale_y_discrete(labels=c(rep("",44)))+  theme(plot.title = element_text(hjust = 0.5)) 


