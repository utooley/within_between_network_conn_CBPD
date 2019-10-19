library(dplyr)
library(reshape2)
library(ggplot2)
data_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/"
data_maintenance_day_dir="~/Box/Mackey_Lab/CBPD (825656)/Data/Data Maintenance Day/CURRENT_2019-09-Data/"

alldata=read.csv(paste0(data_maintenance_day_dir,"CBPD_data_DMD_2019.10.17.csv"))
long_t2_time_gaps_only=read.csv(paste0(data_dir,"long_t2_data_gaps.csv"), )

#make new id to match on, match the longitudinal time gaps from one file to the baseline data in the other
long_t2_time_gaps_only$idlong <- as.character(long_t2_time_gaps_only$id)
long_t2_time_gaps_only$record_id <- sub(pattern = "_2", "",long_t2_time_gaps_only$idlong)
long_t2_time_gaps_only$record_id

all <- left_join(long_t2_time_gaps_only, alldata, by= "record_id")
dim(all)

#Look at distribution of baseline ages and age gaps
all$date_behav.x<- strptime(all$date_behav.x, format="%m/%d/%y")
all$date_behav.y<- strptime(all$date_behav.y, format="%m/%d/%y")

all$gapbetween <- difftime(all$date_behav.x, all$date_behav.y) 
all$gapbetweenmnths <- as.numeric(difftime(all$date_behav.x, all$date_behav.y)/30) #assumption of 30-day month

#age at t1
hist(all$age_behav.y, main="Age at T1", col="lightblue")
summary(all$age_behav.y)
#age at t2
hist(all$age_behav.x, main="Age at T2", col="blue")
summary(all$age_behav.x)

#gap between them
hist(all$gapbetweenmnths, main="T1-T2 (months)", col="darkblue")
summary(all$gapbetweenmnths)
all$date_behav <- all$date_behav.x
all$date_behav <- all$date_behav.y
small <- all %>% select(., record_id, age_behav.x, age_behav.y) %>% melt(.,id_vars="record_id")
small$variable <- "age_behav"
small$record_id2 <- small$record_id
p <- ggplot(data =small, aes(x = value,y=record_id2))
p + geom_line()+xlab(label = "Age (years)")+ylab("")+ ggtitle("Longitudinal data thus far")+ theme(axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks = element_blank())+ theme_classic()+scale_y_discrete(labels=c(rep("",44)))+  theme(plot.title = element_text(hjust = 0.5)) 


