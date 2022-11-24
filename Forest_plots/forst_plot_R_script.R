## remove (almost) everything in the working environment.
## You will get no warning, so don't do this unless you are really sure.
rm(list = ls())

library("forestmodel")
library("survival")
library("dplyr")
library("survminer")

forest_plot <- function(work_dir, csv_name, feat, LN_status, event_name, censor_at, plot_width, plot_height){
  setwd(work_dir)
  tdata <- read.csv(file=csv_name, header=TRUE, sep=",", check.names=TRUE)
  
  ##drop cases with LN=+  ##Lymph.Node.status
  if (LN_status == 0)
    tdata <- tdata[!(tdata$Lymph.Node.status==1),]
  
  tdata$HGC1 = (tdata$HGC)
  tdata$Grade[tdata$Grade==1]<- "1"
  #names(tdata)  ##print the column names
  
  #tdata$HGC[tdata$HGC1<=median(tdata$HGC1)]<-"0"
  #tdata$HGC[tdata$HGC1>median(tdata$HGC1)]<-"1"
  tdata$HGC[tdata$HGC1<=0.92]<-"0"
  tdata$HGC[tdata$HGC1>0.92]<-"1"
  
  tdata$NPI2 = tdata$NPI
  #tdata$NPI[tdata$NPI2<=median(tdata$NPI2)]<-"0"
  #tdata$NPI[tdata$NPI2>median(tdata$NPI2)]<-"1"
  tdata$NPI[tdata$NPI2<=3.3]<-"0"
  tdata$NPI[tdata$NPI2>3.3]<-"1"

  tdata$Tumour_Size2 = tdata$Tumour_Size
  tdata$Tumour_Size[tdata$Tumour_Size2>median(tdata$Tumour_Size2)]<-"1"
  tdata$Tumour_Size[tdata$Tumour_Size2<=median(tdata$Tumour_Size2)]<-"0"
  
  if (event_name=="DMFS"){
    tdata$time = tdata$TTLR..month  ###@@@@@@@@ for DMFS: TTLR..month, for BCSS: Breast.cancer.specific.survival..month  @@@@@@@@@@@@@@@@@
    tdata$event = tdata$Distant.Metastasis ###@@@@@@@ for DMFS: Distant.Metastasis, for BCSS: Survival.Status  @@@@@@@@@@@@@@@@
  }
  
  if (event_name=="BCSS"){
    tdata$time = tdata$Breast.cancer.specific.survival..month  ###@@@@@@@@ for DMFS: TTLR..month, for BCSS: Breast.cancer.specific.survival..month  @@@@@@@@@@@@@@@@@
    tdata$event = tdata$Survival.Status ###@@@@@@@ for DMFS: Distant.Metastasis, for BCSS: Survival.Status  @@@@@@@@@@@@@@@@
  }

  ##censor the data  
  tdata$event[tdata$time > censor_at] <- 0  
  tdata$time[tdata$time > censor_at] <- censor_at
  
  ##convert events to 0, 1
  tdata$event[tdata$event > 1] <- 0
  
  if (feat == "HGC"){
    tdata <- tdata %>%
      transmute(time, event, Age = Age,
                'Tumour size' = factor(Tumour_Size, labels = c("<=1.4 cm",">1.4 cm")),
                LVI = factor(LVI, labels = c("No","Yes")),
                #'Menopause' = factor(Menopausal_status, labels = c("No","Yes")),
                'Multifocality' = factor(Multifocality, labels = c("No","Yes")),
                'Associated DCIS' = factor(Associated_DCIS, labels = c("No","Yes")),
                'Associated LCIS' = factor(Associated_LCIS, labels = c("No","Yes")),
                HGC = factor(HGC, labels = c("Low","High")),
      )
  }
  
  if (feat == "NPI"){
    tdata <- tdata %>%
      transmute(time, event, Age = Age,
                'Tumour size' = factor(Tumour_Size, labels = c("<=1.4 cm",">1.4 cm")),
                LVI = factor(LVI, labels = c("No","Yes")),
                #'Menopause' = factor(Menopausal_status, labels = c("No","Yes")),
                'Multifocality' = factor(Multifocality, labels = c("No","Yes")),
                'Associated DCIS' = factor(Associated_DCIS, labels = c("No","Yes")),
                'Associated LCIS' = factor(Associated_LCIS, labels = c("No","Yes")),
                NPI = factor(NPI, labels = c("Low","High")),
      )
  }
  
  if (feat == "Grade"){
    tdata <- tdata %>%
      transmute(time, event, Age = Age,
                'Tumour size' = factor(Tumour_Size, labels = c("<=1.4 cm",">1.4 cm")),
                LVI = factor(LVI, labels = c("No","Yes")),
                #'Menopause' = factor(Menopausal_status, labels = c("No","Yes")),
                'Multifocality' = factor(Multifocality, labels = c("No","Yes")),
                'Associated DCIS' = factor(Associated_DCIS, labels = c("No","Yes")),
                'Associated LCIS' = factor(Associated_LCIS, labels = c("No","Yes")),
                Grade = factor(Grade, labels = c("1, 2","3")),
      )
  }
  
  panels <- list(
    list(width = 0.01),
    list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
    list(width = 0.1, display = ~level, heading = "Strata"),
    list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
    
    list(width = 0.03, item = "vline", hjust = 0.5),
    list(width = 0.55, item = "forest", hjust = 0.5, heading = "", linetype = "dashed", line_x = 0),
    #list(width = 0.03, item = "vline", hjust = 0.5),
    list(width = 0.12, hjust=1, display = ~ ifelse(reference, "Reference", sprintf(
      "%0.2f (%0.2f, %0.2f)",
      trans(estimate), trans(conf.low), trans(conf.high)
    )), display_na = NA, heading = "HR (95% CI)"),
    list(width = 0.03, item = "vline", hjust = 0.5),
    list(
      width = 0.05,
      display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
      display_na = NA, hjust = 1, heading = "p"
    ),
    list(width = 0.03)
  )

  #mdl = coxph(Surv(time, event) ~ ., tdata)
  #mdl

  print(forest_model(coxph(Surv(time, event) ~ ., tdata), panels))
  ggsave("forest_plot.pdf", width=plot_width, height=plot_height)
}

##### usage of the above function for generating forest plots ######@@@@@@@@@@@@@@@
##the function forest_plot is expecting the event value=1 for an event and anything else will be set to 0 i.e. no event

##params: 
work_dir = 'D:/warwick/tasks/Nottingham_Exemplar/ML/cellular_analysis/R' ##path were the code, files are
csv_name = "HGC_UHCW_g1g2g3_tumor_per2_patchnew.csv" ##name of the csv file which contains the features
feature = "Grade" ##column name of the feature (<"HGC"|"Grade"|"NPI">)
ln_status = 0 ##lymph node status. 0 for LN-, 1 for LN 0-3
event_name = "BCSS" ##<"DMFS"|"BCSS">
censor_at = 120 ##months at which to censor the data (time, event)
plot_width = 9.15 ##width of the plot when saving
plot_height = 3.75 ##width of the plot when saving
forest_plot(work_dir, csv_name, feature, ln_status, event_name, censor_at, plot_width, plot_height)
