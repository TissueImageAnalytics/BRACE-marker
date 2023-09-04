## remove (almost) everything in the working environment.
## You will get no warning, so don't do this unless you are really sure.
rm(list = ls())

library("forestmodel")
library("survival")
library("dplyr")
library("survminer")

forest_plot <- function(work_dir, csv_name, feat, LN_status, event_name, censor_at, plot_width, plot_height, cohort_name){
  setwd(work_dir)
  tdata <- read.csv(file=csv_name, header=TRUE, sep=",", check.names=TRUE)

  ##drop cases with LN=+  ##Lymph.Node.status
  if (LN_status == 0)
    tdata <- tdata[!(tdata$Lymph.Node.status==1),]
  
  tdata$BRACE1 = (tdata$BRACE)
  tdata$Grade[tdata$Grade==1]<- "1"
  
  tdata$BRACE[tdata$BRACE1<=1.28]<-"0"
  tdata$BRACE[tdata$BRACE1>1.28]<-"1"
  
  tdata$Magee2_no_npi = tdata$Magee_new2_no_npi
  if (cohort_name=="NOTT"){
    tdata$Magee2_nonpi[tdata$Magee2_no_npi<=11]<-"0"
    tdata$Magee2_nonpi[tdata$Magee2_no_npi>11]<-"1"}
  
  tdata$Age2 = tdata$Age
  tdata$Age[tdata$Age2>=50]<-"1"
  tdata$Age[tdata$Age2<50]<-"0"
  
  tdata$Tumour_Size2 = tdata$Tumour_Size
  tdata$Tumour_Size[tdata$Tumour_Size2>median(tdata$Tumour_Size2)]<-"1"
  tdata$Tumour_Size[tdata$Tumour_Size2<=median(tdata$Tumour_Size2)]<-"0"
  
  if (event_name=="DMFS"){
    tdata$time = tdata$TTDM..month  ###@@@@@@@@ for DMFS: TTDM..month, for BCSS: Breast.cancer.specific.survival..month  @@@@@@@@@@@@@@@@@
    tdata$event = tdata$Distant.Metastasis ###@@@@@@@ for DMFS: Distant.Metastasis, for BCSS: Survival.Status  @@@@@@@@@@@@@@@@
  }
  
  if (event_name=="BCSS"){
    tdata$time = tdata$Breast.cancer.specific.survival..month  ###@@@@@@@@ for DMFS: TTDM..month, for BCSS: Breast.cancer.specific.survival..month  @@@@@@@@@@@@@@@@@
    tdata$event = tdata$Survival.Status ###@@@@@@@ for DMFS: Distant.Metastasis, for BCSS: Survival.Status  @@@@@@@@@@@@@@@@
  }

  ##censor the data  
  tdata$event[tdata$time > censor_at] <- 0  
  tdata$time[tdata$time > censor_at] <- censor_at
  
  ##convert events to 0, 1
  tdata$event[tdata$event > 1] <- 0
  
  if (feat == "BRACE"){
    if(cohort_name == 'NOTT'){
      if (LN_status == 0){
        tdata <- tdata %>%
        transmute(time, event, 
                'Age at diagnosis' = factor(Age, labels = c("<50",">=50")),
                'Multifocality' = factor(Multifocality, labels = c("No","Yes")),
                'DCIS' = factor(Associated_DCIS, labels = c("No","Yes")),
                'LCIS' = factor(Associated_LCIS, labels = c("No","Yes")),
                LVI = factor(LVI, labels = c("No","Yes")),
                Grade = factor(Grade, labels = c("1, 2","3")),
                Magee2 = factor(Magee2_nonpi, labels = c("Low","High")),
                BRACE = factor(BRACE, labels = c("Low","High")),)
      }
      if (LN_status == 1){
        tdata <- tdata %>%
          transmute(time, event, 
                'Age at diagnosis' = factor(Age, labels = c("<50",">=50")),
                'Multifocality' = factor(Multifocality, labels = c("No","Yes")),
                'DCIS' = factor(Associated_DCIS, labels = c("No","Yes")),
                'LCIS' = factor(Associated_LCIS, labels = c("No","Yes")),
                'Lymph node status' = factor(Lymph.Node.status, labels = c("No","Yes")),
                LVI = factor(LVI, labels = c("No","Yes")),
                Grade = factor(Grade, labels = c("1, 2","3")),
                Magee2 = factor(Magee2_nonpi, labels = c("Low","High")),
                BRACE = factor(BRACE, labels = c("Low","High")),)
      }
    }
    if(cohort_name == 'UHCW'){
      if (LN_status == 0){
        tdata <- tdata %>%
          transmute(time, event, 
                'Age at diagnosis' = factor(Age, labels = c("<50",">=50")),
                'Multifocality' = factor(Multifocality, labels = c("No","Yes")),
                'DCIS' = factor(Associated_DCIS, labels = c("No","Yes")),
                'LCIS' = factor(Associated_LCIS, labels = c("No","Yes")),
                'Tumour size' = factor(Tumour_Size, labels = c("<2 cm",">=2 cm")),
                Grade = factor(Grade, labels = c("1, 2","3")),
                BRACE = factor(BRACE, labels = c("Low","High")),)
      }
      if (LN_status == 1){
        tdata <- tdata %>%
          transmute(time, event, 
                'Age at diagnosis' = factor(Age, labels = c("<50",">=50")),
                'Multifocality' = factor(Multifocality, labels = c("No","Yes")),
                'DCIS' = factor(Associated_DCIS, labels = c("No","Yes")),
                'LCIS' = factor(Associated_LCIS, labels = c("No","Yes")),
                'Tumour size' = factor(Tumour_Size, labels = c("<2 cm",">=2 cm")),
                'Lymph node status' = factor(Lymph.Node.status, labels = c("No","Yes")),
                Grade = factor(Grade, labels = c("1, 2","3")),
                BRACE = factor(BRACE, labels = c("Low","High")),)
      }
    }
  }
  
  panels <- list(
    list(width = 0.01),
    list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
    list(width = 0.1, display = ~level, heading = "Strata"),
    list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
    
    list(width = 0.03, item = "vline", hjust = 0.5),
    list(width = 0.55, item = "forest", hjust = 0.5, heading = "", linetype = "dashed", line_x = 0),
    list(width = 0.12, hjust=1, display = ~ ifelse(reference, "Reference", sprintf(
      "%0.2f (%0.2f, %0.2f)",
      trans(estimate), trans(conf.low), trans(conf.high)
    )), display_na = NA, heading = "HR (95% CI)"),
    list(width = 0.03, item = "vline", hjust = 0.5),
    list(
      width = 0.05,
      display = ~ ifelse(reference, "", format.pval(p.value, digits = 1)),
      display_na = NA, hjust = 1, heading = "P"
    ),
    list(width = 0.03)
  )

  print(forest_model(coxph(Surv(time, event) ~ ., tdata), panels))
  ggsave(paste(cohort_name,  feat, event_name, ln_status, ".pdf"), width=plot_width, height=plot_height)
}

##### usage of the above function for generating forest plots ######@@@@@@@@@@@@@@@
##the function forest_plot is expecting the event value=1 for an event and anything else will be set to 0 i.e. no event

##params: 
work_dir = '/Rcode/' ##path where the code and other files are
cohorts <- list("NOTT","UHCW") ## <"NOTT"|"UHCW">

feature_list <- list("BRACE")
ln_status_list <- list(0, 1) ##lymph node status. 0 for LN-, 1 for LN 0-3
event_name_list <- list("DMFS","BCSS") ##<"DMFS"|"BCSS">
censor_at = 120 ##months at which to censor the data (time, event)
plot_width = 12 ##width of the plot when saving
plot_height = 5 ##width of the plot when saving

for (c in cohorts){
  if (c == "NOTT"){
    csv_name = "NOTT_combi_valid_magee2_onco.csv"
  }
  if (c == "UHCW"){
    csv_name = 'UHCW_BRACE_test_onco_modified.csv'
  }
  
  for (feature in feature_list)
    for (event_name in event_name_list)
      for (ln_status in ln_status_list)
        forest_plot(work_dir, csv_name, feature, ln_status, event_name, censor_at, plot_width, plot_height, c)
}
