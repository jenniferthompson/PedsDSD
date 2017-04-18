####################################################################################################
## Delirium superimposed on dementia data management
## Using database sent from Alessandro, October 2014
## Creates one one-record-per-patient data frame and one longitudinal data frame from DOM and
##  baseline data on Italian geriatric patients assessed for delirium
####################################################################################################

library(Hmisc)
library(tidyr)
library(dplyr)

dsd.data <- read.csv('dsd_13oct2014.csv')
names(dsd.data) <- gsub('_', '.', tolower(names(dsd.data)), fixed = TRUE)

## Add ".g" to rass columns to make consistent with dom columns for reshaping later ##
names(dsd.data) <- gsub('rass\\.', 'rass\\.g\\.', names(dsd.data))

dsd.data <- upData(dsd.data,
  ## Discharge assessments of 999 indicate that patient died; change to NA
  tinetti.dc = ifelse(tinetti.dc == 999, NA, tinetti.dc),
  bi.dc = ifelse(bi.dc == 999, NA, bi.dc),

  ## Create factor versions of dementia type, discharge setting
  type.dem = factor(type.demenza,
                    levels = 0:3,
                    labels = c("Vascular", "Alzheimer's", "Lewy Body", "Other")),
  dc.loc = factor(ifelse(is.na(tinetti.dc), NA, discharge.setting),
                  levels = 0:5,
                  labels = c("Home", "Nursing home", "Other rehab", "Acute care", "Hospice", "Other")),
  
  ## Drop numeric versions of factor variables and variables which get recalculated from raw data ##
  drop = c('type.demenza', 'discharge.setting', 'dom.max', 'del.dur'),
  
  rename = c('study.id' = 'id',
             'bi.disch' = 'bi.dc',
             'tinetti.disch' = 'tinetti.dc',
             'delirium.duration' = 'del.dur',
             'episodes.hypoactive.del' = 'ep.hypo',
             'episodes.hyper.del' = 'ep.hyper'))

## -- Split data into one-record-per-person, longitudinal data sets -- ##
## One record per patient: baseline, discharge values; hypo/hyperactive delirium assessments ##
dsd.oneobs <- dsd.data %>% select(id:ep.hyper, type.dem, dc.loc)

## Transpose DOM assessment data to multiple records per patient, one column per question ##
dsd.long <- dsd.data %>%
  select(id, contains('^dom.[^tot]'), contains('rass')) %>%
    ## select DOM, RASS questions (leave out DOM totals)
  gather(dom.asmt, dom.value, -id) %>%                            ## transpose into columns
  separate(dom.asmt, into = c('dom.q', 'day'), sep = '.g\\.') %>% ## cols for DOM ?, study day
  spread(dom.q, dom.value) %>%                                    ## separate columns for each ?
  mutate(day = as.numeric(day)) %>%                               ## make day column numeric
  arrange(id, day)                                                ## sort by patient, day

## Remove any records where all DOM questions are missing
dsd.long <- dsd.long[rowSums(!is.na(dsd.long[,grep('^dom\\.[0-9]+', names(dsd.long))])) > 0,]

## Recalculate DOM total score for each day
dsd.long$dom.tot <- rowSums(dsd.long[,grep('^dom\\.[0-9]+', names(dsd.long))], na.rm = TRUE)

## Function to determine whether patient had each feature (score not missing and >0) ##
had.feat <- function(dom.score){
  return(ifelse(is.na(dom.score), NA, as.numeric(dom.score > 0)))
}

dsd.long <- dsd.long %>% 
  mutate(sustained.attn = had.feat(dom.1),
         shift.attn = had.feat(dom.2),
         attention = ifelse(sustained.attn == 1 | shift.attn == 1, 1, 0),
         orientation = had.feat(dom.3),
         consciousness = had.feat(dom.4),
         apathy = had.feat(dom.5),
         hypokinesia = had.feat(dom.6),
         incoherence = had.feat(dom.7),
         fluctuation = had.feat(dom.8),
         restlessness = had.feat(dom.9),
         delusions = had.feat(dom.10),
         hallucination = had.feat(dom.11),
         anxiety = had.feat(dom.12),
         del.type = factor(ifelse(is.na(rass), NA, ifelse(rass < 0, 1, ifelse(rass > 0, 2, 3))),
                           levels = 1:3, labels = c('Hypoactive', 'Hyperactive', 'RASS = 0')))

## Create variable for inconsolability - no equivalent in adults ##
dsd.long$inconsolability <- NA

## -- Calculate number of days with delirium assessments, max DOM score per patient -- ##
dsd.sumvars <- dsd.long %>%
  group_by(id) %>%
  summarise(dom.asmts = n(),
            dom.max = max(dom.tot),
            dom.asmts.hypo = sum(del.type == 'Hypoactive', na.rm = TRUE),
            dom.asmts.hyper = sum(del.type == 'Hyperactive', na.rm = TRUE),
            dom.asmts.rass0 = sum(del.type == 'RASS = 0', na.rm = TRUE))

## Merge with oneobs data
dsd.oneobs <- merge(dsd.oneobs, dsd.sumvars, by = 'id')

## -- Add variable labels -- ##
label(dsd.long$id) <- 'Patient ID'
label(dsd.long$day) <- 'Study day'
label(dsd.long$dom.1) <- 'DOM 1: Sustained attention'
label(dsd.long$dom.2) <- 'DOM 2: Shifting attention'
label(dsd.long$dom.3) <- 'DOM 3: Orientation'
label(dsd.long$dom.4) <- 'DOM 4: Consciousness'
label(dsd.long$dom.5) <- 'DOM 5: Apathy'
label(dsd.long$dom.6) <- 'DOM 6: Hypokinesia/psychomotor retardation'
label(dsd.long$dom.7) <- 'DOM 7: Incoherence'
label(dsd.long$dom.8) <- 'DOM 8: Fluctuation in functioning'
label(dsd.long$dom.9) <- 'DOM 9: Restlessness'
label(dsd.long$dom.10) <- 'DOM 10: Delusions'
label(dsd.long$dom.11) <- 'DOM 11: Hallucinations'
label(dsd.long$dom.12) <- 'DOM 12: Anxiety/fear'
label(dsd.long$rass) <- 'RASS'
label(dsd.long$dom.tot) <- 'Total DOM score'
label(dsd.long$sustained.attn) <- 'Sustained attention'
label(dsd.long$shift.attn) <- 'Shfiting attention'
label(dsd.long$attention) <- 'Attention'
label(dsd.long$orientation) <- 'Orientation'
label(dsd.long$consciousness) <- 'Consciousness'
label(dsd.long$apathy) <- 'Apathy'
label(dsd.long$hypokinesia) <- 'Hypokinesia/lethargy'
label(dsd.long$incoherence) <- 'Incoherence'
label(dsd.long$fluctuation) <- 'Fluctuations in functioning'
label(dsd.long$restlessness) <- 'Restlessness'
label(dsd.long$delusions) <- 'Delusions'
label(dsd.long$hallucination) <- 'Hallucinations'
label(dsd.long$anxiety) <- 'Anxiety'
label(dsd.long$inconsolability) <- 'Inconsolability'
label(dsd.long$del.type) <- 'Delirium subtype'

label(dsd.oneobs$id) <- 'Patient ID'
label(dsd.oneobs$age) <- 'Age'
label(dsd.oneobs$los) <- 'Length of stay (days)'
label(dsd.oneobs$iqcode) <- 'IQCODE'
label(dsd.oneobs$cdr) <- 'CDR'
label(dsd.oneobs$bi.pre) <- 'Barthel Index, pre-admission'
label(dsd.oneobs$bi.adm) <- 'Barthel Index, admission'
label(dsd.oneobs$bi.dc) <- 'Barthel Index, discharge'
label(dsd.oneobs$iadl.tot) <- 'IADL, admission'
label(dsd.oneobs$tinetti.adm) <- 'Tinetti score, admission'
label(dsd.oneobs$tinetti.dc) <- 'Tinetti score, discharge'
label(dsd.oneobs$ep.hypo) <- 'Episodes of hypoactive delirium'
label(dsd.oneobs$ep.hyper) <- 'Episodes of hyperactive delirium'
label(dsd.oneobs$type.dem) <- 'Type of dementia'
label(dsd.oneobs$dc.loc) <- 'Discharge setting'
label(dsd.oneobs$dom.asmts) <- 'Number of DOM assessments'
label(dsd.oneobs$dom.max) <- 'Maximum DOM score during study'
label(dsd.oneobs$dom.asmts.hypo) <- 'Number of hypoactive DOM assessments'
label(dsd.oneobs$dom.asmts.hyper) <- 'Number of hyperactive DOM assessments'
label(dsd.oneobs$dom.asmts.rass0) <- 'Number of DOM assessments with RASS = 0'

rm(list = Cs(dsd.data, dsd.sumvars, had.feat))
