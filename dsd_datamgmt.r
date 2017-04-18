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
  select(id, matches('^dom.[^tot]'), contains('rass')) %>%
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

rm(list = Cs(dsd.data, dsd.sumvars, had.feat))
