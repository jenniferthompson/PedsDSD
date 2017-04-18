####################################################################################################
## PsCAM data management for pediatric vs delirium superimposed on dementia project
## Using raw data from PsCAM REDCap database
## Creates one one-record-per-patient data frame and one longitudinal data frame from DOM and
##  baseline data on pediatric patients assessed for delirium
####################################################################################################

library(RCurl)
library(Hmisc)
library(tidyr)
library(dplyr)

curl_handle = getCurlHandle()
curlSetOpt(ssl.verifypeer = FALSE, connecttimeout = 0, timeout = 600, curl = curl_handle)

source('get_pscam_redcap.R')

if(Sys.info()['sysname'] == 'Darwin'){
  source('/Volumes/thomps23/ICUDelirium/PsCAM/pscam_redcap.r')
  mult.ids <- read.csv('/Volumes/thomps23/ICUDelirium/PsCAM/pscam_repeatedpts.csv')
  #   source('pscam_redcap.r')
} else{
  source('/home/thomps23/ICUDelirium/PsCAM/pscam_redcap.r')
  ## Get IDs of patients with multiple hospitalizations
  mult.ids <- read.csv('/home/thomps23/ICUDelirium/PsCAM/pscam_repeatedpts.csv')
}

names(pscam.data) <- tolower(gsub('_+', '.', names(pscam.data)))

## For pts hospitalized multiple times, replace study ID on later hospitalizations with initial ID ##
pscam.data <- merge(pscam.data, mult.ids, by = 'study.id', all = TRUE)
pscam.data$id <- ifelse(is.na(pscam.data$id), pscam.data$study.id, pscam.data$id)

####################################################################################################
## Demographic Data
####################################################################################################

## Take only records with demographic data (assessment 1) ##
pscam.oneobs <- subset(pscam.data, redcap.event.name.factor == 'Assessment 1')

## -- Calculate length of stay -- ##
make.date <- function(dt){ return(as.Date(as.character(dt), format = '%Y-%m-%d')) }

pscam.oneobs$enroll.date <- make.date(pscam.oneobs$enroll.day)
pscam.oneobs$hospdis.date <- make.date(pscam.oneobs$discharge)
pscam.oneobs$death.date <- make.date(pscam.oneobs$death.date)

pscam.oneobs$last.inhosp <- with(pscam.oneobs,{
  as.Date(ifelse(!is.na(hospdis.date), hospdis.date,
          ifelse(!is.na(death.date), death.date, NA)),
          origin = '1970-1-1') })

pscam.oneobs$los <- with(pscam.oneobs, as.numeric(difftime(last.inhosp, enroll.date, units = 'days')))

## -- Patient status at end of hospital stay -- ##
## If discharge date recorded, list that, even if patient withdrew
pscam.oneobs$dc.status <- with(pscam.oneobs, {
  factor(ifelse(!is.na(hospdis.date), 1,
         ifelse(!is.na(death.date), 2,
         ifelse(!is.na(withdraw.factor) & withdraw.factor == 'Yes', 3, 4))),
         levels = 1:4,
         labels = c('Discharged alive', 'Deceased', 'Withdrew; discharge status unknown', 'Unknown')) })

## -- PRISM -- ##
pscam.oneobs <- pscam.oneobs %>%
  mutate(
    ## Cardiovascular/neurologic vital signs ##
    prism.sbp = ifelse(is.na(systbp.infant.factor) & is.na(systbp.child.factor), NA,
                ifelse((!is.na(systbp.infant.factor) & systbp.infant.factor == '45-65') |
                         (!is.na(systbp.child.factor) & systbp.child.factor == '55-75'), 3,
                ifelse((!is.na(systbp.infant.factor) & systbp.infant.factor == '< 45') |
                         (!is.na(systbp.child.factor) & systbp.child.factor == '< 55'), 7, 0))),
    prism.temp = ifelse(is.na(temp.factor), NA, ifelse(temp.factor == 'Normal', 0, 3)),
    prism.mental = ifelse(is.na(coma.factor), NA, ifelse(coma.factor == 'No', 0, 5)),
    prism.hr = ifelse(is.na(hr.infant.factor) & is.na(hr.child.factor), NA,
               ifelse((!is.na(hr.infant.factor) & hr.infant.factor == '215-225') |
                        (!is.na(hr.child.factor) & hr.child.factor == '185 - 205'), 3,
               ifelse((!is.na(hr.infant.factor) & hr.infant.factor == '> 225') |
                        (!is.na(hr.child.factor) & hr.child.factor == '> 205'), 4, 0))),
    prism.pupils = ifelse(is.na(pupils.factor), NA,
                   ifelse(pupils.factor == 'One pupil fixed, pupil > 3 mm', 7,
                   ifelse(pupils.factor == 'Both pupils fixed, pupil > 3 mm', 11, 0))),
    
    ## Acid-base/blood gases ##
    prism.acid = ifelse(is.na(acidosis.factor), NA,
                 ifelse(acidosis.factor == 'pH 7.0 -7.28 or total CO2 5-16.9 mmHg', 2,
                 ifelse(acidosis.factor == 'pH < 7.0 or total CO2 < 5 mmHg', 6, 0))),
    prism.ph = ifelse(is.na(ph.factor), NA,
               ifelse(ph.factor == '7.48 - 7.55', 2,
               ifelse(ph.factor == '>7.55', 3, 0))),
    prism.pco2 = ifelse(is.na(pco2.factor), NA,
                 ifelse(pco2.factor == '50.0 - 75.0', 1,
                 ifelse(pco2.factor == '> 75.0', 3, 0))),
    prism.co2 = ifelse(is.na(co2.factor), NA, ifelse(co2.factor == '> 34.0', 4, 0)),
    prism.pao2 = ifelse(is.na(pao2.factor), NA,
                 ifelse(pao2.factor == '42.0 - 49.9', 3,
                 ifelse(pao2.factor == '< 42.0', 6, 0))),
    
    ## Chemistry ##
    prism.glucose = ifelse(is.na(glucose.sisi.factor), NA,
                    ifelse(glucose.sisi.factor == '> 200', 2, 0)),
    prism.k = ifelse(is.na(potassium.factor), NA, ifelse(potassium.factor == '>6.9', 3, 0)),
    prism.cr = ifelse(is.na(creatinine.factor), NA, ifelse(creatinine.factor == '>0.90', 2, 0)),
    prism.bun = ifelse(is.na(urea.factor), NA, ifelse(urea.factor == '> 14.9', 3, 0)),
    
    ## Hematology ##
    prism.wbc = ifelse(is.na(wbc.factor), NA, ifelse(wbc.factor == '< 3000', 4, 0)),
    prism.plt = ifelse(is.na(platelets.si.factor), NA,
                ifelse(platelets.si.factor == '100,000 - 200,000', 2,
                ifelse(platelets.si.factor == '50,000 - 99,999', 4,
                ifelse(platelets.si.factor == '< 50,000', 5, 0)))),
    prism.pt = ifelse(is.na(pt.ptt.factor), NA,
               ifelse(pt.ptt.factor == 'PT > 22.0 or PTT > 57.0', 3, 0)))

## Calculate total PRISM score ##
prism.vars <- c('sbp', 'temp', 'mental', 'hr', 'pupils', 'acid', 'ph', 'pco2', 'co2', 'pao2',
                'glucose', 'k', 'cr', 'bun', 'wbc', 'plt', 'pt')

pscam.oneobs$prism.score <- rowSums(pscam.oneobs[,paste0('prism.', prism.vars)])
pscam.oneobs <- pscam.oneobs %>%
  mutate(prism.cat = factor(ifelse(is.na(prism.score), NA,
                            ifelse(prism.score > 20, 3,
                            ifelse(prism.score > 10, 2, 1))),
                            levels = 1:3,
                            labels = c('Low risk (0-10)', 'Moderate risk (11-20)', 'High risk (>20)')))

write.csv(pscam.oneobs[sample(1:nrow(pscam.oneobs), size = 7),
                       c('id', paste0('prism.', prism.vars), 'prism.score')],
          file = 'pscam_checkprism.csv', row.names = FALSE)

pscam.oneobs <- pscam.oneobs[,c('id', 'age', 'sex.factor',
                                grep('^admit\\.dx\\.[0-9]+$', names(pscam.oneobs), value = TRUE),
                                'prism.score', 'los', 'dc.status')]


####################################################################################################
## Assessment Data
####################################################################################################

## Create indicators for whether delirium variables "should" be present:
## Old ref. rater form: rater must be filled out with name ##+ level of consciousness not stupor/coma
## New ref. rater form: reference rater marked as present  ##+ level of consciousness not stupor/coma
## UPDATE Assessment can still be done with LOC of stupor or coma

## Variable ref.rater filled out with a name; get list of all text entries indicating no assessment
ref.rater.no <- c('', 'none', 'n/a', 'evaluation not performed', 'no rater available ',
                  'no rater available', 'no reference rater', 'no reference rater ',
                  'rater not available', 'psych did not see', 'psych unavailable ',
                  'psych not present ', 'not available ', 'none ', 'form not available ',
                  'no reference rater', 'missing form', 'none present')

## Patient able to be assessed for delirium per original rater form
pscam.data$del.rater.1 <- with(pscam.data, {
  ifelse(!is.na(ref.rater) & ref.rater %nin% ref.rater.no, # &
           # !is.na(ref.conscious) & ref.conscious < 8,
         TRUE, FALSE) })

## Patient able to be assessed for delirium per updated rater form
pscam.data$del.rater.2 <- with(pscam.data, {
  ifelse(!is.na(ref.rater2) & ref.rater2.factor == 'Yes', # &
           # !is.na(ref.conscious2) & ref.conscious2 < 8,
         TRUE, FALSE) })

## -- Create abnormal(T)/normal(F)/missing(NA) versions of variables involved in symptoms -- ##

## Attention ##
## Original reference rater form
pscam.data$focus.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.focus) | ref.focus.factor == 'Unable to Assess', NA,
         ref.focus.factor == 'Yes,') })
pscam.data$sustain.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.sustain) | ref.sustain.factor == 'Unable to Assess', NA,
         ref.sustain.factor == 'Yes,') })
pscam.data$screen.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.screen) | ref.screen.factor == 'Unable to Assess', NA,
         ref.screen.factor == 'Yes,') })
pscam.data$shift.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.shift) | ref.shift.factor == 'Unable to Assess', NA,
         ref.shift.factor == 'Yes,') })

## New reference rater form
pscam.data$attention.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.attention2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.attention2.factor == 'Yes') })
pscam.data$sustained.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.susattn2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.susattn2.factor == 'Yes') })
pscam.data$stimuli.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.stim2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.stim2.factor == 'Yes') })
pscam.data$shift.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.shiftattn2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.shiftattn2.factor == 'Yes') })

vars.sustained.attn <- Cs(focus.1, sustain.1, screen.1, attention.2, sustained.2, stimuli.2)
vars.shift.attn <- Cs(shift.1, shift.2)
vars.attention <- c(vars.sustained.attn, vars.shift.attn)

## Orientation ##
## Old form
pscam.data$orient.person.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.orient.person) | ref.orient.person.factor == 'Unable to Assess', NA,
         ref.orient.person.factor == 'Yes,') })
pscam.data$orient.place.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.orient.place) | ref.orient.place.factor == 'Unable to Assess', NA,
         ref.orient.place.factor == 'Yes,') })
pscam.data$orient.time.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.orient.time) | ref.orient.time.factor == 'Unable to Assess', NA,
         ref.orient.time.factor == 'Yes,') })
pscam.data$orient.situation.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.orient.situation) | ref.orient.situation.factor == 'Unable to Assess',
         NA,
         ref.orient.situation.factor == 'Yes,') })

## New form
pscam.data$consistentpref.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.consispref2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.consispref2.factor == 'Yes') })
pscam.data$orient.person.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.orientation2.1) | ref.orientation2.0 == 1 |
           is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.orientation2.1 == 1) })
pscam.data$orient.place.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.orientation2.2) | ref.orientation2.0 == 1 |
           is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.orientation2.2 == 1) })

vars.orientation <- Cs(orient.person.1, orient.place.1, orient.time.1, orient.situation.1,
                       consistentpref.2, orient.person.2, orient.place.2)

## Consciousness ##
pscam.data$conscious.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.conscious), NA, ref.conscious != 4) })
pscam.data$conscious.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.conscious2), NA, ref.conscious2 != 4) })

vars.consciousness <- Cs(conscious.1, conscious.2)

## Apathy ##
## Old form
pscam.data$interest.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.interest) | ref.interest.factor == 'Unable to Assess', NA,
         ref.interest.factor == 'Yes,') })
pscam.data$inappropriate.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.inappropriate) | ref.inappropriate.factor == 'Unable to Assess', NA,
         ref.inappropriate.factor == 'Yes,') })
pscam.data$speech.latency.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.speech.0) | ref.speech.4 == 1, NA, ref.speech.0 == 1) })
pscam.data$speech.amount.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.speech.1) | ref.speech.4 == 1, NA, ref.speech.1 == 1) })
pscam.data$spontaneity.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.spontaneity) | ref.spontaneity.factor == 'Unable to Assess', NA,
         ref.spontaneity.factor == 'Yes,') })

## New form
pscam.data$interact.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.interact2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.interact2.factor == 'Yes') })
pscam.data$socsmile.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.socsmile2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.socsmile2.factor == 'Yes') })
pscam.data$peek.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.peek2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.peek2.factor == 'No') })
pscam.data$inappropriate.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.inapp) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.inapp.factor == 'Yes') })
pscam.data$speech.amount.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.speech2) | ref.speech2.factor != 'Present' |
           is.na(ref.speech22.2) | ref.speech22.7 == 1 | is.na(ref.conscious2) | ref.conscious2 > 7,
         NA,
         ref.speech22.2 == 1) })
pscam.data$spontaneity.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.speech2) | ref.speech2.factor != 'Present' |
           is.na(ref.speech22.4) | ref.speech22.7 == 1 | is.na(ref.conscious2) | ref.conscious2 > 7,
         NA,
         ref.speech22.4 == 1) })
pscam.data$speech.latency.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.speech2) | ref.speech2.factor != 'Present' |
           is.na(ref.speech22.6) | ref.speech22.7 == 1 | is.na(ref.conscious2) | ref.conscious2 > 7,
         NA,
         ref.speech22.6 == 1) })

vars.apathy <- Cs(interest.1, inappropriate.1, speech.latency.1, speech.amount.1, spontaneity.1,
                  interact.2, socsmile.2, peek.2, inappropriate.2, speech.amount.2, spontaneity.2,
                  speech.latency.2)

## Lethargy ##
pscam.data$lethargy.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.conscious), NA, ref.conscious > 4) })
pscam.data$lethargy.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.conscious2), NA, ref.conscious2 > 4) })

vars.lethargy <- Cs(lethargy.1, lethargy.2)

## Incoherence ##
## Old form
pscam.data$thought.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.thought) | ref.thought.factor == 'Unable to Assess', NA,
         ref.thought.factor == 'Yes,') })
pscam.data$receptive.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | (!is.na(ref.receptive.3) & ref.receptive.3 == 1), NA,
         (ref.receptive.0 == 1 | ref.receptive.1 == 1 | ref.receptive.2 == 1)) })

## New form
pscam.data$speech.change.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.speech2) | ref.speech2.factor != 'Present' |
           is.na(ref.speech22.1) | ref.speech22.7 == 1 | is.na(ref.conscious2) | ref.conscious2 > 7,
         NA,
         ref.speech22.1 == 1) })
pscam.data$receptive.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         (ref.receptive2.3 == 1 & !is.na(ref.command.reason2.1) & ref.command.reason2.1 == 1)) })

vars.incoherence <- Cs(thought.1, receptive.1, speech.change.2, receptive.2)

## Fluctuation in functioning ##
## Old form
pscam.data$mentalstate.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.mentalstate), NA, ref.mentalstate.factor == 'Acute Change') })
pscam.data$mentalpattern.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.mentalpattern), NA, ref.mentalpattern.factor == 'Fluctuating') })
pscam.data$naps.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.naps) | ref.naps.factor == 'Unable to Assess', NA,
         ref.naps.factor == 'No') })
pscam.data$daynight.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.daynight) | ref.daynight.factor == 'Unable to Assess', NA,
         ref.daynight.factor == 'Yes,') })
pscam.data$nocturnal.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.nocturnal) | ref.nocturnal.factor == 'Unable to Assess', NA,
         ref.nocturnal.factor == 'No') })

## New form
pscam.data$mentalstate.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.mentalstate2), NA, ref.mentalstate2.factor == 'Acute Change') })
pscam.data$mentalpattern.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.mentalpattern2), NA, ref.mentalpattern2.factor == 'Fluctuating') })
pscam.data$naps.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.naps2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.naps2.factor == 'No') })
pscam.data$daynight.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.daynight2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.daynight2.factor == 'Yes') })
pscam.data$nocturnal.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.nocturnal2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.nocturnal2.factor == 'No') })

vars.fluctuation <- Cs(mentalstate.1, mentalpattern.1, naps.1, daynight.1, nocturnal.1,
                       mentalstate.2, mentalpattern.2, naps.2, daynight.2, nocturnal.2)

## Restlessness ##
## Old form
pscam.data$agitated.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.conscious), NA, ref.conscious < 4) })
pscam.data$euphoria.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.euphoria) | ref.euphoria.factor == 'Unable to Assess', NA,
         ref.euphoria.factor == 'Yes,') })

## New form
pscam.data$agitated.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.conscious2), NA, ref.conscious2 < 4) })
pscam.data$energy.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.conscious2) | ref.conscious2 > 7 | is.na(ref.excess2), NA,
         ref.excess2.factor == 'Yes') })

vars.restless <- Cs(agitated.1, euphoria.1, agitated.2, energy.2)

## Hallucinations ##
## Old form
pscam.data$auditory.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.auditory) | ref.auditory.factor == 'Unable to Assess', NA,
         ref.auditory.factor == 'Yes,') })
pscam.data$visual.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.visual) | ref.visual.factor == 'Unable to Assess', NA,
         ref.visual.factor == 'Yes,') })
pscam.data$gust.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.gustatory) | ref.gustatory.factor == 'Unable to Assess', NA,
         ref.gustatory.factor == 'Yes,') })
pscam.data$hyperacusis.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.hyperacusis) | ref.hyperacusis.factor == 'Unable to Assess', NA,
         ref.hyperacusis.factor == 'Yes,') })
pscam.data$disturbance.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.disturbance) | ref.disturbance.factor == 'Unable to Assess', NA,
         ref.disturbance.factor == 'Yes,') })

## New form
pscam.data$hallucinations.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.hallucinations2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.hallucinations2.factor == 'Yes') })
pscam.data$atypresp.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.atypresp2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.atypresp2.factor == 'Yes') })
pscam.data$hyperacusis.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.hyperacusis2) | ref.hyperacusis2.factor == 'N/A' |
           is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.hyperacusis2.factor == 'Yes') })
pscam.data$soothe.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.soothe2) | ref.soothe2.factor == 'N/A' |
           is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.soothe2.factor == 'Yes') })

vars.hallucinations <- Cs(auditory.1, visual.1, gust.1, hyperacusis.1, disturbance.1,
                          hallucinations.2, atypresp.2, hyperacusis.2, soothe.2)

## Anxiety/fear ##
pscam.data$anxiety.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.irritability) | ref.irritability.factor == 'Unable to Assess', NA,
         ref.irritability.factor == 'Yes,') })
pscam.data$anxiety.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.irrang2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.irrang2.factor == 'Yes') })

vars.anxiety <- Cs(anxiety.1, anxiety.2)

## Inconsolability ##
pscam.data$inconsol.1 <- with(pscam.data, {
  ifelse(!del.rater.1 | is.na(ref.inconsolability) | ref.inconsolability.factor == 'Unable to Assess',
         NA,
         ref.inconsolability.factor == 'Yes,') })
pscam.data$inconsol.2 <- with(pscam.data, {
  ifelse(!del.rater.2 | is.na(ref.inconsol2) | is.na(ref.conscious2) | ref.conscious2 > 7, NA,
         ref.inconsol2.factor == 'Yes') })

vars.inconsolability <- Cs(inconsol.1, inconsol.2)

### --- Determine whether each symptom is present: at least one question answered and abnormal --- ###
pscam.data$sustained.attn <-
  ifelse(rowSums(is.na(pscam.data[,vars.sustained.attn])) == length(vars.sustained.attn), NA,
         rowSums(pscam.data[,vars.sustained.attn], na.rm = TRUE) > 0)
pscam.data$shift.attn <-
  ifelse(rowSums(is.na(pscam.data[,vars.shift.attn])) == length(vars.shift.attn), NA,
         rowSums(pscam.data[,vars.shift.attn], na.rm = TRUE) > 0)
pscam.data$attention <-
  ifelse(rowSums(is.na(pscam.data[,vars.attention])) == length(vars.attention), NA,
         rowSums(pscam.data[,vars.attention], na.rm = TRUE) > 0)
pscam.data$orientation <-
  ifelse(rowSums(is.na(pscam.data[,vars.orientation])) == length(vars.orientation), NA,
         rowSums(pscam.data[,vars.orientation], na.rm = TRUE) > 0)
pscam.data$consciousness <-
  ifelse(rowSums(is.na(pscam.data[,vars.consciousness])) == length(vars.consciousness), NA,
         rowSums(pscam.data[,vars.consciousness], na.rm = TRUE) > 0)
pscam.data$apathy <-
  ifelse(rowSums(is.na(pscam.data[,vars.apathy])) == length(vars.apathy), NA,
         rowSums(pscam.data[,vars.apathy], na.rm = TRUE) > 0)
pscam.data$hypokinesia <-
  ifelse(rowSums(is.na(pscam.data[,vars.lethargy])) == length(vars.lethargy), NA,
         rowSums(pscam.data[,vars.lethargy], na.rm = TRUE) > 0)
pscam.data$incoherence <-
  ifelse(rowSums(is.na(pscam.data[,vars.incoherence])) == length(vars.incoherence), NA,
         rowSums(pscam.data[,vars.incoherence], na.rm = TRUE) > 0)
pscam.data$fluctuation <-
  ifelse(rowSums(is.na(pscam.data[,vars.fluctuation])) == length(vars.fluctuation), NA,
         rowSums(pscam.data[,vars.fluctuation], na.rm = TRUE) > 0)
pscam.data$restlessness <-
  ifelse(rowSums(is.na(pscam.data[,vars.restless])) == length(vars.restless), NA,
         rowSums(pscam.data[,vars.restless], na.rm = TRUE) > 0)
pscam.data$hallucination <-
  ifelse(rowSums(is.na(pscam.data[,vars.hallucinations])) == length(vars.hallucinations), NA,
         rowSums(pscam.data[,vars.hallucinations], na.rm = TRUE) > 0)
pscam.data$anxiety <-
  ifelse(rowSums(is.na(pscam.data[,vars.anxiety])) == length(vars.anxiety), NA,
         rowSums(pscam.data[,vars.anxiety], na.rm = TRUE) > 0)
pscam.data$inconsolability <-
  ifelse(rowSums(is.na(pscam.data[,vars.inconsolability])) == length(vars.inconsolability), NA,
         rowSums(pscam.data[,vars.inconsolability], na.rm = TRUE) > 0)

## Delusions not present in peds data; create completely blank variable ##
pscam.data$delusions <- NA

## -- Variable indicating whether patient had delirium per reference rater (combines two forms) -- ##
pscam.data$ref.delirium <- with(pscam.data, {
  ifelse(is.na(ref.delirium.present) &
           (is.na(ref.delirium.present2) | ref.delirium.present2.factor == 'UTA'), NA,
         (!is.na(ref.delirium.present) & ref.delirium.present.factor == 'Yes') |
           (!is.na(ref.delirium.present2) & ref.delirium.present2.factor == 'Yes')) })

pscam.data$del.type <- with(pscam.data, {
  factor(ifelse(is.na(ref.delirium) | !ref.delirium |
                 (is.na(ref.delirium.present.type) & is.na(ref.deliriumtype2)), NA,
         ifelse((!is.na(ref.delirium.present.type) & ref.delirium.present.type.factor == 'HYPOactive') |
                  (!is.na(ref.deliriumtype2) & ref.deliriumtype2.factor == 'Hypoactive'), 1,
         ifelse((!is.na(ref.delirium.present.type) & ref.delirium.present.type.factor == 'HYPERactive') |
                  (!is.na(ref.deliriumtype2) & ref.deliriumtype2.factor == 'Hyperactive'), 2, 3))),
         levels = 1:3, labels = c('Hypoactive', 'Hyperactive', 'Mixed')) })

## Get study day ##
pscam.data$day <- as.numeric(gsub('^Assessment', '', as.character(pscam.data$redcap.event.name.factor)))

## -- Create subset of delirious assessments ONLY to match DSD database, with only relevant variables -- ##
# pscam.del <- subset(pscam.data, !is.na(ref.delirium) & ref.delirium, select = c(id, redcap.event.name.factor, ))
pscam.long <- pscam.data[!is.na(pscam.data$ref.delirium) & pscam.data$ref.delirium,
                         c('id', 'day', 'ref.delirium', 'del.type', vars.attention, vars.orientation,
                           vars.consciousness, vars.apathy, vars.lethargy, vars.incoherence,
                           vars.fluctuation, vars.restless, vars.hallucinations, vars.anxiety,
                           'sustained.attn', 'shift.attn', 'attention', 'orientation',
                           'consciousness', 'apathy', 'hypokinesia', 'incoherence', 'fluctuation',
                           'restlessness', 'delusions', 'hallucination', 'anxiety',
                           'inconsolability')]

## Get number of delirious assessments per patient ##
num.del.asmts <- pscam.long %>%
  group_by(id) %>%
  summarise(del.asmts = sum(ref.delirium),
            del.asmts.hypo = sum(del.type == 'Hypoactive', na.rm = TRUE),
            del.asmts.hyper = sum(del.type == 'Hyperactive', na.rm = TRUE),
            del.asmts.mixed = sum(del.type == 'Mixed', na.rm = TRUE))

pscam.oneobs <- merge(pscam.oneobs, num.del.asmts, by = 'id', all = TRUE)

## Restrict oneobs data to patients with delirium assessments, needed variables ##
pscam.oneobs <- pscam.oneobs %>%
  filter(!is.na(del.asmts)) %>%
  select(id, age, los, starts_with('admit.dx'), prism.score, dc.status, contains('del.asmts')) %>%
  mutate(prism.cat = factor(ifelse(is.na(prism.score), NA,
                            ifelse(prism.score <= 10, 1,
                            ifelse(prism.score <= 20, 2, 3))),
                            levels = 1:3,
                            labels = c('Low risk (<=10)',
                                       'Moderate risk (11-20)',
                                       'High risk (<20)')))

rm(list = Cs(mult.ids, num.del.asmts, pscam.data, curl_handle, prism.vars,
             ref.rater.no, vars.anxiety, vars.apathy, vars.attention, vars.consciousness,
             vars.fluctuation, vars.hallucinations, vars.incoherence, vars.inconsolability,
             vars.lethargy, vars.orientation, vars.restless, vars.shift.attn, vars.sustained.attn,
             get.api.data, make.date))