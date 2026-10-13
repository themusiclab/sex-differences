### data preparation and recoding
require(tidyverse)
is.ordinal <- function(d, ncat=7){
  o.vars <- NULL
  for(i in seq(along=colnames(d))){
    if(length(table(d[,i])) <= ncat){
      o.vars[i] <- names(d)[i]}
  }
  o.vars <- as.character(na.omit(o.vars))
}

data_prep <- function(l.age, u.age, lessons.l.age, lessons.u.age, repeats, HI){
  dat <- readRDS("../data/clean-person-level-data-filtered.rds")

  dat <- dat %>%  mutate(
    lessons_enjoy_n = case_when(
      lessons_enjoy == "I disliked them a lot" ~ 1L,
      lessons_enjoy =="I disliked them a little" ~ 2L,
      lessons_enjoy =="I liked them a little" ~ 3L,
      lessons_enjoy =="I liked them a lot" ~ 4L))

  dat <- dat %>%  mutate(lessons_peers_n = case_when(  
    lessons_peers == "They were a lot better than me" ~ 1L, 
    lessons_peers == "They were a little better than me" ~ 2L,
    lessons_peers == "We had about the same level of musical skills" ~ 3L,
    lessons_peers == "I was a little better than them" ~ 4L,
    lessons_peers == "I was a lot better than them" ~ 5L))

  dat <- dat %>% mutate(music_tap_n = case_when(   
    music_tap == "No" ~ 0L, 
    music_tap == "Yes" ~ 1L,
    music_tap == "I'm not sure" ~ NA_integer_))

  dat <- dat %>% mutate(out_of_tune_n = case_when(  
    out_of_tune == "No" ~ 0L, 
    out_of_tune == "Yes" ~ 1L,
    out_of_tune == "I'm not sure" ~ as.integer(NA)))

  dat <- dat %>% mutate(music_lessons_n = case_when(  
    music_lessons == "No" ~ 0L, 
    music_lessons == "Yes" ~ 1L))

  dat <- dat %>% mutate(music_listen_time_n = case_when(  
    music_listen_time == "No time at all" ~ 1L, 
    music_listen_time == "1-5 minutes" ~ 2L,
    music_listen_time == "6-10 minutes" ~ 3L,
    music_listen_time == "11-15 minutes" ~ 4L,  
    music_listen_time == "16-30 minutes" ~ 5L,
    music_listen_time == "31-60 minutes" ~ 6L,
    music_listen_time == "1-2 hours" ~ 7L,
    music_listen_time == "2-4 hours" ~ 8L,
    music_listen_time == "More than 4 hours" ~ 9L))

  dat <- dat %>% mutate(music_make_time_n = case_when(  
    music_make_time == "No time at all" ~ 1L, 
    music_make_time == "1-5 minutes" ~ 2L,
    music_make_time == "6-10 minutes" ~ 3L,
    music_make_time == "11-15 minutes" ~ 4L,  
    music_make_time == "16-30 minutes" ~ 5L,
    music_make_time == "31-60 minutes" ~ 6L,
    music_make_time == "1-2 hours" ~ 7L,
    music_make_time == "2-4 hours" ~ 8L,
    music_make_time == "More than 4 hours" ~ 9L))

  dat <- dat %>% mutate(music_skill_n = case_when( 
    music_skill == "I have no skill at all" ~ 1L, 
    music_skill == "I'm a novice" ~ 2L,
    music_skill == "I have some skill" ~ 3L,
    music_skill == "I have a lot of skill" ~ 4L,
    music_skill == "I'm an expert" ~ 5L))

  dat <- dat %>% mutate(music_p_sing_n = case_when( 
    music_p_sing == "Once every 3 days or less" ~ 1L, 
    music_p_sing == "Once every day or two" ~ 2L,
    music_p_sing == "2-3 times a day" ~ 3L,
    music_p_sing == "4-7 times a day" ~ 4L,
    music_p_sing == "8 or more times a day" ~ 5L,
    music_p_sing == "I don't know" ~ NA_integer_))

  dat <- dat %>% mutate(music_feel_n = case_when( 
    music_feel == "No, never" ~ 1L, 
    music_feel == "Yes, but rarely" ~ 2L,
    music_feel == "Yes, sometimes" ~ 3L,
    music_feel == "Yes, often" ~ 4L,
    music_feel == "I’m not sure" ~ NA_integer_))

  dat <- dat %>% mutate(lessons_why = case_when( 
    lessons_why == "A parent made me start music lessons." ~ "Parent", 
    lessons_why ==  "My friends were starting music lessons so I did too." ~ "Friends",
    lessons_why ==  "Music lessons were a school requirement." ~ "School",
    lessons_why ==  "I really wanted to start music lessons." ~ "Myself"))

  dat <- dat %>% mutate(education_n = case_when( 
    education == "Some elementary/middle school (primary school)" ~ 1L, 
    education == "Completed elementary/middle school (primary school)" ~ 2L,
    education == "Some high school" ~ 3L,
    education == "Completed high school (secondary school)" ~ 4L,
    education == "Some undergrad (higher education)" ~ 5L,
    education == "Completed undergrad degree (~3-5 years higher education)" ~ 6L,
    education == "Some graduate school" ~ 7L,
    education == "Completed graduate school" ~ 8L))

  dat <- dat %>% mutate(income_n = case_when( 
    income == "Under $10,000" ~ 1L, 
    income == "$10,000 to $19,999" ~ 2L,
    income == "$20,000 to $29,999" ~ 3L,
    income == "$30,000 to $39,999" ~ 4L,
    income == "$40,000 to $49,999" ~ 5L,
    income == "$50,000 to $74,999" ~ 6L,
    income == "$75,000 to $99,999" ~ 7L,
    income == "$100,000 to $150,000" ~ 8L,
    income == "Over $150,000" ~ 9L,
    income == "I'd prefer not to say" ~ NA_integer_))

  dat <- dat %>% mutate(theory_training_n = case_when( 
    theory_training == "No music theory training" ~ 1L, 
    theory_training == "A little music theory training" ~ 2L,
    theory_training == "Some music theory training" ~ 3L,
    theory_training == "A moderate amount of music theory training" ~ 3L,
    theory_training == "A lot of music theory training" ~ 4L))

  dat <- dat %>% mutate(familiarity_tradmusic_n = case_when( 
    familiarity_tradmusic == "I've never heard traditional music" ~ 1L, 
    familiarity_tradmusic == "I'm a little familiar with traditional music" ~ 2L,
    familiarity_tradmusic == "I'm somewhat familiar with traditional music" ~ 3L,
    familiarity_tradmusic == "I'm very familiar with traditional music" ~ 4L))

  dat <- dat %>% mutate(workspace_n = case_when( 
    workspace == "I am in a very noisy place" ~ 1L, 
    workspace == "I am in a somewhat noisy place" ~ 2L,
    workspace == "I am in a somewhat quiet place" ~ 3L,
    workspace == "I am in a very quiet place" ~ 4L))

  dat <- dat %>% mutate(hearing_imp_n = case_when( 
    hearing_imp == "Yes" ~ 0L, 
    hearing_imp == "No" ~ 1L,
    hearing_imp == "I don't know" ~ NA_integer_))

  dat$music_listen_n <- cut_interval(dat$music_listen, n=6)
  levels(dat$music_listen_n) <- 1:6
  dat$music_listen_n <- as.numeric(dat$music_listen_n)
  
  dat$gender_n <- as.factor(dat$gender)
  
  dat <- dat[dat$age>=l.age & dat$age<= u.age,]
  
  dat$lessons_age_n <- dat$lessons_age
  dat$lessons_age_n[dat$lessons_age_n > dat$age] <- NA
  
  lessons.u.age <- dat$age
  dat$lessons_age_n[dat$lessons_age_n < lessons.l.age | dat$lessons_age_n > lessons.u.age] <- NA
  
  dat <- dat[dat$hearing_imp != HI, ]
  dat <- dat[dat$taken_before == repeats, ]
  
  dat <- as.data.frame(dat[!is.na(dat$cabat_ability) & !is.na(dat$mdt_ability) & !is.na(dat$mpt_ability), c("user_id","country","lessons_enjoy_n",
                                                                                           "lessons_peers_n",
                                                                                           "theory_training_n",
                                                                                           "music_tap_n",
                                                                                           "out_of_tune_n",
                                                                                           "music_lessons_n",
                                                                                           "music_skill_n",
                                                                                           "music_listen_n",
                                                                                           "lessons_age",
                                                                                           "lessons_age_n",
                                                                                           "music_p_sing_n",
                                                                                           "music_feel_n",
                                                                                           "music_listen_time_n",
                                                                                           "familiarity_tradmusic_n",
                                                                                           "education_n",
                                                                                           "income_n",
                                                                                           "gender_n",
                                                                                           "age",
                                                                                           "hearing_imp_n",
                                                                                           "workspace_n",
                                                                                           "cabat_ability",
                                                                                           "mdt_ability",
                                                                                           "mpt_ability")])
save(dat,file="../results/dat.RData")
dat
}
