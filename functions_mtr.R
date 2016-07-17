library(magrittr) # pipe operator
library(rjson) # read files in json format
library(lubridate) # dates
library(readr) # read from files
library(stringr) # string manipulation
library(plyr) # matrices and lists manipulation
library(dplyr) # data frame manipulation
library(tidyr) # tidy the data set
library(data.table) # bind multiple data frames

from_raw <- function(name, replacement) {
    ## read a raw data file, manipulate, and output a data frame
    df <- readr::read_csv(name) %>% 
        mutate(created_at = timestamp %>% 
                   as.POSIXct(format = "%Y-%m-%d %H:%M:%S %z", ### no need origin when specify format
                              tz = "Asia/Hong_Kong"), ### change to correct time zone
               text = text %>% 
                   stringr::str_to_title() %>% ### To Capitalised Titles
                   stringr::str_replace_all(replacement)) %>% 
        select(text, created_at) %>% 
        arrange(created_at)
}

operational_date <- function(df) {
    ## deduce the operational date from created_at, which is slightly different to calendar date, because of the running times of MTR
    df1 <- df %>% 
        mutate(hour = hour(created_at), 
               op_date = as.Date(created_at, tz = "Asia/Hong_Kong"))
    df2 <- df1 %>% filter(hour < 5) %>% ### op_date should be prev day
        mutate(op_date = op_date - 1)
    df3 <- df1 %>% filter(hour >= 5)
    bind_rows(df2, df3) %>% 
        arrange(created_at) %>% 
        select(-hour)
}

filter_incident_with_keywords <- function(df) {
    df %>%
        mutate(logi.signal = text %>% str_detect("Signal Failure") | text %>% str_detect("Signaling Fault") | text %>% str_detect("Point Failure"),
               logi.train = text %>% str_detect("Faulty Train"),
               logi.power = text %>% str_detect("Power Failure"),
               logi.track = text %>% str_detect("Faulty Track") | text %>% str_detect("Cracked Track"),
               logi.door = text %>% str_detect("Door") | text %>% str_detect("Platform Gate"),
               logi.electric = text %>% str_detect("Irregular Electricity"),
               logi.engine = text %>% str_detect("Engineering Works"),
               logi.unknown = text %>% str_detect("Due To") & text %>% str_detect("Unknown"),
               logi.delay = text %>% str_detect("Major Delay") | text %>% str_detect("Minor Delay") | text %>% str_detect("Severe Delay"), # may be useful
               logi.wtoo = text %>% str_detect("Wrong Type Of Object"),
               logi.exit = text %>% str_detect("Exit"),
               logi.dog = text %>% str_detect("Dog")) %>%
        filter((logi.signal | logi.train | logi.power | logi.track | logi.door | logi.electric | logi.engine | logi.unknown) & !(logi.wtoo | logi.exit | logi.dog)) %>%
        tidyr::gather(key = type, value = boolean, # line 1 to get type
                      logi.signal, logi.train, logi.power, logi.track, logi.door, logi.electric, logi.engine, logi.unknown) %>%
        filter(boolean) %>% # line 2
        select(-boolean, -contains("logi.")) %>%
        mutate(type = type %>% str_sub(6, -1))
}

extract_incident <- function(df) {
    ## identify the relevant tweets
    df %>% 
        mutate(logi.line = outer(text, lines, str_detect) %>% apply(1, any),
               logi.station = outer(text, stations, str_detect) %>% apply(1, any),
               logi.time.start = text %>% str_detect("^[0-9]")) %>% # w/ timestamps at the beginning
        filter((logi.line | logi.station) & logi.time.start) %>%
        select(text, created_at, op_date) %>% 
        filter_incident_with_keywords %>% 
        arrange(created_at)
}

remove_use_or_using <- function(df) {
    ## remove anything after use or using (in incident tweets) to avoid confusion when looking up
    df1 <- df %>% 
        filter(text %>% str_detect("Use")) %>%
        mutate(position = text %>% str_locate("Use") %>% magrittr::extract(,1),
               text = text %>% str_sub(1L, position - 2L)) %>% 
        select(-position)
    df2 <- df %>% 
        filter(!(text %>% str_detect("Use")))
    df <- bind_rows(df1, df2)
    df1 <- df %>% 
        filter(text %>% str_detect("Using")) %>% 
        mutate(position = text %>% str_locate("Using") %>% magrittr::extract(,1),
               text = text %>% str_sub(1L, position - 2L)) %>%
        select(-position)
    df2 <- df %>% 
        filter(!(text %>% str_detect("Using")))
    bind_rows(df1, df2)
}    

extract_good_service <- function(df) {
    ## identify the tweets with good service
    df %>% 
        mutate(logi.line = outer(text, lines, str_detect) %>% apply(1, any),
               logi.station = outer(text, stations, str_detect) %>% apply(1, any),
               logi.time.start = text %>% str_detect("^[0-9]"),
               logi.good = text %>% str_detect("Good Service"),
               logi.petition = text %>% str_detect("Petition")) %>%
        filter((logi.line | logi.station) & logi.time.start & logi.good & !logi.petition) %>%
        ## has a line/station, starts with a time stamp, includes "good service", and excludes "petition"
        select(-starts_with("logi")) %>%
        arrange(created_at)
}

occurrence_time <- function(df) {
    ## find the time of occurrence, and split original tweet into time and text
    ## op_date not to be used because of complicated relationship b/w occurred_at and created_at
    df %>% tbl_df %>%
        mutate(colon = text %>% str_sub(3L, 3L) == ":",
               pos = ifelse(colon, 6L, 5L), # position of splitter
               hour = text %>% str_sub(1L, 2L), 
               minute = text %>% str_sub(pos - 2L, pos - 1L),
               text = text %>% str_sub(pos + 1L, -1L) %>% str_trim,
               occurred_at0 = update(created_at, hour = hour, minute = minute, second = 0),
               timediff0 = as.numeric(created_at - occurred_at0, unit = "mins"),
               occurred_at = ifelse(timediff0 < -100, occurred_at0 - duration(1, "day"), occurred_at0) %>%
                   as.POSIXct(origin = as.Date("1970-01-01"),
                              tz = "Asia/Hong_Kong"), # improve coding here!
               timediff = as.numeric(created_at - occurred_at, unit = "mins")) %>%
        select(-colon, -pos, -occurred_at0, -timediff0, -hour, -minute)
}

operational_date1 <- function(df) {
    ## this can be used to check if the same op_date is deduced by created_at or occurred_at
    df1 <- df %>% 
        mutate(hour = hour(occurred_at), 
               op_date = as.Date(occurred_at, tz = "Asia/Hong_Kong"))
    df2 <- df1 %>% filter(hour < 5) %>% ### op_date should be prev day
        mutate(op_date = op_date - 1)
    df3 <- df1 %>% filter(hour >= 5)
    bind_rows(df2, df3) %>% 
        arrange(occurred_at) %>% 
        select(-hour)
}

augment_plus_lines <- function(df) {
    ## augment the tweets with "A + B Lines" to "A Line + B Lines" for detectability
    df %>% 
        mutate(text = ifelse(str_count(text, "\\+") == 1 & str_detect(text, "Lines"),
                             str_replace(text, "\\+", "Line \\+"),
                             text))
}

extract_impute_line <- function(df) {
    ## extract or impute the lines corresponding to each relevant record
    list.line <- outer(df$text, lines, str_detect) %>% plyr::alply(1, function(x) lines[x]) ### is there a better way?
    ## the i-th element in the list a vector of stations in the i-th tweet (if any)
    station.lengths <- list.line %>% sapply(length) ### each element is length of corresponding vector of stations
    if (all(station.lengths > 0)) { ### all tweets have line(s) extracted
        df0 <- df %>% mutate(line = list.line)
    }
    else {
        subset.with <- df %>% mutate(line = list.line) %>% filter(station.lengths > 0)
        subset.wout <- df %>% filter(station.lengths == 0)
        list.station <- outer(subset.wout$text, stations, str_detect) %>% plyr::alply(1, function(x) stations[x])
        imputed.line <- list.station %>% lapply(function(x) data_frame(station = x)) %>% 
            lapply(left_join, both, by = "station") %>% 
            lapply(count, line) %>% 
            lapply(top_n, 1) %>% 
            lapply(select, line) %>% 
            lapply(unlist, use.names = F) ### again, is there a better way?
        subset.wout %<>% mutate(line = imputed.line)
        df0 <- bind_rows(subset.with, subset.wout) %>% 
            arrange(created_at)
    }
    df0 
}

interpose_good_service <- function(df) {
    ## from a matched data frame to one with incident tweets and good service tweets one after the other
    ## assuming that the tweets all have the same op_date & line
    df0 <- df %>% # remove records with both tweets being good service but with different time stamps
        arrange(occurred_at.x, occurred_at.y) %>%
        distinct %>% 
        filter(!(str_detect(text.x, "Good Service") & str_detect(text.y, "Good Service") & occurred_at.x != occurred_at.y))
    df1 <- df0 %>% filter(!is.na(occurred_at.y))
    df2 <- df0 %>% filter(is.na(occurred_at.y)) %>% slice(1) # just first row
    ind0 <- str_detect(df1$text.x, "Good Service") %>% which # indices
    m <- ind0 %>% length
    if (m == 0) {# all are incident tweets
        df3 <- df1
    }
    else {
        n <- df1 %>% nrow
        ind1 <- c(1, ind0[-m] + 1)
        if (ind0[m] != n) # the last index is not n i.e. there are orphan tweets at the end
            ind1 <- c(ind1, (ind0[m] + 1):n)
        ind1 %<>% unique
        df3 <- df1 %>% slice(ind1)
    }
    bind_rows(df2, df3)
}

dnbinom_rqk <- function(x, r1, r2, q1, q2, k, n) {
    ## pmf of incident count in 1-chgpt model, in r & q simultaneously; vectorised over all model pars
    dnbinom(x, r1, q1) * k / n + dnbinom(x, r2, q2) * (n - k) /n
}

