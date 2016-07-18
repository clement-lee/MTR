### 00) workspace preliminary
options(width = 200) # changes when console screen does so, on windows
remove(list = ls())
print(base::date())
options(width = 150)
source("functions_mtr.R")
### functions assumed to be dplyr/base/self-written unless otherwise specified



### 01) data preliminary
alias <- c("Kowloon Tong" = "Kow Tong", 
           "East Tsim Sha Tsui" = "East TST", 
           "Mong Kok East" = "MK East", 
           "Kowloon Bay" = "Kow Bay", 
           "Tsuen Wan West" = "TW West",
           "Sha Tin Wai" = "ST Wai") # necessary changes - see from_raw()
both <- readr::read_csv("stations_and_lines.csv") %>% ### data frame of lines and stations
    mutate(station = station %>% 
               stringr::str_to_title() %>% 
               stringr::str_replace_all(alias),
           line = line %>% 
               stringr::str_to_title())
stations <- both$station %>% unique
lines <- c(both$line %>% unique, "Light Rail")
typo <- c("Tsueng Kwan O" = "Tseung Kwan O",
          "Ma On Than" = "Ma On Shan",
          "Kuwn Tong" = "Kwun Tong",
          "Sinal" = "Signal",
          "Final Failure" = "Signal Failure",
          "Signals Failure" = "Signal Failure",
          "Powet" = "Power",
          "Fault Track" = "Faulty Track")
df.lines <- readr::read_csv("lines.csv") %>%
    select(line = 1, line.chi = 2) # for lines in chinese
df.types <- readr::read_csv("types.csv") %>%
    select(type = 1, type.chi = 2) # for types in chinese



### 02) read in the data
df1 <- from_raw("tweets.csv", c(alias, typo)) %>% 
    operational_date # based on created_at



### 03) extract tweets with incidents
df2.incidents <- df1 %>% 
    extract_incident %>% 
    remove_use_or_using %>% 
    augment_plus_lines %>% 
    extract_impute_line %>% 
    occurrence_time %>% # why not earlier? coz we need tweets WITH timestamps
    tidyr::unnest(line) # for joining, see below



### 04) extract good service tweets
df2.gdservice <- df1 %>% 
    extract_good_service %>% 
    augment_plus_lines %>% 
    extract_impute_line %>% 
    occurrence_time %>% 
    tidyr::unnest(line) # a list can't be the by variable in joining



### 05) check that op_date is same whether derived from created_at or occurred_at
df2.both <- df2.incidents %>% 
    select(-type) %>%
    bind_rows(df2.gdservice) %>%
    distinct %>%
    left_join(df2.incidents, by = c("text", "created_at", "op_date", "occurred_at", "timediff", "line")) # join all just to make sure
df2.both %>%
    select(-op_date) %>% # remove the one based on created_at
    operational_date1 %>% # based on occurred_at
    all.equal(df2.both) # should be true
## works for tbl_df even if order of columns and/or rows is diff.!



### 06) match incidents with good service
df3.match <- df2.both %>%
    left_join(df2.gdservice, by = c("op_date", "line")) %>%
    mutate(downtime = as.numeric(occurred_at.y - occurred_at.x, unit = "mins")) %>% # time diff. b/w tweets
    filter(is.na(created_at.y) | downtime >= 0) # either no match or good service tweet occurring AFTER incident tweet



### 07) further manipulation after joining
### (2016-04-25) closest to what we want: 1 incident-good service pair per record.
df4.match <- df3.match %>%
    ## i) nest/unnest to group/remove in-between tweets (next 3 lines)
    tidyr::nest(-op_date, -line) %>%
    mutate(filtered = sapply(data, interpose_good_service, simplify = F)) %>%
    tidyr::unnest(filtered) %>%
    ## ii) remove standalone good service tweets in text.x (next 3 lines)
    mutate(text = text.x) %>% # for function next line
    select(-type) %>% 
    filter_incident_with_keywords %>% # checked, will get the same types back
    select(-text) %>%
    ## iii) nest line & create new line variable to accommodate AE + TC 
    tidyr::nest(line) %>%
    mutate(check = data %>% sapply(all.equal, data_frame(line = c("Tung Chung Line", "Airport Express"))) %>% as.character,
           line = ifelse(check == "TRUE", list(data_frame(line = c("Tung Chung Line + Airport Express"))), data)) %>%
    select(-data, -check) %>% 
    tidyr::unnest(line) %>%
    ## iv) dates manipulation - for modelling & visualisaion
    mutate(month = op_date %>% cut(breaks = "month") %>% as.Date, # alt: month1 = op_date %>% lubridate::floor_date("month")
           month.mid = month + days(14), # mid-month, quick fix for position adjustment
           week = op_date %>% cut(breaks = "week") %>% as.Date) %>% 
    left_join(df.lines, by = "line") %>%
    left_join(df.types, by = "type")
write_csv(df4.match, "incidents.csv")
## Cross check with df3.match for any future changes, to make sure there are no unwanted omission.
## Issues: incidents with NA good service tweets can't have downtime computed.



### 08) examine records with same op_date & text.x which are yet nested
df4.match %>% 
    count(op_date, text.x) %>% 
    filter(n > 1) %>% 
    select(-n) %>% 
    inner_join(df4.match, by = c("op_date", "text.x"))
## interestingly, most are genuinely separate incidents (same type & line!)



### 09) weekly number of incidents (can do the same for monthly counts)
df0.weeks <- seq(df4.match$week %>% min, df4.match$week %>% max, by = 7) %>%
    data_frame(week = .) %>% # complete list of weeks
    ## i) all
    left_join(df4.match %>% count(week), by = "week") %>%
    mutate(all = ifelse(is.na(n), 0L, n)) %>%
    select(-n) %>%
    ## ii) signal
    left_join(df4.match %>% filter(type == "signal") %>% count(week), by = "week") %>% 
    mutate(signal = ifelse(is.na(n), 0L, n)) %>%
    select(-n) %>% 
    ## iii) train
    left_join(df4.match %>% filter(type == "train") %>% count(week), by = "week") %>%
    mutate(train = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>%
    ## iv) East Rail Line
    left_join(df4.match %>% filter(line == "East Rail Line") %>% count(week), by = "week") %>% 
    mutate(east.rail = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>% 
    ## v) Kwun Tong Line
    left_join(df4.match %>% filter(line == "Kwun Tong Line") %>% count(week), by = "week") %>% 
    mutate(kwun.tong = ifelse(is.na(n), 0L, n)) %>% 
    select(-n)
## write a function to do the above!
write_csv(df0.weeks, "ts_weeks.csv") # move if we model weekly counts
n0.weeks <- df0.weeks %>% nrow



### 10) monthly number of incidents (similar to df0.weeks)
n0.months <- ((difftime(df4.match$month %>% max,
                      df4.match$month %>% min, units = "days") %>% 
    as.numeric) / (365 / 12)) %>% round + 1
df0.months <- (df4.match$month %>% min + months(0:(n0.months - 1))) %>% 
    data_frame(month = .) %>% # complete list of months
    left_join(df4.match %>% count(month), by = "month") %>% 
    mutate(all = ifelse(is.na(n), 0L, n)) %>% 
    select(-n)
write_csv(df0.months, "ts_months.csv") # move if we model monthly counts



### 11) daily number of incidents (similar to df0.weeks & df0.months)
df0.days <- seq(df4.match$op_date %>% min, df4.match$op_date %>% max, by = 1L) %>% 
    data_frame(op_date = .) %>% ## complete vector of days
    ## i) all
    left_join(df4.match %>% count(op_date), by = "op_date") %>%
    mutate(all = ifelse(is.na(n), 0L, n)) %>%
    select(-n) %>%
    ## ii) signal
    left_join(df4.match %>% filter(type == "signal") %>% count(op_date), by = "op_date") %>% 
    mutate(signal = ifelse(is.na(n), 0L, n)) %>%
    select(-n) %>% 
    ## iii) train
    left_join(df4.match %>% filter(type == "train") %>% count(op_date), by = "op_date") %>%
    mutate(train = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>%
    ## iv) East Rail Line
    left_join(df4.match %>% filter(line == "East Rail Line") %>% count(op_date), by = "op_date") %>% 
    mutate(east.rail = ifelse(is.na(n), 0L, n)) %>% 
    select(-n) %>% 
    ## v) Kwun Tong Line
    left_join(df4.match %>% filter(line == "Kwun Tong Line") %>% count(op_date), by = "op_date") %>% 
    mutate(kwun.tong = ifelse(is.na(n), 0L, n)) %>% 
    select(-n)
n0.days <- df0.days %>% nrow












