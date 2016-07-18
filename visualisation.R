### Principle: we shouldn't do any further data manipulation like counting here. 
### Any summaries can and should be done via stat functions in ggplot2.
### Unless otherwise specified, any function used is in ggplot2/base.
### This means even dplyr functions have to be specified & kept to minimum.

remove(list = ls())
#source("modelling.R") # run the whole script, which in turn runs extraction.R
if (TRUE) { # load only necessity instead of running extraction & modelling
    source("functions_mtr.R")
    df4.match <- read_csv("incidents.csv") # from extraction
    df1.days <- read_csv("ts_days.csv") # from modelling
    x <- df1.days$all
    t0 <- seq_along(x)
    df0.weeks <- read_csv("ts_weeks.csv") # from extraction
    df0.months <- read_csv("ts_months.csv") # from extraction
    df1.counts <- read_csv("counts_days.csv") # from modelling
}
library(ggplot2)
theme_set(theme_bw(24))
dev.new(width = 5.97, height = 7.45)





### 01) Most basic bar charts
## bar chart by type
df4.match %>% 
    ggplot(aes(type)) + 
    geom_bar() + 
    labs(title = "Number of incidents by type",
         y = "Number of incidents") + 
    coord_flip()

## bar chart by type
jpeg(filename = "by_type_chi.jpg", 1920, 1080, res = 200, quality = 100)
df4.match %>% 
    ggplot(aes(type.chi)) + 
    geom_bar() +
    theme(axis.title.y = element_blank()) + # y is AFTER flipping
    labs(y = "\u4E8B\u6545\u5B97\u6578") + # y is before flipping
    coord_flip()
dev.off()

## bar chart by line
df4.match %>% 
    ggplot() + 
    geom_bar(aes(line)) + 
    labs(title = "Number of incidents by line",
         y = "Number of incidents") +
    coord_flip()

### bar chart by line (in chinese)
jpeg(filename = "by_line_chi.jpg", 1920, 1080, res = 200, quality = 100)
df4.match %>%
    ggplot(aes(line.chi)) +
    geom_bar() +
    theme(axis.title.y = element_blank()) + # y is AFTER flipping
    labs(y = "\u4E8B\u6545\u5B97\u6578") + # y is before flipping
    coord_flip()
dev.off()





### 02) Overall temporal analysis
## basic with graphics::hist
df4.match %>% 
    magrittr::extract2("month.mid") %>% 
    hist(., breaks = unique(.))

## histogram over time
df4.match %>% 
    ggplot(aes(month.mid)) + 
    geom_histogram(bins = 63) + 
    labs(x = "Time", y = "Count")
## dates are continuous (even though we discretise them), 
## so use histogram, not bar chart (the latter will fail when fill/stacking)

## better still, a time series plot
df0.months %>% 
    ggplot(aes(month, all)) + 
    geom_line() + 
    geom_point() + 
    labs(title = "Monthly number of all incidents",
         x = "Time", y = "Number of incidents")
df0.weeks %>% # monthly counterpart
    ggplot(aes(week, all)) + 
    geom_line() + 
    geom_point() +
    labs(title = "Weekly number of all incidents",
         x = "Time", y = "Number of incidents")
df1.days %>% # essentially df1.days + w/ model estimates
    ggplot(aes(op_date, all)) + 
    geom_line() + 
    geom_point() +
    labs(title = "Daily number of all incidents",
         x = "Time", y = "Number of incidents")

## time series plots (in chinese)
jpeg(filename = "by_month_chi.jpg", 1920, 1080, res = 200, quality = 100)
df0.months %>% 
    ggplot(aes(month, all)) + 
    geom_line() + 
    geom_point() + 
    labs(x = "\u5E74\u4EFD", y = "\u4E8B\u6545\u5B97\u6578")
dev.off()

jpeg(filename = "by_week_chi.jpg", 1920, 1080, res = 200, quality = 100)
df0.weeks %>% 
    ggplot(aes(week, all)) + 
    geom_line() + 
    geom_point() +
    labs(x = "\u5E74\u4EFD", y = "\u4E8B\u6545\u5B97\u6578")
dev.off()

jpeg(filename = "by_day_chi.jpg", 1920, 1080, res = 200, quality = 100)
df1.days %>% # essentially df1.days + w/ model estimates
    ggplot(aes(op_date, all)) + 
    geom_line() + 
    geom_point() +
    labs(x = "\u5E74\u4EFD", y = "\u4E8B\u6545\u5B97\u6578")
dev.off()

## ACF (to show independence)
df0.months %>% 
    magrittr::extract2("all") %>% 
    acf(50) # slight yearly effect
df0.weeks %>% 
    magrittr::extract2("all") %>%
    acf(100)
df1.days %>% 
    magrittr::extract2("all") %>% 
    acf(1000)





### 03) Type vs time: monthly number of incidents by type
## raster/tile/rectangle
df4.match %>% 
    dplyr::count(month.mid, type) %>% 
    ggplot(aes(month.mid, type, fill = n)) + 
    geom_tile() + # geom_raster() fails big here
    scale_fill_gradientn(name = NULL,
                         colours = topo.colors(12) %>% rev,
                         breaks = c(1, 4, 7, 10, 13)) + 
    labs(title = "Monthly number of incidents by type", 
         x = "Time", y = "Type") +
    coord_cartesian(xlim = c(as.Date("2011-03-01"), as.Date("2016-02-01")))

## raster/tile/rectangle (in chinese)
jpeg(filename = "month_by_type_chi.jpg", 1920, 1080, res = 200, quality = 100)
df4.match %>% 
    dplyr::count(month.mid, type.chi) %>% 
    ggplot(aes(month.mid, type.chi, fill = n)) + 
    geom_tile() +
    scale_fill_gradientn(name = NULL,
                         colours = topo.colors(12) %>% rev,
                         breaks = c(1, 4, 7, 10, 13)) + 
    labs(x = "\u5E74\u4EFD") +
    theme(axis.title.y = element_blank()) + # y is AFTER flipping
    coord_cartesian(xlim = c(as.Date("2011-03-01"), as.Date("2016-02-01")))
dev.off()

## stacked histogram over time
df4.match %>% 
    ggplot(aes(month.mid)) + 
    geom_histogram(aes(fill = type), bins = 63) + # geom_bar() fails big here
    labs(x = "Time", y = "Count") # not very easily interpretable






### 04) Line vs time: monthly number of incidents by line
## raster/tile/rectangle
df4.match %>% 
    dplyr::count(month.mid, line) %>%
    ggplot(aes(month.mid, line)) + # flexible location of fill
    geom_tile(aes(fill = n)) +
    scale_fill_gradientn(name = NULL,
                         colours = topo.colors(12) %>% rev) + 
    labs(title = "Monthly number of incidents by line",
         x = "Time", y = "Line") +
    coord_cartesian(xlim = c(as.Date("2011-03-01"), as.Date("2016-02-01")))

## raster/tile/rectangle by line (in chinese)
jpeg(filename = "month_by_line_chi.jpg", 1920, 1080, res = 200, quality = 100)
df4.match %>% 
    dplyr::count(month.mid, line.chi) %>%
    ggplot(aes(month.mid, line.chi)) + # flexible location of fill
    geom_tile(aes(fill = n)) +
    scale_fill_gradientn(name = NULL,
                         colours = topo.colors(12) %>% rev) + 
    labs(x = "\u5E74\u4EFD") +
    theme(axis.title.y = element_blank()) + # y is AFTER flipping
    coord_cartesian(xlim = c(as.Date("2011-03-01"), as.Date("2016-02-01")))
dev.off()

## stacked histogram over time
df4.match %>% 
    ggplot(aes(month.mid)) + 
    geom_histogram(aes(fill = line), bins = 63) +
    labs(x = "Time", y = "Count") # not very easily interpretable





### 05) efficiency
## of MTR
df4.match %>% 
    dplyr::filter(downtime != 0 & !is.na(downtime)) %>%
    ggplot(aes(downtime)) + 
    geom_histogram(binwidth = 10, center = 10 / 2) + 
    coord_cartesian(xlim = c(-10, 800)) + 
    labs(title = "Histogram of incident downtime",
         x = "Minutes", y = "Number of incidents")

## of MTR (in chinese)
jpeg(filename = "downtime_chi.jpg", 1920, 1080, res = 200, quality = 100)
df4.match %>% 
    dplyr::filter(downtime != 0 & !is.na(downtime)) %>%
    ggplot(aes(downtime)) + 
    geom_histogram(binwidth = 10, center = 10 / 2) + 
    coord_cartesian(xlim = c(-10, 800)) + 
    labs(x = "\u5206\u9418", y = "\u4E8B\u6545\u5B97\u6578")
dev.off()

## of mtrupdate
df4.match %>% 
    ggplot(aes(timediff.x)) +  
    geom_histogram(binwidth = 2, center = 2 / 2) + 
    coord_cartesian(xlim = c(-10, 100)) + 
    labs(title = "Histogram of reaction time of mtrupdate",
         x = "Minutes", y = "Number of incidents")

## of mtrupdate (in chinese)
jpeg(filename = "reaction_time_chi.jpg", 1920, 1080, res = 200, quality = 100)
df4.match %>% 
    ggplot(aes(timediff.x)) +  
    geom_histogram(binwidth = 2, center = 2 / 2) + 
    coord_cartesian(xlim = c(-10, 100)) + 
    labs(x = "\u5206\u9418", y = "\u4E8B\u6545\u5B97\u6578")
dev.off()





### 06) Bar chart of daily number of incidents, for 2 types and 2 lines
df1.days %>% 
    ggplot(aes(signal)) + # flexible position of aes - see below
    geom_bar() +
    labs(title = "Bar chart of daily number of signal faults",
         x = "Number of incidents", y = "Count of days")
df1.days %>% 
    ggplot + 
    geom_bar(aes(train)) + # flexible position of aes - see above
    labs(title = "Bar chart of daily number of faulty trains",
         x = "Number of incidents", y = "Count of days")
df1.days %>% 
    ggplot(aes(east.rail)) +
    geom_bar() +
    labs(title = "Bar chart of daily number of incidents, East Rail Line",
         x = "Number of incidents", y = "Count of days")
df1.days %>% 
    ggplot + 
    geom_bar(aes(kwun.tong)) + 
    labs(title = "Bar chart of daily number of incidents, Kwun Tong Line",
         x = "Number of incidents", y = "Count of days")
## not really comparable b/w them (and with aes(all))





### 07) Bar chart of daily number of incidents, with fitted results
df1.days %>% 
    ggplot(aes(all)) + 
    geom_bar() + # counted when visualising
    labs(title = "Bar chart of daily number of all incidents",
         x = "Number of incidents", y = "Count of days")
## above and below produce the SAME plot!
## this is where objects obtained in modelling.R kick in
df1.counts %>% # counted when modelling
    ggplot(aes(count, all)) + 
    geom_bar(stat = "identity") + # this is the trick
    labs(title = "Bar chart of daily number of all incidents",
         x = "Number of incidents", y = "Count of days")
## after fitting, compare w/ above - same, but add estimated values
df1.counts %>% 
    tidyr::gather(type, value, all, all.est.rqk) %>% 
    ggplot(aes(count, value)) + 
    geom_bar(aes(fill = type), stat = "identity", position = "dodge") + 
    labs(title = "Bar chart of daily number of all incidents",
         x = "Number of incidents", y = "Count of days") + 
    scale_fill_discrete(name = NULL, 
                        labels = c("Observed", "Estimated"),
                        breaks = c("all", "all.est.rqk"))

## in chinese
jpeg(filename = "bar_observed_estimated_chi.jpg", 1920, 1080, res = 200, quality = 100)
df1.counts %>%
    filter(count != 0) %>% 
    tidyr::gather(type, value, all, all.est.rqk) %>% 
    ggplot(aes(count, value)) + 
    geom_bar(aes(fill = type), stat = "identity", position = "dodge") +
    labs(x = "\u6BCF\u65E5\u4E8B\u6545\u5B97\u6578", y = "\u65E5\u6578") +
    scale_fill_discrete(name = NULL, 
                        labels = c("\u6578\u64DA", "\u6A21\u578B\u4F30\u7B97"),
                        breaks = c("all", "all.est.rqk"))
dev.off()




### 08) type & line (vs time) - how


