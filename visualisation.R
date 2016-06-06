### Principle: we shouldn't do any further data manipulation like counting here. 
### Any summaries can and should be done via stat functions in ggplot2.

remove(list = ls())
source("modelling.R") # run the whole script, which in turn runs extraction.R
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

## bar chart by line
df4.match %>% 
    ggplot() + 
    geom_bar(aes(line)) + 
    labs(title = "Number of incidents by line",
         y = "Number of incidents") +
    coord_flip()





### 02) Overall temporal analysis
## basic with graphics::hist
df4.match %>% 
    extract2("month.mid") %>% 
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

## ACF (to show independence)
df0.months %>% 
    magrittr::extract2("all") %>% 
    acf(50) # slight yearly effect
df0.weeks %>% 
    magrittr::extract2("all") %>%
    acf(100)





### 03) Type vs time: monthly number of incidents by type
## raster/tile/rectangle
df4.match %>% 
    count(month.mid, type) %>%
    ggplot(aes(month.mid, type, fill = n)) + 
    geom_tile() + # geom_raster() fails big here
    scale_fill_gradientn(name = NULL,
                         colours = topo.colors(12) %>% rev,
                         breaks = c(1, 4, 7, 10, 13)) + 
    labs(title = "Monthly number of incidents by type", 
         x = "Time", y = "Type") +
    coord_cartesian(xlim = c(as.Date("2011-03-01"), as.Date("2016-02-01")))

## stacked histogram over time
df4.match %>% 
    ggplot(aes(month.mid)) + 
    geom_histogram(aes(fill = type), bins = 63) + # geom_bar() fails big here
    labs(x = "Time", y = "Count") # not very easily interpretable






### 04) Line vs time: monthly number of incidents by line
## raster/tile/rectangle
df4.match %>% 
    count(month.mid, line) %>%
    ggplot(aes(month.mid, line)) + # flexible location of fill
    geom_tile(aes(fill = n)) +
    scale_fill_gradientn(name = NULL,
                         colours = topo.colors(12) %>% rev) + 
    labs(title = "Monthly number of incidents by line",
         x = "Time", y = "Line") +
    coord_cartesian(xlim = c(as.Date("2011-03-01"), as.Date("2016-02-01")))

## stacked histogram over time
df4.match %>% 
    ggplot(aes(month.mid)) + 
    geom_histogram(aes(fill = line), bins = 63) +
    labs(x = "Time", y = "Count") # not very easily interpretable




### 05) efficiency
## of MTR
df4.match %>% 
    filter(downtime != 0 & !is.na(downtime)) %>%
    ggplot(aes(downtime)) + 
    geom_histogram(binwidth = 10, center = 10 / 2) + 
    coord_cartesian(xlim = c(-10, 800)) + 
    labs(title = "Histogram of incident downtime",
         x = "Minutes", y = "Number of incidents")

## of mtrupdate
df4.match %>% 
    ggplot(aes(timediff.x)) +  
    geom_histogram(binwidth = 2, center = 2 / 2) + 
    coord_cartesian(xlim = c(-10, 100)) + 
    labs(title = "Histogram of reaction time of mtrupdate",
         x = "Minutes", y = "Number of incidents")





### 06) Bar chart of weekly number of incidents, for 2 types and 2 lines
df0.weeks %>% 
    ggplot(aes(signal)) + # flexible position of aes - see below
    geom_bar() +
    labs(title = "Bar chart of weekly number of signal faults",
         x = "Number of incidents", y = "Count of weeks")
df0.weeks %>% 
    ggplot + 
    geom_bar(aes(train)) + # flexible position of aes - see above
    labs(title = "Bar chart of weekly number of faulty trains",
         x = "Number of incidents", y = "Count of weeks")
df0.weeks %>% 
    ggplot(aes(east.rail)) +
    geom_bar() +
    labs(title = "Bar chart of weekly number of incidents, East Rail Line",
         x = "Number of incidents", y = "Count of weeks")
df0.weeks %>% 
    ggplot + 
    geom_bar(aes(kwun.tong)) + 
    labs(title = "Bar chart of weekly number of incidents, Kwun Tong Line",
         x = "Number of incidents", y = "Count of weeks")
## not really comparable b/w them (and with aes(all))





### 07) Bar chart of weekly number of incidents, with fitted results
df0.weeks %>% 
    ggplot(aes(all)) + 
    geom_bar() + # counted when visualising
    labs(title = "Bar chart of weekly number of all incidents",
         x = "Number of incidents", y = "Count of weeks")
## above and below produce the SAME plot!
## this is where objects obtained in modelling.R kick in
df0.counts %>% # counted when extracting
    ggplot(aes(count, all)) + 
    geom_bar(stat = "identity") + # this is the trick
    labs(title = "Bar chart of weekly number of all incidents",
         x = "Number of incidents", y = "Count of weeks")
## after fitting, compare w/ above - same, but add estimated values
df0.counts %>% 
    gather(type, value, all, all.est.nb) %>% 
    ggplot(aes(count, value)) + 
    geom_bar(aes(fill = type), stat = "identity", position = "dodge") + 
    labs(title = "Bar chart of weekly number of all incidents",
         x = "Number of incidents", y = "Count of weeks") + 
    scale_fill_discrete(name = NULL, 
                        labels = c("Observed", "Estimated"),
                        breaks = c("all", "all.est.nb"))





### 08) type & line (vs time) - how


