# Analysis of MTR incidents

Statistical analysis of MTR incidents in Hong Kong is carried out on data from tweets by MTR service update, using the following R script files:

1. **extraction.R** cleans the raw data extracted (from a .csv file) to obtain objects usable for analysis and visualisation. 
2. The main object is a data table with each row corresponding to one incident, with info such as:
 * operational date, line, type
 * start tweet: content, time of occurrence, time of tweet creation
 * end tweet: content, time of occurrence, time of tweet creation
3. **modelling.R** contains the codes for fitting the Poisson & negative Binomial models, with possibly linear trend and/or changepoint, to the weekly number of incidents.
4. **visualisation.R** visualises different aspects of the data, esp. of the main data table, via R package ggplot2.
5. Flow: **visualisation.R** sources **modelling.R**, which in turn sources **extraction.R**, which in turn requires **stations_and_lines.csv** & **functions_mtr.R**.

For the latest news by MTR service update, which is both the data collector and provider, you can follow their [Twitter](https://twitter.com/mtrupdate "@mtrupdate"). For other information, please visit their [Facebook](https://www.facebook.com/mtrupdate/) or [swiftzer homepage](http://checkfare.swiftzer.net/).