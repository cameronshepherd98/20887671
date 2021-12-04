# Purpose

Purpose of this work folder.

Ideally store a minimum working example data set in data folder.

Add binary files in bin, and closed R functions in code. Human Readable
settings files (e.g. csv) should be placed in settings/

``` r
rm(list = ls()) # Clean your environment:
gc() # garbage collection - It can be useful to call gc after a large object has been removed, as this may prompt R to return memory to the operating system.
```

    ##          used (Mb) gc trigger (Mb) limit (Mb) max used (Mb)
    ## Ncells 408172 21.8     837302 44.8         NA   657805 35.2
    ## Vcells 767369  5.9    8388608 64.0      16384  1802254 13.8

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.4     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.3     ✓ stringr 1.4.0
    ## ✓ readr   2.0.1     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
list.files('code/', full.names = T, recursive = T) %>% .[grepl('.R', .)] %>% as.list() %>% walk(~source(.))
```

``` r
#Texevier::create_template_html(directory = "Questions", template_name = "Question1")

#Texevier::create_template_html(directory = "Questions", template_name = "Question2")

#Texevier::create_template_html(directory = "Questions", template_name = "Question3")

#Texevier::create_template_html(directory = "Questions", template_name = "Question4")

#Texevier::create_template_html(directory = "Questions", template_name = "Question5")

#Texevier::create_template_html(directory = "Questions", template_name = "Question6")
```

\#Question 1
