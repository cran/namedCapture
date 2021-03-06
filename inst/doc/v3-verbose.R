## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
subject <- '198.214.42.14 - - [21/Jul/1995:14:31:46 -0400] "GET /images/ HTTP/1.0" 200 17688
lahal.ksc.nasa.gov - - [24/Jul/1995:12:42:40 -0400] "GET /images/USA-logosmall.gif HTTP/1.0" 200 234
199.171.112.23 - - [02/Jul/1995:02:30:34 -0400] "GET /images/KSC-logosmall.gif HTTP/1.0" 200 1204
gate3.fmr.com - - [05/Jul/1995:13:51:39 -0400] "GET /shuttle/countdown/ HTTP/1.0" 200 3985
curly02.slip.yorku.ca - - [10/Jul/1995:23:11:49 -0400] "GET /shuttle/missions/sts-70/sts-70-patch-small.gif HTTP/1.0" 200 5026
boson.epita.fr - - [15/Jul/1995:11:27:49 -0400] "GET /shuttle/missions/sts-71/movies/sts-71-mir-dock.mpg HTTP/1.0" 200 946425
134.153.50.9 - - [13/Jul/1995:11:02:50 -0400] "GET /icons/text.xbm HTTP/1.0" 200 527
port00.ventura.rain.org - - [23/Jul/1995:09:11:06 -0400] "GET /shuttle/countdown/ HTTP/1.0" 200 4324
128.159.145.91 - - [14/Jul/1995:10:38:04 -0400] "GET /statistics/images/getstats_big.gif HTTP/1.0" 200 6777
slo.eei.upmc.edu - - [25/Jul/1995:09:33:01 -0400] "GET /images/KSC-logosmall.gif HTTP/1.0" 200 1204
206.13.med.umich.edu - - [14/Jul/1995:09:11:28 -0400] "GET /shuttle/resources/orbiters/challenger-logo.gif HTTP/1.0" 200 4179'
subject.vec <- strsplit(subject, split="\n")[[1]]

result.list <- list()

## namedCapture 10 lines, each group name specified only once.
result.list$namedCapture <- namedCapture::str_match_variable(
  subject.vec,
  "\\[",
  time="[^]]*", function(x)as.POSIXct(x, format="%d/%b/%Y:%H:%M:%S %z"),
  "\\]",
  ' "GET ',
  list(
    "[^ ]+[.]",
    filetype='[^ .?"]+', tolower
  ), "?")

## free-spacing 14 lines, 2 repetitions of each group name.
options(namedCapture.engine="PCRE")
result.list$freespace <- namedCapture::str_match_named(
  subject.vec, '(?x) #activates PCRE_EXTENDED mode.
\\[
(?P<time>[^]]*) 
\\]
[ ]"GET[ ]
(?:
  [^ ]+[.]
  (?P<filetype>[^ .?"]+)
)?
', list(
  time=function(x)as.POSIXct(x, format="%d/%b/%Y:%H:%M:%S %z"),
  filetype=tolower))

## rex 17 lines, 3 repetitions of each group name.
library(rex)
library(dplyr)
result.list$rex <- re_matches(
  subject.vec,
  rex(
    "[",
    capture(name = "time",
            none_of("]") %>% zero_or_more()),
    "]",
    space, double_quote, "GET", space,
    maybe(
      non_spaces, ".",
      capture(name = 'filetype',
              none_of(space, ".", "?", double_quote) %>% one_or_more())
    )
  )
) %>%
  mutate(filetype = tolower(filetype),
         time = as.POSIXct(time, format="%d/%b/%Y:%H:%M:%S %z"))

with(result.list, identical(rex, namedCapture))
with(result.list, identical(freespace, namedCapture))


