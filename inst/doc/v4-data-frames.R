## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
## First define data.
(sacct.df <- data.frame(
  position=c(
    "chr10:213,054,000-213,055,000",
    "chrM:111,000-222,000",
    "this will not match",
    NA, # neither will this.
    "chr1:110-111 chr2:220-222"), # two possible matches.
  JobID=c(
    "13937810_25",
    "13937810_25.batch",
    "13937810_25.extern",
    "14022192_[1-3]",
    "14022204_[4]"),
  stringsAsFactors=FALSE))
remove.commas <- function(x)gsub(",", "", x)
long.list <- list()

## namedCapture: 29 lines of code.
range.list <- list(
  "\\[",
  task1="[0-9]+", as.integer,
  "(?:-",#begin optional end of range.
  taskN="[0-9]+", as.integer,
  ")?", #end is optional.
  "\\]")
task.list <- list(
  "(?:",#begin alternate
  task="[0-9]+", as.integer,
  "|",#either one task(above) or range(below)
  range.list,
  ")")#end alternate
to.int <- function(x)as.integer(remove.commas(x))
(long.list$namedCapture <- namedCapture::df_match_variable(
  sacct.df,
  JobID=list(
    job="[0-9]+", as.integer,
    "_",
    task.list,
    "(?:[.]",
    type=".*",
    ")?"),
  position=list(
    chrom="chr.*?",
    ":",
    chromStart=".*?", to.int,
    "-",
    chromEnd="[0-9,]*", to.int)))

## tidyr: 46 lines of code.
range.vec <- c(
  "\\[",
  task1="[0-9]+", 
  "(?:-",#begin optional end of range.
  taskN="[0-9]+", 
  ")?", #end is optional.
  "\\]")
task.vec <- c(
  "(?:",#begin alternate
  task="[0-9]+", 
  "|",#either one task(above) or range(below)
  range.vec,
  ")")#end alternate
regex.list <- list(
  JobID=c(
    job="[0-9]+", 
    "_",
    task.vec,
    "(?:[.]",
    type=".*",
    ")?"),
  position=c(
    chrom="chr.*?",
    ":",
    chromStart=".*?",
    "-",
    chromEnd="[0-9,]*"))
tidyr.input <- transform(
  sacct.df,
  position=remove.commas(position))
tidyr.df.list <- list(sacct.df)
for(col.name in names(regex.list)){
  regex.vec <- regex.list[[col.name]]
  is.group <- names(regex.vec)!=""
  format.vec <- ifelse(is.group, "(%s)", "%s")
  group.vec <- sprintf(format.vec, regex.vec)
  regex <- paste(group.vec, collapse="")
  group.names <- names(regex.vec)[is.group]
  result <- tidyr::extract(
    tidyr.input, col.name, group.names, regex, convert=TRUE)
  to.save <- result[, group.names, drop=FALSE]
  names(to.save) <- paste0(col.name, ".", group.names)
  tidyr.df.list[[col.name]] <- to.save
}
names(tidyr.df.list) <- NULL
long.list$tidyr <- do.call(cbind, tidyr.df.list)

## Make sure the results are the same.
t(sapply(long.list, names))
t(sapply(long.list, sapply, class))
long.list$tidyr$JobID.type <- ifelse(
  is.na(long.list$tidyr$JobID.type),
  "",
  long.list$tidyr$JobID.type)
with(long.list, identical(tidyr, namedCapture))


## -----------------------------------------------------------------------------
## First define data.
(sacct.df <- data.frame(
  position=c(
    "chr10:213,054,000-213,055,000",
    "chrM:111,000-222,000",
    "this will not match",
    NA, # neither will this.
    "chr1:110-111 chr2:220-222"), # two possible matches.
  JobID=c(
    "13937810_25",
    "13937810_25.batch",
    "13937810_25.extern",
    "14022192_[1-3]",
    "14022204_[4]"),
  stringsAsFactors=FALSE))
short.list <- list()

## tidyr alternate (13 lines total)
e <- function(col.name, group.names, pattern){
  result <- tidyr::extract(
    sacct.df, col.name, group.names, pattern, convert=TRUE)
  to.save <- result[, group.names, drop=FALSE]
  names(to.save) <- paste0(col.name, ".", group.names)
  to.save
}
short.list$tidyr <- do.call(cbind, list(
  sacct.df,
  e("JobID", c("job", "task", "task1", "taskN", "type"),
    "([0-9]+)_(?:([0-9]+)|\\[([0-9]+)(?:-([0-9]+))?\\])(?:[.](.*))?"),
  e("position", c("chrom", "chromStart", "chromEnd"),
    "(chr.*?):(.*?)-([0-9,]*)")))

## namedCapture alternate (7 lines total)
(short.list$namedCapture <- namedCapture::df_match_variable(
  sacct.df,
  JobID="(?P<job>[0-9]+)_(?:(?P<task>[0-9]+)|\\[(?P<task1>[0-9]+)(?:-(?P<taskN>[0-9]+))?\\])(?:[.](?P<type>.*))?",
  position="(?P<chrom>chr.*?):(?P<chromStart>.*?)-(?P<chromEnd>[0-9,]*)"))
for(N in names(short.list$namedCapture)){
  short.list$namedCapture[[N]] <- type.convert(short.list$namedCapture[[N]], as.is=TRUE)
}

## Make sure the results are the same.
t(sapply(short.list, names))
t(sapply(short.list, sapply, class))
short.list$tidyr$JobID.type <- ifelse(
  is.na(short.list$tidyr$JobID.type),
  "",
  short.list$tidyr$JobID.type)
with(short.list, identical(tidyr, namedCapture))


## -----------------------------------------------------------------------------
range.list <- list(
  "\\[",
  task1="[0-9]+", as.integer,
  list(
    "-",#begin optional end of range.
    taskN="[0-9]+", as.integer
  ), "?", #end is optional.
  "\\]")
namedCapture::df_match_variable(sacct.df, JobID=range.list)

range.pat <- paste0(
  "\\[",
  "(?<task1>[0-9]+)", 
  "(?:",
  "-",#begin optional end of range.
  "(?<taskN>[0-9]+)",
  ")?", #end is optional.
  "\\]")
rematch2::bind_re_match(sacct.df, JobID, range.pat)
task.list <- list(
  "_",
  list(
    task="[0-9]+", as.integer,
    "|",#either one task(above) or range(below)
    range.list))
namedCapture::df_match_variable(sacct.df, JobID=task.list)

task.pat <- paste0(
  "_",
  "(?:",
  "(?<task>[0-9]+)", 
  "|", #either one task(above) or range(below)
  range.pat,
  ")")
rematch2::bind_re_match(sacct.df, JobID, task.pat)

job.list <- list(
  job="[0-9]+", as.integer,
  task.list,
  list(
    "[.]",
    type=".*"
  ), "?")
(job.namedCapture <- namedCapture::df_match_variable(sacct.df, JobID=job.list))

job.pat <- paste0(
  "(?<job>[0-9]+)", 
  task.pat,
  "(?:",
  "[.]",
  "(?<type>.*)",
  ")?")
(job.rematch2 <- rematch2::bind_re_match(sacct.df, JobID, job.pat))

pos.namedCapture <- namedCapture::df_match_variable(
  job.namedCapture, position=list(
    chrom="chr.*?",
    ":",
    chromStart=".*?", to.int,
    "-",
    chromEnd="[0-9,]*", to.int))
str(pos.namedCapture)

pos.rematch2 <- rematch2::bind_re_match(
  job.rematch2,  position, paste0(
    "(?<chrom>chr.*?)",
    ":",
    "(?<chromStart>.*?)", 
    "-",
    "(?<chromEnd>[0-9,]*)"))
str(pos.rematch2)


## -----------------------------------------------------------------------------
converted.rematch2 <- transform(
  pos.rematch2,
  JobID.job=to.int(job),
  JobID.task1=to.int(task1),
  JobID.taskN=to.int(taskN),
  JobID.task=to.int(task),
  JobID.type=type,
  position.chrom=chrom,
  position.chromStart=to.int(chromStart),
  position.chromEnd=to.int(chromEnd),
  stringsAsFactors=FALSE)
some.rematch2 <- converted.rematch2[, names(pos.namedCapture)]
identical(some.rematch2, pos.namedCapture)

