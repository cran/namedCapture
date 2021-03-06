library(testthat)
library(namedCapture)
context("variable args syntax")
source(system.file("test_engines.R", package="namedCapture", mustWork=TRUE), local=TRUE)

subject <- c(
  ten="chr10:213,054,000-213,055,000",
  chrNA="chrNA:111,000-222,000",
  no.match="foo bar",
  missing=NA,
  two="chr1:110-111 chr2:220-222")

test_engines("str_match_variable returns character matrix", {
  computed <- str_match_variable(
    subject,
    chrom="chr.*?",
    ":",
    chromStart=".*?",
    "-",
    chromEnd="[0-9,]*")
  expected <- cbind(
    chrom=c("chr10", "chrNA", NA, NA, "chr1"),
    chromStart=c("213,054,000", "111,000", NA, NA, "110"),
    chromEnd=c("213,055,000", "222,000", NA, NA, "111"))
  rownames(expected) <- names(subject)
  expect_identical(computed, expected)
})

keep.digits <- function(x)as.integer(gsub("[^0-9]", "", x))
test_engines("str_match_variable returns data.frame", {
  computed <- str_match_variable(
    subject,
    chrom="chr.*?",
    ":",
    chromStart=".*?", keep.digits,
    "-",
    chromEnd="[0-9,]*", keep.digits)
  expected <- data.frame(
    chrom=c("chr10", "chrNA", NA, NA, "chr1"),
    chromStart=as.integer(c(213054000, 111000, NA, NA, 110)),
    chromEnd=as.integer(c(213055000, 222000, NA, NA, 111)),
    stringsAsFactors=FALSE)
  rownames(expected) <- names(subject)
  expect_equivalent(computed, expected)
})

test_engines("named function is an error", {
  expect_error({
    str_match_variable(
      subject,
      chrom="chr.*?",
      ":",
      chromStart=".*?", fun=keep.digits,
      "-",
      chromEnd="[0-9,]*", keep.digits)
  }, "functions must not be named, problem: fun")
})

test_engines("str_match_all_variable returns character matrix", {
  computed <- str_match_all_variable(
    subject,
    chrom="chr.*?",
    ":",
    chromStart=".*?",
    "-",
    chromEnd="[0-9,]*")
  r <- function(chrom, chromStart, chromEnd){
    cbind(chrom=chrom, chromStart=chromStart, chromEnd=chromEnd)
  }
  expected <- rbind(
    r("chr10", "213,054,000", "213,055,000"),
    r("chrNA", "111,000", "222,000"),
    r("chr1", "110", "111"),
    r("chr2", "220", "222"))
  expect_identical(computed, expected)
})

test_engines("str_match_all_variable removes missing subjects", {
  computed <- str_match_all_variable(
    subject,
    "(?P<na>NA)")
  ## There should be only one NA (not two) because chrNA matches but
  ## the missing NA subject should be removed.
  expect_identical(computed, cbind(na="NA"))
})

test_engines("str_match_all_variable returns data.frame", {
  conversion.list <- list(chromStart=keep.digits, chromEnd=keep.digits)
  computed <- str_match_all_variable(
    subject,
    chrom="chr.*?",
    ":",
    chromStart=".*?", keep.digits,
    "-",
    chromEnd="[0-9,]*", keep.digits)
  expected <- rbind(
    data.frame(
      chrom="chr10", chromStart=213054000L, chromEnd=213055000L,
      stringsAsFactors=FALSE),
    data.frame(
      chrom="chrNA", chromStart=111000L, chromEnd=222000L,
      stringsAsFactors=FALSE),
    data.frame(
      chrom=c("chr1", "chr2"),
      chromStart=as.integer(c("110", "220")),
      chromEnd=as.integer(c("111", "222")),
      stringsAsFactors=FALSE))
  expect_identical(computed, expected)
})

test_engines("str_match_variable errors for one argument", {
  expect_error({
    str_match_variable("foo")
  }, "pattern must have at least one argument")
})

test_engines("str_match_all_variable errors for one argument", {
  expect_error({
    str_match_all_variable("foo")
  }, "pattern must have at least one argument")
})

test_engines("str_match_variable errors for multi-dim patterns", {
  expect_error({
    str_match_variable("foo", c("bar", "baz"))
  }, "patterns must be character vectors of length 1")
})

test_engines("str_match_all_variable errors for multi-dim patterns", {
  expect_error({
    str_match_all_variable("foo", c("bar", "baz"))
  }, "patterns must be character vectors of length 1")
})

test_engines("str_match_variable errors for 0-length patterns", {
  expect_error({
    str_match_variable("foo", character())
  }, "patterns must be character vectors of length 1")
})

test_engines("str_match_all_variable errors for 0-length patterns", {
  expect_error({
    str_match_all_variable("foo", character())
  }, "patterns must be character vectors of length 1")
})

test_engines("str_match_variable errors for non char/fun args", {
  expect_error({
    str_match_variable("foo", "bar", 1)
  }, "arguments must be", fixed=TRUE)
})

test_engines("str_match_all_variable errors for non char/fun args", {
  expect_error({
    str_match_all_variable("foo", "bar", 1)
  }, "arguments must be", fixed=TRUE)
})

test_engines("str_match_variable errors for two funs in a row", {
  expect_error({
    str_match_variable("foo", g="bar", as.integer, as.numeric)
  },
  "too many functions; up to one function may follow each named pattern")
})

test_engines("str_match_all_variable errors for two funs in a row", {
  expect_error({
    str_match_all_variable("foo", g="bar", as.integer, as.numeric)
  },
  "too many functions; up to one function may follow each named pattern")
})

test_engines("str_match_variable errors for fun at start", {
  expect_error({
    str_match_variable("foo", as.numeric)
  },
  "too many functions; up to one function may follow each named pattern")
})

test_engines("str_match_all_variable errors for fun at start", {
  expect_error({
    str_match_all_variable("foo", as.numeric)
  },
  "too many functions; up to one function may follow each named pattern")
})

test_engines("str_match_variable errors for NA pattern", {
  expect_error({
    str_match_variable("foo", g="bar", NA_character_, "baz")
  }, "patterns must not be missing/NA")
})

test_engines("str_match_all_variable errors for NA pattern", {
  expect_error({
    str_match_all_variable("foo", g="bar", NA_character_, "baz")
  }, "patterns must not be missing/NA")
})

range.pattern <- list(
  "[[]",
  task1="[0-9]+", as.integer,
  "(?:-",#begin optional end of range.
  taskN="[0-9]+", as.integer,
  ")?", #end is optional.
  "[]]")
full.pattern <- list(
  job="[0-9]+", as.integer,
  "_",
  "(?:",#begin alternate
  task="[0-9]+", as.integer,
  "|",#either one task(above) or range(below)
  range.pattern,
  ")",#end alternate
  "(?:[.]",
  type=".*",
  ")?")
subject.vec <- c(
  "13937810_25",
  "13937810_25.batch",
  "13937810_25.extern",
  "14022192_[1-3]",
  "14022204_[4]")
all.args <- list(subject.vec, full.pattern)
test_engines("nested lists are OK", {
  task.df <- do.call(str_match_variable, all.args)
  expect_identical(
    names(task.df),
    c("job", "task", "task1", "taskN", "type"))
  expect_identical(task.df$job, as.integer(c(
    13937810, 13937810, 13937810, 14022192, 14022204)))
  expect_identical(task.df$task, as.integer(c(
    25, 25, 25, NA, NA)))
  expect_identical(task.df$task1, as.integer(c(
    NA, NA, NA, 1, 4)))
  expect_identical(task.df$taskN, as.integer(c(
    NA, NA, NA, 3, NA)))
  expect_identical(task.df$type, c(
    "", "batch", "extern", "", ""))
})


trackDb.txt.gz <- system.file("extdata", "trackDb.txt.gz", package="namedCapture")
trackDb.vec <- readLines(trackDb.txt.gz)

test_engines("nested capture groups works", {
  name.pattern <- list(
    cellType=".*?",
    "_",
    sampleName=list(as.factor,
                    "McGill",
                    sampleID="[0-9]+", as.integer),
    dataType="Coverage|Peaks",
    "|",
    "[^\n]+")
  match.df <- namedCapture::str_match_all_variable(
    trackDb.vec,
    "track ",
    name=name.pattern,
    "(?:\n[^\n]+)*",
    "\\s+bigDataUrl ",
    bigDataUrl="[^\n]+")
  expect_is(match.df, "data.frame")
  expect_identical(
    names(match.df),
    c("cellType", "sampleName", "sampleID", "dataType", "bigDataUrl"))
  expect_is(match.df$sampleName, "factor")
  expect_is(match.df$sampleID, "integer")
})

subject.vec <- c(
  "chr10:213,054,000-213,055,000",
  "chrM:111,000",
  "this will not match",
  NA, # neither will this.
  "chr1:110-111 chr2:220-222") # two possible matches.
chr.pos.df <- str_match_variable(
  subject.vec,
  chrom="chr.*?",
  ":",
  chromStart="[0-9,]+", keep.digits,
  list(
    "-",
    chromEnd="[0-9,]+", keep.digits
  ), "?")
test_engines("un-named list interpreted as non-capturing group", {
  expect_identical(
    chr.pos.df$chromStart,
    as.integer(c(213054000, 111000, NA, NA, 110)))
  expect_identical(
    chr.pos.df$chromEnd,
    as.integer(c(213055000, NA, NA, NA, 111)))
})

matching.subjects <- c(
  "chr10:213,054,000-213,055,000",
  "chrM:111,000",
  "chr1:110-111 chr2:220-222") # two possible matches.
test_engines("str subject no error if nomatch.error=TRUE and all matches", {
  match.df <- str_match_variable(
    matching.subjects, nomatch.error=TRUE,
    chrom="chr.*?",
    ":",
    chromStart="[0-9,]+", keep.digits,
    list(
      "-",
      chromEnd="[0-9,]+", keep.digits
    ), "?")
  expect_identical(
    match.df$chromEnd,
    as.integer(c(213055000, NA, 111)))
})
test_engines("str subject stop if nomatch.error=TRUE and no match", {
  expect_error({
    str_match_variable(
      subject.vec, nomatch.error=TRUE,
      chrom="chr.*?",
      ":",
      chromStart="[0-9,]+", keep.digits,
      list(
        "-",
        chromEnd="[0-9,]+", keep.digits
      ), "?")
  }, "subjects printed above did not match regex below")
})

test_engines("df subject no error if nomatch.error=TRUE and all matches", {
  subject.df <- data.frame(subject.col=matching.subjects, stringsAsFactors=FALSE)
  match.df <- df_match_variable(
    subject.df,
    subject.col=list(
      nomatch.error=TRUE,
      chrom="chr.*?",
      ":",
      chromStart="[0-9,]+", keep.digits,
      list(
        "-",
        chromEnd="[0-9,]+", keep.digits
      ), "?"))
  expect_identical(
    match.df$subject.col.chromEnd,
    as.integer(c(213055000, NA, 111)))
})
test_engines("df subject stop if nomatch.error=TRUE and no match", {
  subject.df <- data.frame(subject.vec, stringsAsFactors=FALSE)
  expect_error({
    df_match_variable(
      subject.df,
      subject.vec=list(
        nomatch.error=TRUE,
        chrom="chr.*?",
        ":",
        chromStart="[0-9,]+", keep.digits,
        list(
          "-",
          chromEnd="[0-9,]+", keep.digits
        ), "?"))
  }, "subjects printed above did not match regex below")
})

(foo.mat <- str_match_variable(
  c("foo", "foobar", "fooba"),
  first="foo",
  list("b", second="ar"), "?"))
test_engines("un-named list interpreted as non-capturing group foo subject", {
  expect_identical(foo.mat[, "first"], c("foo", "foo", "foo"))
  expect_identical(foo.mat[, "second"], c("", "ar", ""))
})

subject <- "foo55bar"
test_engines("str_match_variable returns mat with only one group = name", {
  out.mat <- namedCapture::str_match_variable(
    subject,
    name="[0-9]+")
  exp.mat <- matrix(character(), 1, 0, dimnames=list("55"))
  expect_identical(out.mat, exp.mat)
})
test_engines("str_match_variable returns df with only one group = name", {
  out.df <- namedCapture::str_match_variable(
    subject,
    name="[0-9]+", as.integer)
  exp.df <- data.frame(row.names="55")
  expect_identical(out.df, exp.df)
})


