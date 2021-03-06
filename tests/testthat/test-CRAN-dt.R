library(testthat)
library(namedCapture)
library(data.table)
context("dt")
source(system.file("test_engines.R", package="namedCapture", mustWork=TRUE), local=TRUE)

subject.dt <- data.table(
  JobID=c(
    "13937810_25",
    "13937810_25.batch",
    "13937810_25.extern",
    "14022192_[1-3]",
    "14022204_[4]"),
  position=c(
    "chr10:213,054,000-213,055,000",
    "chrNA:111,000-222,000",
    "foo bar",
    NA,
    "chr1:110-111 chr2:220-222"),
  stringsAsFactors=FALSE)
range.pattern <- list(
  "[[]",
  task1="[0-9]+", as.integer,
  "(?:-",#begin optional end of range.
  taskN="[0-9]+", as.integer,
  ")?", #end is optional.
  "[]]")
test_engines("df_match_variable returns data.table", {
  match.dt <- df_match_variable(
    subject.dt,
    JobID=list(
      job="[0-9]+", as.integer,
      "_",
      "(?:",#begin alternate
      task="[0-9]+", as.integer,
      "|",#either one task(above) or range(below)
      range.pattern,
      ")",#end alternate
      "(?:[.]",
      type=".*",
      ")?"),
    position=list(
      chrom="chr.*?",
      ":",
      chromStart=".*?",
      "-",
      chromEnd="[0-9,]*"))
  expect_identical(names(match.dt), c(
    "JobID", "position",
    "JobID.job", "JobID.task", "JobID.task1", "JobID.taskN", "JobID.type",
    "position.chrom", "position.chromStart", "position.chromEnd"))
  expect_identical(match.dt$JobID.job, as.integer(c(
    13937810, 13937810, 13937810, 14022192, 14022204)))
  expect_identical(match.dt$JobID.task, as.integer(c(
    25, 25, 25, NA, NA)))
  expect_identical(match.dt$JobID.task1, as.integer(c(
    NA, NA, NA, 1, 4)))
  expect_identical(match.dt$JobID.taskN, as.integer(c(
    NA, NA, NA, 3, NA)))
  expect_identical(match.dt$JobID.type, c(
    "", "batch", "extern", "", ""))
  expect_identical(match.dt$position.chrom, c(
    "chr10", "chrNA", NA, NA, "chr1"))
  expect_identical(match.dt$position.chromStart, c(
    "213,054,000", "111,000", NA, NA, "110"))
  expect_identical(match.dt$position.chromEnd, c(
    "213,055,000", "222,000", NA, NA, "111"))
  expect_is(match.dt, "data.table")
})

