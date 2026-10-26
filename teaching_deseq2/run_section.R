# =============================================================================
#  run_section.R: run DESeq2_TNBC_vs_nonTNBC.R one section at a time in a
#  plain terminal R session (no RStudio), e.g. on a cluster compute node.
#
#    source("run_section.R")   # once, after starting R in this folder
#    section(2)                # runs section "2. Read the data" and shows
#                              # every line with its output, like Ctrl+Enter
#    sections()                # lists the section numbers and titles
#
#  Everything runs in your normal workspace, so objects made in one section
#  (counts, dds, res, ...) are there for the next one and for your own typing.
# =============================================================================

.script <- "DESeq2_TNBC_vs_nonTNBC.R"

# A section header is "# N. Title" with a "# ====" line directly above and below.
.section_starts <- function() {
    x <- readLines(.script)
    bar <- grepl("^# =+$", x)
    i <- which(grepl("^# [0-9]+\\. ", x) & c(FALSE, head(bar, -1)) & c(tail(bar, -1), FALSE))
    data.frame(n = as.integer(sub("^# ([0-9]+)\\..*", "\\1", x[i])),
               title = sub("^# ", "", x[i]), from = i - 1,
               to = c(i[-1] - 2, length(x)))
}

sections <- function() {
    s <- .section_starts()
    writeLines(s$title)
    invisible(s)
}

section <- function(n) {
    s <- .section_starts()
    if (!n %in% s$n) stop("No section ", n, ". Type sections() to see them.")
    lines <- readLines(.script)[s$from[s$n == n]:s$to[s$n == n]]
    tmp <- tempfile(fileext = ".R")
    writeLines(lines, tmp)
    # keep.source = TRUE prints the script's comments too (they are the lesson)
    source(tmp, echo = TRUE, max.deparse.length = Inf, spaced = TRUE, local = globalenv(),
           keep.source = TRUE)
    invisible(NULL)
}

message("Ready. Type sections() for the list, then section(1), section(2), ...")
