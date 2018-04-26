.splitAt <- function(x, pos)
    unname(split(x, cumsum(seq_along(x) %in% pos)))

.chr_vector_to_df <- function(Text, ...) {
    con <- textConnection(paste(Text, collapse = "\n"))
    tryCatch(
        read.table(con,
                   stringsAsFactors = FALSE,
                   na.strings = c("", "NA"),
                   check.names = FALSE,
                   fill = TRUE,
                   ...),
        finally = close(con)
    )
}
