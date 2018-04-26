
#' Read LXB files
#'
#' @param path path to lxb files. It can be a single file, a directory or a zipped directory
#' @param ... arguments that will be passed to \code{lxb::readLxb}.
#'
#' @return
#' A dataframe with all the data from the lxb files. A new column "Well"
#' is added to indicate the well of origin for every line.
#'
#' @details
#' In case a zip file is provided, it will be unziped into a \code{tmpdir} which
#' will be deleted \code{on.exit}.
#'
#' @author Teo Sakel
#' @export
#'
#' @examples
read_lxb <- function(path, ...) {

    file_ext <- tools::file_ext(path)
    if (file_ext == "lxb") {
        data <- as.data.frame(lxb::readLxb(path, ...))
        well <- stringr::str_extract(path, "[a-zA-Z][0-9]{1,2}.lxb")
        well <- sub("\\.lxb$", "", well)
        well <- data.frame(Well = rep(well, nrow(data)),
                           stringsAsFactors = FALSE)
        return(cbind(well, data))
    } else if (file_ext == "zip") {
        dpath <- tempdir()
        on.exit(unlink(dpath, recursive = TRUE))
        unzip(path, exdir = dpath)
        path <- dpath
    } else if (file_ext == "" && !dir.exists(path)) {
        # is should be already a dir
        stop(sprintf("%s has no file extension and is not a directory.", path))
    }

    .read_lxb_dir(path, ...)
}

.read_lxb_dir <- function(path, ...) {
    # TODO: Bad idea to use recursive = TRUE?
    data <- lxb::readLxb(list.files(path, pattern = "*.lxb", full.names = TRUE,
                                    recursive = TRUE), ...)
    wells <- data.frame(
        Well = rep(names(data), vapply(data, nrow, FUN.VALUE = 0L)),
        stringsAsFactors = FALSE
    )

    cbind(wells, as.data.frame(do.call(rbind, data)))
}