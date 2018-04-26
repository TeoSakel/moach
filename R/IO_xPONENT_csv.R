#' Read a Luminex-xPONENT generated csv file
#'
#' @param path The path to a Luminex xponent csv file
#' @param sep field delimiter used
#'
#' @return A list of lists with all the information from the file.
#'     At the top level there are two lists
#'     \itemize{
#'         \item \code{BatchHeader} that stores all the information before the
#'         \emph{Results} section as well as the CRC and
#'         \item \code{AssayData} that stores all the information in the
#'         \emph{Results} section
#'     }
#'
#'     \code{BatchHeader} stores data as key-value pairs with keys taken from
#'     the first column of the csv file and values from the rest. Most entries
#'     are singled valued but some (like \code{CALInfo}) are nested lists. Check
#'     references for more details on csv format.
#'
#'     \code{AssayData} stores the \emph{Results} as dataframes, except for
#'     \code{Run} data which is a list of dataframes (Audit Logs and Warning/Errors).
#'     There are 4 dataframes in \code{AssayData}:
#'     \itemize{
#'         \item \strong{Exprs} Expression data, ie \emph{raw} data, including \strong{Medians, Counts, Net MFI}
#'         \item \strong{Average} Averages for replicate wells
#'         \item \strong{Well} Information about each well (currently only the \emph{Dilution Factor})
#'         \item \strong{Bead} Information about the beads, eg BeadID
#'     }
#'     again for more information about the generated data check the references.
#'
#' @details The format details (eg comma or semicolon delimiter) are infered by the function.
#'     There is a known edge-case, for some semicolon-separated files,
#'     where some calibration information is misplaced. This is the **only** case where you should
#'     mess with the csv file, otherwise, you should **not mess** with it.
#'
#' @references \itemize{
#'     \item http://www.appliedcytometry.com/Technotes/System_Operation/Technote_32/Luminex_100_IS_User_Manual_(IS2.3).pdf
#'     \item https://www.luminexcorp.com/blog/its-all-about-the-stats/
#' }
#'
#' @author Teo Sakel
#' @export
#'
#' @examples
read_xponent_csv <- function(path, sep = c("infer", ",", ";")) {
    #
    sep <- match.arg(sep)
    # Read and Clean Lines
    Lines <- readLines(path)
    if (sep == "infer")
        sep <- ifelse(stringr::str_count(Lines[1], ",") > 0, ",", ";")
    csv_regex <- paste0(sep, "(?=([^\"]*\"[^\"]*\")*[^\"]*$)")
    Lines <- stringr::str_replace(Lines, paste0(sep, "*$"), "") # remove trailing delimiters
    Lines <- stringr::str_split(Lines, csv_regex)
    Lines <- lapply(Lines, stringr::str_replace_all, "^\"|\"$", "") # unquote
    Lines <- vapply(Lines, paste, "", collapse = "\t")  # rejoin with tabs. A little stupid to be honest...

    # Organize Lines into Blocks
    empty_lines <- which(Lines == "") + 1L
    Lines <- .splitAt(Lines, empty_lines)
    Lines[vapply(Lines, function(x) all(x == ""), FALSE)] <- NULL  # remove empty blocks
    header_end <- grep("^Results", vapply(Lines, "[", "", 1L))
    header_lines <- Lines[1:(header_end - 1L)]
    results_lines <- Lines[(header_end + 1L):length(Lines)]

    Xcsv <- list(
        BatchHeader = .parse_Xcsv_header(header_lines),
        AssayData = .parse_Xcsv_results(results_lines)
    )
    # special case
    Xcsv[["BatchHeader"]][["CRC"]] <- tryCatch(
        as.character(results_lines[[length(results_lines)]][2]),
        error = function(e) NA_character_
    )

    class(Xcsv) <- c("xponentCSV", class(Xcsv))
    return(Xcsv)
}


# SubMethods ----------------------------------------------------------------------------------

.parse_Xcsv_header <- function(Lines) {
    fields <- vapply(Lines, "[", "", 1)

    # Parse Batch Info - first two fields are Batch Info
    batch_info <- .parse_Xcsv_BatchInfo(c(Lines[[1]], Lines[[2]]))

    # Parse Calibration Info
    lastcal <- grep("^Most Recent Calibration and Verification Results:", fields)
    calinfo <- grep("^CALInfo:", fields)
    CALInfo <- .parse_Xcsv_CALInfo(list(unlist(Lines[lastcal]), unlist(Lines[calinfo])))

    # Parse Sample Field
    # 2 fields: Samples: <NumSamples> MinEvents: <MinEvents> "Per X"
    sample_field <- Lines[[length(Lines)]]
    sample_field <- stringr::str_split(sample_field[grep("^Samples", sample_field)], "\t")
    sample_field <- sample_field[[1]]
    Samples <- as.integer(sample_field[2])
    MinEvents <- paste(sample_field[4:length(sample_field)], collapse = " ")

    return(list(
        BatchInfo = batch_info,
        CALInfo = CALInfo,
        Samples = Samples,
        MinEvents = MinEvents
    ))
}

.parse_Xcsv_BatchInfo <- function(Lines) {
    # BatchInfo come in two forms:
    #     field: value
    #     list_field: [field: value]
    # this function creates a list with nested sub-lists where necessary.

    list_fields <- c("ProtocolPlate", "ProtocolMicrosphere")
    parse_list_fields <- function(x) {
        y <- as.list(x[seq(2, length(x), 2)])
        names(y) <- x[seq(1, length(x), 2)]
        y[[3]] <- as.integer(y[[3]])  # Plates/Count
        y
        # lapply(y, type.convert, na.strings = c("", "NA"), as.is = TRUE)
    }

    parse_simple_fields <- function(x) {
        y <- ifelse(any(stringr::str_detect(x, "\\S+")),
                    paste(x, collapse = " "), NA_character_)
        stringr::str_replace(y, "\\s+", " ")
    }

    # Clean up
    Lines <- stringr::str_split(Lines, "\t")  # split into list of vectors
    fields <- vapply(Lines, "[", "", 1L)
    fields <- stringr::str_replace_all(fields, "\\s+", "")
    # remove empty lines
    non_empty_fields <- fields != ""
    Lines <- Lines[non_empty_fields]
    fields <- fields[non_empty_fields]
    simple_fields <- setdiff(fields, list_fields)
    # Final Form!
    BatchInfo <- lapply(Lines, function(x) x[2:max(2, length(x))]) # drop field names
    names(BatchInfo) <- fields
    BatchInfo[simple_fields] <- lapply(BatchInfo[simple_fields], parse_simple_fields)
    BatchInfo[list_fields] <- lapply(BatchInfo[list_fields], parse_list_fields)

    return(BatchInfo)
}

.parse_Xcsv_CALInfo <- function(Lines) {
    # I am not sure about the format. Added a lot of tryCatch to avoid trouble...

    LastCAL <- NA
    if (!is.null(Lines[[1]])) {
        # TODO: list instead of data.frame?
        LastCAL <- .chr_vector_to_df(Lines[[1]][-1], sep = "\t", header = FALSE,
                                     col.names = c("Description", "Results"),
                                     colClasses = "character")
    }

    Classification_Calibrator_Extended <- NA
    Reporter_Calibrator <- NA
    if (!is.null(Lines[[2]])) {
        CALInfo <- Lines[[2]][-1]
        cce <- which(CALInfo == "Classification Calibrator Extended")
        rcal <- which(CALInfo == "Reporter Calibrator")
        Classification_Calibrator_Extended <- tryCatch(
            .chr_vector_to_df(CALInfo[(cce + 1):(rcal - 1)],
                              sep = "\t", header = TRUE),
            error = function(e) NA
        )

        Reporter_Calibrator <- tryCatch({
            rcal <- CALInfo[(rcal + 1):length(CALInfo)]
            lots <- grep("^Lot", rcal)
            rcal <- .splitAt(rcal, lots)
            lapply(rcal, .chr_vector_to_df, sep = "\t", header = TRUE)
        }, error = function(e) NA)

    }

    return(list(
        LastCalibration = LastCAL,
        ClassificationCalibratorExtended = Classification_Calibrator_Extended,
        ReporterCalibrator = Reporter_Calibrator
    ))

}

.parse_Xcsv_results <- function(Lines) {
    # TODO: Break into smaller functions?
    # TODO: Export as tibble?

    # Parameters
    id_vars   <- c("Location", "Sample", "Analyte")
    exprs_var <- c("Median", "Count", "Net MFI", "Result",
                   "Range", "% Recovery", "Comments")
    aver_vars <- c("Avg Net MFI", "Avg Result", "Avg Range",
                   "%CV Replicates", "Standard Expected Concentration")
    well_vars <- c("Dilution Factor")
    bead_vars <- c("Units", "Per Bead Count", "Analysis Types")
    run_vars  <- c("Audit Logs", "Warnings/Errors")

    extract_Location <- function(x) stringr::str_extract(x, "[A-Z]+[0-9]+")

    # Prepare Data
    fields <- vapply(Lines, "[", "", 1)  # first line of the form: "DataType:\t<DataType>"
    fields <- stringr::str_split(fields, "\t")
    fields <- vapply(fields, "[", "", 2)
    names(Lines) <- fields

    # Parse Well Variables
    exprs_var <- intersect(exprs_var, fields)
    Exprs <- lapply(exprs_var,
        function(type) {
            df <- .chr_vector_to_df(Lines[[type]][-1], header = TRUE, sep = "\t")
            df <- df[, setdiff(colnames(df), "Total Events"), drop = FALSE]
            df[["Location"]] <- extract_Location(df[["Location"]])
            analytes <- setdiff(colnames(df), id_vars[1:2])
            # TODO: replace gather_ (drop tidyr dependency?)
            df <- tidyr::gather_(df, id_vars[3], type, analytes)
            colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "")
            df
        }
    )
    Exprs <- Reduce(
        function(...) merge.data.frame(by = id_vars, sort = FALSE, ...),
        Exprs
    )

    # Parse Average Data
    aver_vars <- intersect(aver_vars, fields)
    Average <- NA
    if (length(aver_vars) > 0) {
        Average <- lapply(aver_vars,
            function(type) {
                df <- .chr_vector_to_df(Lines[[type]][-1], header = TRUE, sep = "\t")
                analytes <- setdiff(colnames(df), id_vars[-3])
                # TODO: replace gather_ (drop tidyr dependency?)
                df <- tidyr::gather_(df, id_vars[3], type, analytes)
                colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "")
                df
            }
        )
        Average <- Reduce(
            function(...) merge.data.frame(by = id_vars[-1], sort = FALSE, ...),
            Average
        )
    }

    # Parse Well     Data
    well_vars <- intersect(well_vars, fields)
    Wells <- NA
    if (length(well_vars) > 0) {
        Wells <- lapply(well_vars,
            function(type) {
                df <- .chr_vector_to_df(Lines[[type]][-1], header = TRUE, sep = "\t")
                df[["Location"]] <- extract_Location(df[["Location"]])
                colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "")
                df
            }
        )
        Wells <- Reduce(
            function(...) merge.data.frame(by = id_vars[-3], sort = FALSE, ...),
            Wells
        )
    }

    # Parse Bead Data
    bead_vars <- intersect(bead_vars, fields)
    Beads <- lapply(bead_vars,
        function(type) {
            lines <- Lines[[type]]
            cols <- stringr::str_split(lines[2:(length(lines) - 1)], "\t")
            col_names <- vapply(cols, "[", "", 1)
            col_names <- stringr::str_replace_all(col_names, "[\\s:]+", "")
            df <- lapply(cols, function(x) x[2:max(2, length(x))])
            names(df) <- col_names
            if ("BeadID" %in% col_names)
                df[["BeadID"]] <- as.integer(df[["BeadID"]])
            if (type == "Per Bead Count") {
                df[[3]] <- as.integer(df[[3]])
                names(df)[3] <- "PerBeadCount"
            }

#            df[[type]] <- c(df[[type]], rep(NA, length(df$BeadID) - length(df[[type]])))
            return(data.frame(df, stringsAsFactors = FALSE))
    })
    Beads <- Reduce(
        function(...) merge.data.frame(sort = FALSE, ...),
        Beads
    )

    # RunData
    run_vars <- intersect(run_vars, fields)
    Run <- lapply(run_vars,
        function(type) {
            df <- .chr_vector_to_df(Lines[[type]][-1], header = TRUE, sep = "\t",
                                    colClasses = "character")
            if ("Location" %in% colnames(df))
                df[["Location"]] <- extract_Location(df[["Location"]])
            colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "")
            df
        }
    )
    names(Run) <- run_vars

    return(list(
        Exprs = Exprs,
        Average = Average,
        Wells = Wells,
        Beads = Beads,
        Run = Run
    ))
}

