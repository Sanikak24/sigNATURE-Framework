#!/usr/bin/env Rscript

# Shared setup for the CD4 and CD8 reference atlases used by the BCC and NSCLC
# workflows.

options(timeout = max(3600, getOption("timeout")))

reference_dir <- file.path("data", "reference")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

download_reference <- function(name, url, destination, expected_md5) {
  if (file.exists(destination) && file.info(destination)$size > 0) {
    message("Skipping download; nonempty file already exists: ", destination)
  } else {
    partial_file <- paste0(destination, ".part")

    if (file.exists(partial_file)) {
      message("Partial download found; resuming from: ", partial_file)
    }

    message("Downloading ", name, " from ", url)

    curl_path <- Sys.which("curl")
    if (nzchar(curl_path)) {
      download_status <- system2(
        command = curl_path,
        args = c(
          "--location",
          "--fail",
          "--retry", "5",
          "--retry-delay", "10",
          "--retry-all-errors",
          "--continue-at", "-",
          "--output", shQuote(partial_file),
          shQuote(url)
        )
      )

      if (!identical(download_status, 0L)) {
        stop(
          "Failed to download ", name,
          "; curl returned status ", download_status,
          ". Partial download retained at: ", partial_file
        )
      }
    } else {
      message("System curl is unavailable; using download.file().")
      download_status <- tryCatch(
        download.file(
          url = url,
          destfile = partial_file,
          mode = "wb",
          method = "libcurl",
          quiet = FALSE
        ),
        error = function(error) {
          stop(
            "Failed to download ", name, ": ", conditionMessage(error),
            ". Partial download retained at: ", partial_file
          )
        }
      )

      if (!identical(download_status, 0L)) {
        stop(
          "Failed to download ", name,
          "; download.file returned status ", download_status,
          ". Partial download retained at: ", partial_file
        )
      }
    }

    if (!file.exists(partial_file) || file.info(partial_file)$size == 0) {
      stop(
        "Downloaded ", name,
        " is missing or empty. Expected partial file: ", partial_file
      )
    }

    downloaded_md5 <- unname(tools::md5sum(partial_file))
    if (!identical(tolower(downloaded_md5), tolower(expected_md5))) {
      stop(
        "MD5 mismatch for downloaded ", name, ". Expected ",
        expected_md5, " but received ", downloaded_md5,
        ". Partial download retained at: ", partial_file
      )
    }

    if (file.exists(destination) && unlink(destination) != 0) {
      stop("Could not remove empty destination file: ", destination)
    }

    if (!file.rename(partial_file, destination)) {
      stop("Validated ", name, " could not be moved to: ", destination)
    }
    message("Saved validated ", name, " to: ", destination)
  }

  observed_md5 <- unname(tools::md5sum(destination))
  if (!identical(tolower(observed_md5), tolower(expected_md5))) {
    stop(
      "MD5 mismatch for ", destination, ". Expected ",
      expected_md5, " but received ", observed_md5, "."
    )
  }

  message("MD5 validation passed for: ", destination)
  invisible(destination)
}

# MD Anderson T Cell Map reference objects:
# https://singlecell.mdanderson.org/TCM/
reference_inputs <- list(
  list(
    name = "CD8 reference",
    url = "https://singlecell.mdanderson.org/TCM/download/CD8",
    destination = file.path(reference_dir, "CD8_Obj_for_mapping.rds"),
    md5 = "52c6daf010f16020a8b36fb40bb95618"
  ),
  list(
    name = "CD4 reference",
    url = "https://singlecell.mdanderson.org/TCM/download/CD4",
    destination = file.path(reference_dir, "CD4_Obj_for_mapping.rds"),
    md5 = "69b4164d2672b39d70c61e574bcdcaa4"
  )
)

for (input in reference_inputs) {
  download_reference(
    name = input$name,
    url = input$url,
    destination = input$destination,
    expected_md5 = input$md5
  )
}

message("Shared reference-atlas setup complete.")
