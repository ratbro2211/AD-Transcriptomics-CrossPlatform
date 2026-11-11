# Download public GEO files whose MD5 matches the supplied Drive inputs.
# This does NOT convert normalized RNA-seq values into raw counts.
if (!file.exists("supplied_project/README.md"))
  stop("Run from the repository root.")
options(timeout = max(600, getOption("timeout", 60)))
files <- data.frame(
  path = c(
    "supplied_project/micro_array/GSE132903_matrix.txt.gz",
    "supplied_project/GSE104704/GSE104704_RNA-Seq_Table.txt.gz"),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132903/matrix/GSE132903_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104704/suppl/GSE104704_RNA-Seq_Table.txt.gz"),
  md5 = c("692ce03f44a52a25611357a2aa3e37d8",
          "e08e5d04bff8189081c767ae68529512"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(files))) {
  path <- files$path[i]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    if (unname(tools::md5sum(path)) == files$md5[i]) {
      message("Already verified: ", path)
      next
    }
    stop("An existing file has a different checksum; it was not overwritten: ", path)
  }
  tmp <- tempfile(fileext = ".gz")
  status <- tryCatch(download.file(files$url[i], tmp, mode = "wb", method = "libcurl"),
                     error = function(e) stop("Download failed: ", conditionMessage(e)))
  if (status != 0L || unname(tools::md5sum(tmp)) != files$md5[i]) {
    unlink(tmp)
    stop("Download/checksum verification failed: ", path)
  }
  if (!file.copy(tmp, path, overwrite = FALSE)) stop("Cannot copy verified download to ", path)
  unlink(tmp)
  message("Downloaded and verified: ", path)
}
source <- "supplied_project/GSE104704/GSE104704_RNA-Seq_Table.txt.gz"
duplicate <- "supplied_project/RNA_seq/GSE104704/GSE104704_RNA-Seq_Table.txt.gz"
dir.create(dirname(duplicate), recursive = TRUE, showWarnings = FALSE)
if (file.exists(duplicate)) {
  if (unname(tools::md5sum(duplicate)) != files$md5[2])
    stop("Existing duplicate has a different checksum; not overwritten.")
} else if (!file.copy(source, duplicate, overwrite = FALSE)) {
  stop("Could not restore the duplicate RNA-seq input path.")
}
message("Public inputs restored. RNA-seq supplementary values are normalized, NOT raw counts.")
