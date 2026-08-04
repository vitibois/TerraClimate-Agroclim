# Splitting a raster's pixels into batches ("tiles") for extraction.
#
# The extraction pipelines (scripts/03_... and scripts/04_...) process every
# land pixel of a global TerraClimate raster, but do so in ~100 batches
# rather than all at once, to keep memory usage and per-batch runtime
# manageable. This file factors out that (identical) batching logic.

#' Split all non-NA cells of a raster into `n_tiles` roughly equal batches.
#'
#' @param r a `SpatRaster` (terra), used only for its cell coordinates/ids.
#' @param n_tiles target number of tiles (default 100).
#'
#' @return a data.frame with one row per pixel and columns `x`, `y`,
#'   `pixID` (the raster cell number, i.e. the index used to write values
#'   back with `r[pixID] <- ...`) and `tile` (integer tile id, `1:n_tiles`).
make_pixel_tiles <- function(r, n_tiles = 100) {
  coords <- terra::crds(r, na.rm = TRUE)
  pixID <- terra::cellFromXY(r, coords)
  n <- nrow(coords)
  breaks <- c(seq.int(0, n, by = round(n / (n_tiles - 1))), n)
  tile <- cut(seq_len(n), breaks, labels = FALSE)
  data.frame(x = coords[, "x"], y = coords[, "y"], pixID = pixID, tile = tile)
}

#' Write the list of tile ids to a text file (for scripts that iterate over
#' tiles read back from disk, e.g. when compiling per-tile CSVs into rasters).
#'
#' @param tiles a data.frame as returned by [make_pixel_tiles()].
#' @param path output file path.
write_tile_list <- function(tiles, path) {
  tile_ids <- sort(unique(tiles$tile))
  last_row <- sapply(tile_ids, function(id) max(which(tiles$tile == id)))
  write.table(data.frame(ID = tile_ids, Breaks = last_row), path, row.names = FALSE)
}
