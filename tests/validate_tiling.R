# Regression check for R/tiling.R: compares make_pixel_tiles() pixel ids and
# tile assignment against the original inline logic from
# scripts/03_.../04_... (both variants used across the legacy scripts).
#
# Run with:
#   Rscript tests/validate_tiling.R

source("R/tiling.R")
library(terra)

set.seed(1)
r <- rast(nrows = 37, ncols = 53, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
values(r) <- rnorm(ncell(r))
r[sample(ncell(r), round(ncell(r) * 0.2))] <- NA # simulate ocean/NA pixels

## ---- legacy variant A (scripts/01-...v02.R style) --------------------------
r.coords_a <- crds(r)
pixID_a <- terra::extract(r, r.coords_a, cells = TRUE)[, 1]

## ---- legacy variant B (scripts/01_bis-..._v02.R style) ---------------------
r.coords_b <- crds(r, na.rm = TRUE)
pixID_b <- terra::cellFromXY(r, r.coords_b)

## ---- new shared helper ------------------------------------------------------
tiles <- make_pixel_tiles(r, n_tiles = 100)

stopifnot(
  "pixID from the two legacy variants must be identical" =
    identical(pixID_a, pixID_b),
  "new helper pixID must match legacy pixID" =
    identical(tiles$pixID, pixID_b),
  "new helper coordinates must match legacy coordinates" =
    isTRUE(all.equal(tiles$x, r.coords_b[, "x"])) &&
    isTRUE(all.equal(tiles$y, r.coords_b[, "y"])),
  "every pixel must be assigned to exactly one tile" =
    !anyNA(tiles$tile),
  "tile ids must be contiguous 1:n_tiles (or fewer for small rasters)" =
    all(sort(unique(tiles$tile)) == seq_len(max(tiles$tile)))
)

cat("OK: make_pixel_tiles() matches both legacy pixID computations, tiles ",
    max(tiles$tile), " groups for ", nrow(tiles), " pixels.\n", sep = "")
