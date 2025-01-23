# trendfilter (development version)

## Breaking changes

- Add an option to standardize `y` (default `TRUE`). This improves computation
  significantly
- Move to x.x.x.9xxx version numbering.
- Add News
- `lambda` must be an argument to `cv_trendfilter()`. Otherwise, passing it 
  optionally fails to forward it correctly through `...` in the CV fits.
