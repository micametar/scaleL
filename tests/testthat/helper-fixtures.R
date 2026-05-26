make_synthetic_data <- function() {
  data.frame(
    study_id   = c("A", "A", "B", "C"),
    measure_id = c("m1", "m2", "m1", "m1"),
    mean       = c(3.4, 3.8, 5.1, 3.7),
    sd         = c(0.7, 0.6, 1.1, 0.8),
    n          = c(120, 120, 90, 200),
    L          = c(5, 5, 7, NA),
    scale_min  = c(1, 1, 1, 1),
    stringsAsFactors = FALSE
  )
}
