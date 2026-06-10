# Plot Method for Kaplan-Meier Curves

Draws the Kaplan-Meier curves of a `kmcurve_fast` object. With
`rmst = TRUE` the area under each curve up to `tau` is shaded to
visualize the restricted mean survival time. With `hr = TRUE` a second
panel shows the smoothed time-varying hazard ratio together with the
constant Cox proportional-hazards estimate and, optionally, the average
hazard ratio, in the style of stylized nonproportional-hazards figures.
For a single-group object only one survival curve is drawn; `hr = TRUE`
is ignored because a hazard ratio requires two groups.

## Usage

``` r
# S3 method for class 'kmcurve_fast'
plot(
  x,
  hr = FALSE,
  rmst = FALSE,
  tau = NULL,
  bw = NULL,
  n_grid = 200L,
  conf.int = TRUE,
  hr_max = NULL,
  cox_line = TRUE,
  ahr_line = FALSE,
  mark_time = TRUE,
  conf.level = 0.95,
  conf.type = "log",
  col_control = "#1b9e77",
  col_treat = "#d95f02",
  col_hr = "#3060c0",
  col_cox = "black",
  col_ahr = "grey40",
  lty_control = 1,
  lty_treat = 1,
  lty_hr = 1,
  lwd_control = 2,
  lwd_treat = 2,
  lwd_hr = 2,
  base_size = 1,
  legend_pos = "topright",
  xlab = "Time",
  ylab = "Survival probability",
  main_km = "Kaplan-Meier",
  main_hr = "Hazard ratio",
  group_labels = NULL,
  ...
)
```

## Arguments

- x:

  A `kmcurve_fast` object from
  [`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md).

- hr:

  Logical. If `TRUE`, add a hazard-ratio panel to the right.

- rmst:

  Logical. If `TRUE`, shade the area under each curve up to `tau` and
  annotate the restricted mean survival time difference.

- tau:

  Upper limit used for the restricted mean survival time and the average
  hazard ratio. Defaults to the smaller of the two groups' largest
  observed times.

- bw:

  Kernel bandwidth for the smoothed hazard ratio. Defaults to a value
  derived from the spread of the event times.

- n_grid:

  Number of grid points for the smoothed hazard ratio.

- conf.int:

  Logical. If `TRUE`, add pointwise confidence limits to the survival
  curves.

- hr_max:

  Upper limit of the hazard-ratio axis. Defaults to a value derived from
  the estimated hazard ratio.

- cox_line:

  Logical. If `TRUE`, add the constant Cox hazard ratio as a horizontal
  line in the hazard-ratio panel.

- ahr_line:

  Logical. If `TRUE`, add the average hazard ratio as a horizontal line
  in the hazard-ratio panel.

- mark_time:

  Logical. If `TRUE`, mark censoring times on the curves.

- conf.level:

  Confidence level for the survival confidence limits.

- conf.type:

  Type of survival confidence limits: one of `"plain"`, `"log"`
  (default), or `"log-log"`, matching
  [`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md).

- col_control, col_treat, col_hr, col_cox, col_ahr:

  Colors for the control curve, the treatment curve, the smoothed hazard
  ratio, the Cox line, and the average-hazard-ratio line.

- lty_control, lty_treat:

  Line types for the control and treatment curves.

- lty_hr:

  Line type for the smoothed hazard ratio.

- lwd_control, lwd_treat:

  Line widths for the control and treatment curves.

- lwd_hr:

  Line width for the lines in the hazard-ratio panel.

- base_size:

  Relative size of all text.

- legend_pos:

  Keyword giving the legend position in both panels, passed to
  [`legend`](https://rdrr.io/r/graphics/legend.html) (for example
  `"topright"` or `"bottomleft"`).

- xlab, ylab:

  Axis labels for the survival panel.

- main_km, main_hr:

  Panel titles.

- group_labels:

  Optional length-two character vector of legend labels for the control
  and treatment groups.

- ...:

  Further arguments, currently ignored.

## Value

The input `x`, invisibly.

## Details

The smoothed hazard ratio is obtained by kernel smoothing the
Nelson-Aalen increments of each group. For group `j` the kernel hazard
estimate is the sum over event times of
`(1 / bw) K((t - t_i) / bw) (d_i / Y_i)`, with an Epanechnikov kernel
`K`, and the hazard ratio is the treatment estimate divided by the
control estimate. The bandwidth `bw` controls the width of the averaging
window and materially affects the result, so it is exposed as an
argument. To limit boundary bias the hazard ratio is drawn over an
interior time window, and points where the control hazard estimate is
effectively zero are dropped.

The constant Cox hazard ratio is taken from
[`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md),
the average hazard ratio from
[`ahr_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md),
and the restricted mean survival time difference from
[`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md).

## See also

[`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md),
[`print.kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/print.kmcurve_fast.md)
