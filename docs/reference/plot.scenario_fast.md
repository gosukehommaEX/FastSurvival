# Plot Method for Scenario Sets

Draws, for each scenario in a `scenario_fast` object, the analytic
control and experimental survival curves together with the
piecewise-constant hazard ratio on a secondary axis, in the style of the
stylized nonproportional-hazards figures used in the design literature.
The curves are computed analytically from the piecewise-exponential
parameters, so no simulation is involved.

## Usage

``` r
# S3 method for class 'scenario_fast'
plot(
  x,
  which = NULL,
  tmax = NULL,
  hr_max = NULL,
  n_time = 481L,
  mfrow = NULL,
  col_control = "#1b9e77",
  col_treat = "#d95f02",
  col_hr = "#3060c0",
  lty_control = 1,
  lty_treat = 1,
  lty_hr = 2,
  lwd_control = 2,
  lwd_treat = 2,
  lwd_hr = 1.5,
  base_size = 1,
  xlab = "Time",
  ylab = "Survival probability",
  legend_panel = 1L,
  legend_pos = "topright",
  ...
)
```

## Arguments

- x:

  A `scenario_fast` object from
  [`gen_scenario_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/gen_scenario_fast.md).

- which:

  Optional integer or character vector selecting a subset of scenarios
  to draw.

- tmax:

  Upper limit of the time axis. Defaults to a value derived from the
  scenarios.

- hr_max:

  Upper limit of the hazard-ratio axis. Defaults to a value derived from
  the scenarios.

- n_time:

  Number of points on the time grid.

- mfrow:

  Optional length-two panel layout `c(rows, cols)`. Defaults to a
  near-square grid.

- col_control, col_treat, col_hr:

  Colors for the control survival curve, the experimental survival
  curve, and the hazard ratio.

- lty_control, lty_treat, lty_hr:

  Line types for the control survival curve, the experimental survival
  curve, and the hazard ratio.

- lwd_control, lwd_treat, lwd_hr:

  Line widths for the control survival curve, the experimental survival
  curve, and the hazard ratio.

- base_size:

  Relative size of all text (titles, axis labels, axis annotation, and
  the legend). The default `1` reproduces the standard sizes; values
  above one enlarge the text and values below one shrink it.

- xlab, ylab:

  Axis labels for the survival panels.

- legend_panel:

  Index (within the drawn subset) of the panel that carries the legend.

- legend_pos:

  Keyword giving the legend position, passed to
  [`legend`](https://rdrr.io/r/graphics/legend.html) (for example
  `"topright"` or `"bottomleft"`).

- ...:

  Further arguments, currently ignored.

## Value

The input `x`, invisibly.

## Details

Each scenario is shown in its own panel arranged on a grid. The hazard
ratio is rescaled into the survival plotting region and read on the
right-hand axis; a reference line marks a hazard ratio of one. A
single-arm scenario (no treatment specification) is drawn as a single
survival curve with no hazard ratio.

## See also

[`gen_scenario_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/gen_scenario_fast.md),
[`print.scenario_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/print.scenario_fast.md)
