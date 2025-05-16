library(tidyverse)
library(kableExtra)
library(survminer)
library(survival)
library(patchwork)
library(tictoc)

### Examples using the toy data
# Parameter set
mst.C = 12
HR = 0.5
n.j = 5
tau = 5
E = 8
seed = 123

# Set seed
set.seed(seed)

# Total sample size
N = 2 * n.j

# Assumed hazard rates
hazard.C = log(2) / mst.C
hazard.T = HR * hazard.C

# toy.data
toy.data = tibble(
  j = factor(rep(c('T', 'C'), each = n.j), levels = c('T', 'C')),
  u.ij = runif(N, min = 0, max = tau),
  t.ij = '+'(
    (j == 'T') * rexp(n.j, rate = hazard.T),
    (j == 'C') * rexp(n.j, rate = hazard.C)
  ),
  Total = u.ij + t.ij
) %>%
  arrange(Total) %>%
  mutate(
    analysis.time = Total[E],
    event = ifelse(Total > analysis.time, 0, 1),
    t.ij = ifelse(Total > analysis.time, pmax(0, analysis.time - u.ij), t.ij),
    j = relevel(j, ref = 'C')
  ) %>%
  select(
    -c(Total, analysis.time)
  ) %>%
  arrange(t.ij) %>%
  mutate(
    j = ifelse(j == 'T', 1, 0)
  ) %>%
  mutate(
    O.T = event * j,
    O.C = event * (1 - j),
    E.T = event * rev(cumsum(rev(j))) / (n():1),
    E.C = event * ((n():1) - rev(cumsum(rev(j)))) / (n():1),
    V.T = '/'(
      ((n():1) - rev(cumsum(rev(j)))) * rev(cumsum(rev(j))) * event * ((n():1) - event),
      (n():1) ^ 2 * ((n():1) - 1)
    )
  ) %>%
  mutate(
    j = relevel(factor(ifelse(j == 1, 'T', 'C'), levels = c('T', 'C')), ref = 'C'),
    across(where(is.numeric), coalesce, 0)
  )

# Create table1
table1 = toy.data %>%
  mutate_if(is.numeric, round, 2) %>%
  arrange(t.ij) %>%
  select(
    j, u.ij, t.ij, event
  ) %>%
  kbl(align = 'c', format = 'latex', booktabs = T)

# Kaplan-Meier plot
kaplan.Meier.data = survfit(Surv(t.ij, event) ~ j, data = toy.data)
kaplan.Meier.plot = ggsurvplot(
  kaplan.Meier.data,
  dat = toy.data,
  legend.title = '',
  legend.labs = c('C', 'T'),
  size = 2.5,
  xlab = 'Survival time',
  ylab = 'Survival probability',
  censor = T,
  censor.shape = '|',
  censor.size = 15,
  ggtheme = theme_light()
)$plot +
  scale_color_manual(
    values = c(
      'C' = '#939597',
      'T' = '#D91E49'
    ),
    labels =  c('C', 'T')
  ) +
  theme(
    strip.text.x = element_text(size = 50),
    strip.text.y = element_text(size = 50),
    text = element_text(size = 50),
    legend.key.width = unit(4, 'cm'),
    legend.text = element_text(size = 50),
    legend.position = 'bottom'
  )
ggsave(
  file = paste('figure1', 'eps', sep = '.'),
  plot = kaplan.Meier.plot,
  dpi = 800,
  width = 16,
  height = 16
)

## Estimate hazard ratio
# Cox model
HR.cox = coxph(Surv(t.ij, event) ~ j, data = toy.data) %>%
  coefficients() %>%
  exp() %>%
  as.numeric()
# Pearson-year method
HR.PY = toy.data %>%
  summarise(
    '/'(
      sum(event[j == 'T']) / sum(t.ij[j == 'T']),
      sum(event[j == 'C']) / sum(t.ij[j == 'C'])
    )
  ) %>%
  pull()
# Median survival time method
HR.MST = toy.data %>%
  summarise(
    median(t.ij[j == 'C']) / median(t.ij[j == 'T'])
  ) %>%
  pull()

# Create table2
table2 = toy.data %>%
  mutate_if(is.numeric, round, 2) %>%
  arrange(t.ij) %>%
  kbl(align = 'c', format = 'latex', booktabs = T)

# Pike method
HR.Pike = toy.data %>%
  summarise(
    (sum(O.T) * sum(E.C)) / (sum(O.C) * sum(E.T))
  ) %>%
  pull()
# Peto method
HR.Peto = toy.data %>%
  summarise(
    exp((sum(O.T) - sum(E.T)) / sum(V.T))
  ) %>%
  pull()
# Log-rank test method
HR.LR = toy.data %>%
  summarise(
    Z.LR = (sum(O.T) - sum(E.T)) / sqrt(sum(V.T)),
    HR = exp(2 * Z.LR / sqrt(E))
  ) %>%
  pull(HR)

### Operating characteristics comparison
nsim = 10000
## Scenario
# Parameter set
Scenario = LETTERS[1:6]
mst.C = 6
HR = c(0.3, 1, 0.5, 1, 0.7, 1)
n.j = c(22, 22, 59, 59, 204, 204)
tau = 12
E = c(29, 29, 88, 88, 331, 331)
seed = 123

results = lapply(seq(Scenario), function(i) {
  ## Create dataset
  # Set seed
  set.seed(seed)

  # Total sample size
  N = 2 * n.j[i]

  # Assumed hazard rates
  hazard.C = log(2) / mst.C
  hazard.T = HR[i] * hazard.C

  # Dataset
  dataset = tibble(
    sim = rep(seq(nsim), each = N),
    j = rep(factor(rep(c('T', 'C'), each = n.j[i]), levels = c('T', 'C')), nsim),
    u.ij = runif(N * nsim, min = 0, max = tau),
    t.ij = '+'(
      (j == 'T') * rexp(n.j[i] * nsim, rate = hazard.T),
      (j == 'C') * rexp(n.j[i] * nsim, rate = hazard.C)
    ),
    Total = u.ij + t.ij
  ) %>%
    group_by(sim) %>%
    arrange(Total) %>%
    mutate(
      analysis.time = Total[E[i]],
      event = ifelse(Total > analysis.time, 0, 1),
      t.ij = ifelse(Total > analysis.time, pmax(0, analysis.time - u.ij), t.ij),
      j = relevel(j, ref = 'C')
    ) %>%
    select(
      -c(Total, analysis.time)
    ) %>%
    arrange(t.ij) %>%
    mutate(
      j = ifelse(j == 'T', 1, 0)
    ) %>%
    mutate(
      O.T = event * j,
      O.C = event * (1 - j),
      E.T = event * rev(cumsum(rev(j))) / (n():1),
      E.C = event * ((n():1) - rev(cumsum(rev(j)))) / (n():1),
      V.T = '/'(
        ((n():1) - rev(cumsum(rev(j)))) * rev(cumsum(rev(j))) * event * ((n():1) - event),
        (n():1) ^ 2 * ((n():1) - 1)
      )
    ) %>%
    mutate(
      j = relevel(factor(ifelse(j == 1, 'T', 'C'), levels = c('T', 'C')), ref = 'C'),
      across(where(is.numeric), coalesce, 0)
    )

  # Computing time comparison
  tic.clearlog()
  tic('total time')
  # Cox model
  tic('HR.cox')
  HR.cox = dataset %>%
    reframe(
      Cox.result = list(coxph(Surv(t.ij, event) ~ j)),
      HR = as.numeric(exp(Cox.result[[1]][['coefficients']]))
    ) %>%
    pull(HR)
  toc(log = TRUE)
  # Median survival time method
  tic('HR.MST')
  HR.MST = dataset %>%
    reframe(
      median(t.ij[j == 'C']) / median(t.ij[j == 'T'])
    ) %>%
    pull()
  toc(log = TRUE)
  # Pearson-year method
  tic('HR.PY')
  HR.PY = dataset %>%
    reframe(
      '/'(
        sum(event[j == 'T']) / sum(t.ij[j == 'T']),
        sum(event[j == 'C']) / sum(t.ij[j == 'C'])
      )
    ) %>%
    pull()
  toc(log = TRUE)
  # Pike method
  tic('HR.Pike')
  HR.Pike = dataset %>%
    reframe(
      (sum(O.T) * sum(E.C)) / (sum(O.C) * sum(E.T))
    ) %>%
    pull()
  toc(log = TRUE)
  # Peto method
  tic('HR.Peto')
  HR.Peto = dataset %>%
    reframe(
      exp((sum(O.T) - sum(E.T)) / sum(V.T))
    ) %>%
    pull()
  toc(log = TRUE)
  # Log-rank test method
  tic('HR.LR')
  HR.LR = dataset %>%
    reframe(
      Z.LR = (sum(O.T) - sum(E.T)) / sqrt(sum(V.T)),
      HR = exp(2 * Z.LR / sqrt(E[i]))
    ) %>%
    pull(HR)
  toc(log = TRUE)
  toc(log = TRUE)
  time.log = tic.log(format = TRUE)
  computing.time = matrix(unlist(time.log), ncol = 1)

  # Evaluate bias
  bias = tibble(
    Scenario = paste('Scenario', Scenario[i], sep = ' '),
    #Scenario = Scenario[i],
    Method = factor(
      rep(c('MST', 'PY', 'Pike', 'Peto', 'LR'), each = nsim),
      levels = c('MST', 'PY', 'Pike', 'Peto', 'LR')
    ),
    HR.Cox = rep(HR.cox, 5),
    HR = c(HR.MST, HR.PY, HR.Pike, HR.Peto, HR.LR),
    Bias = HR - HR.Cox
  )

  # Return result for a scenario
  return(
    list(computing.time = computing.time, bias = bias)
  )
})

# Summary of computing time
summary.computing.time = Reduce(
  'cbind',
  lapply(seq(Scenario), function(i) {
    results[[i]][['computing.time']]
  })
)

# Summary of bias
summary.bias = Reduce(
  'bind_rows',
  lapply(seq(Scenario), function(i) {
    results[[i]][['bias']]
  })
)

# Boxplot of the bias
Box.plot.bias = ggplot(
  summary.bias,
  aes(x = Method, y = Bias, color = Method)
) +
  geom_boxplot(size = 1.5) +
  theme_bw() +
  facet_grid(
    . ~ Scenario
  ) +
  scale_color_manual(
    values = c(
      'MST' = '#A62B4E',
      'PY' = '#F0B323',
      'Pike' = '#004C97',
      'Peto' = '#658D1B',
      'LR' = '#A20067'
    ),
    labels =  c('MST', 'PY', 'Pike', 'Peto', 'LR')
  ) +
  theme(
    strip.text.x = element_text(size = 50),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_text(size = 50),
    text = element_text(size = 50),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 50),
    legend.title = element_blank(),
    legend.position = 'bottom',
    panel.spacing = unit(0.5, 'lines')
  )
ggsave(
  file = paste('figure2', 'eps', sep = '.'),
  plot = Box.plot.bias,
  dpi = 800,
  width = 32,
  height = 16
)
