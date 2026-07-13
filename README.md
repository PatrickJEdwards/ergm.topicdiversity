# ergm.topicdiversity

A minimal standalone extension package defining the binary bipartite ERGM term
`b1topicdiversity()`.

## Statistic

For every mode-1 actor, the term sums the continuous topic compositions of the
mode-2 neighbors, normalizes the resulting portfolio, calculates its inverse
Simpson effective number of topics, and sums those actor-level counts over all
mode-1 actors.

For actor `i`, the measure is:

```text
D_i = 1 / sum_k p_ik^2
```

where `p_ik` is the proportion of the actor's aggregate topic mass assigned to
topic `k`. This is mathematically the same formula as the Laakso-Taagepera
effective number of parties, with topics replacing parties.

The default, `subtract.one = FALSE`, returns the ordinary inverse Simpson
effective topic count for every nonisolate and zero for isolates. With
`subtract.one = TRUE`, the contribution is `effective topics - 1` for
nonisolates and zero for isolates. The shifted version centers a perfectly
one-topic portfolio at zero. Because an individual EDM can itself have a mixed
topic composition, its first tie can still add positive shifted diversity.

## Install

The package contains compiled C code, so Windows requires a working Rtools
installation compatible with the installed R version. From the parent
directory:

```r
install.packages(c("devtools", "roxygen2", "testthat"))
devtools::document("ergm.topicdiversity")
devtools::install(
  "ergm.topicdiversity",
  dependencies = FALSE,
  upgrade = FALSE,
  reload = FALSE
)
```

After changing C code on Windows, restart R before reinstalling so the old DLL
is not locked.

## Use

Each network must store the same named topic attributes on its mode-2 vertices.
For example:

```r
topic_attributes <- paste0("umbrella_topic_", 1:17)

summary(
  network_list[[1]] ~
    b1topicdiversity(topic_attributes, subtract.one = FALSE)
)
```

The coefficient name will be:

```text
b1topicdiversity.inverse_simpson
```

or, with `subtract.one = TRUE`:

```text
b1topicdiversity.inverse_simpson_minus_1
```

The intended dependent-term interaction is:

```r
b1topicdiversity(topic_attributes, subtract.one = FALSE) :
  b1factor("openlyLGBT_d", levels = I(1))
```

Current `ergm` versions reject interactions involving dyad-dependent terms by
default. Pass `term.options = list(interact.dependent = "silent")` through the
`control.ergm()` or `ergm_model()`/`ergmMPLE()` call path used by the modified
BTERGM package after validating that the interaction statistic is the intended
actor-specific product.

## Test

```r
devtools::test("ergm.topicdiversity")
```

The unit tests include an uneven topic portfolio specifically chosen to
distinguish inverse Simpson diversity from Shannon effective diversity.
