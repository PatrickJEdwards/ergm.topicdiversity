test_that("b1topicdiversity matches a manual calculation", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")

  # Three mode-1 actors and four mode-2 EDMs.
  y <- matrix(
    c(
      1, 1, 0, 0,
      0, 1, 1, 0,
      0, 0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )

  nw <- network::network(
    y,
    matrix.type = "bipartite",
    directed = FALSE
  )

  theta <- rbind(
    c(1.0, 0.0),
    c(0.5, 0.5),
    c(0.0, 1.0),
    c(0.8, 0.2)
  )

  nw %v% "topic_1" <- c(rep(NA_real_, 3), theta[, 1])
  nw %v% "topic_2" <- c(rep(NA_real_, 3), theta[, 2])

  topic_mass <- y %*% theta
  degree <- rowSums(y)
  p <- sweep(topic_mass, 1, pmax(degree, 1), "/")

  plogp <- matrix(0, nrow = nrow(p), ncol = ncol(p))
  positive <- p > 0
  plogp[positive] <- p[positive] * log(p[positive])

  entropy <- -rowSums(plogp)
  effective <- ifelse(degree > 0, exp(entropy), 0)
  expected <- sum(effective)

  observed <- summary(
    nw ~ b1topicdiversity(
      c("topic_1", "topic_2"),
      subtract.one = FALSE
    )
  )

  expect_equal(unname(observed), expected, tolerance = 1e-10)
})
