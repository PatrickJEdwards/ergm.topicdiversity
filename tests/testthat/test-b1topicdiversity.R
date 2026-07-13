test_that("b1topicdiversity matches inverse Simpson calculations", {
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

  inverse_simpson <- ifelse(
    degree > 0,
    1 / rowSums(p^2),
    0
  )

  # MP1 and MP2 both have portfolio (0.75, 0.25), so each contributes
  # 1 / (0.75^2 + 0.25^2) = 1.6. MP3 is an isolate. Total = 3.2.
  expect_equal(unname(inverse_simpson), c(1.6, 1.6, 0), tolerance = 1e-12)

  observed <- summary(
    nw ~ b1topicdiversity(
      c("topic_1", "topic_2"),
      subtract.one = FALSE
    )
  )

  expect_equal(unname(observed), 3.2, tolerance = 1e-10)

  observed_minus_1 <- summary(
    nw ~ b1topicdiversity(
      c("topic_1", "topic_2"),
      subtract.one = TRUE
    )
  )

  expect_equal(unname(observed_minus_1), 1.2, tolerance = 1e-10)
})

test_that("b1topicdiversity has the correct tie-toggle change", {
  y_before <- matrix(
    c(
      1, 1, 0,
      0, 0, 0
    ),
    nrow = 2,
    byrow = TRUE
  )

  y_after <- y_before
  y_after[1, 3] <- 1

  theta <- rbind(
    c(1, 0),
    c(1, 0),
    c(0, 1)
  )

  make_network <- function(y) {
    nw <- network::network(
      y,
      matrix.type = "bipartite",
      directed = FALSE
    )
    nw %v% "topic_1" <- c(rep(NA_real_, nrow(y)), theta[, 1])
    nw %v% "topic_2" <- c(rep(NA_real_, nrow(y)), theta[, 2])
    nw
  }

  nw_before <- make_network(y_before)
  nw_after <- make_network(y_after)

  stat_before <- summary(
    nw_before ~ b1topicdiversity(c("topic_1", "topic_2"))
  )

  stat_after <- summary(
    nw_after ~ b1topicdiversity(c("topic_1", "topic_2"))
  )

  # MP1 changes from (1, 0), effective count 1, to (2/3, 1/3),
  # effective count 1 / ((2/3)^2 + (1/3)^2) = 1.8.
  expect_equal(
    unname(stat_after - stat_before),
    0.8,
    tolerance = 1e-10
  )
})
