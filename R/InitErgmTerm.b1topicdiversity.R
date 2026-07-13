# File R/InitErgmTerm.b1topicdiversity.R
#
# This package extends the Statnet ergm term API. It is distributed under
# GPL-3 and retains the Statnet attribution described at
# https://statnet.org/attribution .

#' @templateVar name b1topicdiversity
#'
#' @title Shannon Effective Topic Diversity of Mode-1 Neighborhoods
#'
#' @description
#' This binary bipartite ERGM term adds one statistic equal to the sum,
#' across first-mode actors, of the Shannon effective number of topics
#' represented among the second-mode vertices to which each actor is tied.
#'
#' For actor i, let theta[j,k] denote mode-2 vertex j's proportion on topic k.
#' The actor's accumulated topic mass is s[i,k] = sum_j y[i,j] theta[j,k].
#' After normalizing by actor degree d[i], its portfolio proportions are
#' p[i,k] = s[i,k] / d[i]. The effective number of topics is
#' exp(-sum_k p[i,k] log(p[i,k])).
#'
#' @usage
#' # binary: b1topicdiversity(topics, subtract.one = FALSE)
#'
#' @param topics A mode-2 vertex-attribute specification resolving to a
#'   numeric matrix. Rows must correspond to mode-2 vertices and columns to
#'   topics. In practice, this can be a character vector containing the names
#'   of topic-proportion vertex attributes or a formula such as
#'   `~ cbind(topic_1, topic_2, topic_3)`.
#' @param subtract.one Logical. If `FALSE` (the default), a nonisolated actor
#'   contributes the ordinary Shannon effective topic count D, while an
#'   isolate contributes zero. If `TRUE`, a nonisolated actor contributes
#'   D - 1. This centers a perfectly one-topic portfolio at zero and removes
#'   the unavoidable baseline value of one from every nonisolate. A single
#'   mixed-topic EDM can still contribute more than zero under this shifted
#'   version.
#'
#' @details
#' Topic rows are normalized internally to sum to one. Therefore, all supplied
#' topic attributes must be finite and nonnegative, and every mode-2 vertex
#' must have positive total topic mass.
#'
#' The term is dyad-dependent because the change associated with an MP-EDM tie
#' depends on the MP's existing topic portfolio.
#'
#' @name b1topicdiversity-ergmTerm
#' @aliases b1topicdiversity
#' @concept bipartite
#' @concept undirected
#' @concept binary
NULL

InitErgmTerm.b1topicdiversity <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(
    nw,
    arglist,
    directed = FALSE,
    bipartite = TRUE,
    varnames = c("topics", "subtract.one"),
    vartypes = c(ERGM_VATTR_SPEC, "logical"),
    required = c(TRUE, FALSE),
    defaultvalues = list(NULL, FALSE)
  )

  if (length(a$subtract.one) != 1L || is.na(a$subtract.one)) {
    ergm_Init_stop("`subtract.one` must be one non-missing logical value.")
  }

  theta <- ergm_get_vattr(
    a$topics,
    nw,
    accept = "numeric",
    bip = "b2",
    multiple = "matrix"
  )

  theta <- as.matrix(theta)
  storage.mode(theta) <- "double"

  n_mode1 <- nw %n% "bipartite"
  n_mode2 <- network.size(nw) - n_mode1

  if (nrow(theta) != n_mode2) {
    ergm_Init_stop(
      "`topics` must provide one row for every mode-2 vertex."
    )
  }

  if (ncol(theta) < 2L) {
    ergm_Init_stop("`topics` must contain at least two topic columns.")
  }

  if (any(!is.finite(theta))) {
    ergm_Init_stop("All topic proportions must be finite and non-missing.")
  }

  if (any(theta < 0)) {
    ergm_Init_stop("Topic proportions cannot be negative.")
  }

  row_mass <- rowSums(theta)

  if (any(row_mass <= 0)) {
    ergm_Init_stop(
      "Every mode-2 vertex must have positive total topic mass."
    )
  }

  # Normalize each EDM's supplied topic weights to a proper composition.
  theta <- theta / row_mass

  statistic_name <- if (isTRUE(a$subtract.one)) {
    "b1topicdiversity.effective_minus_1"
  } else {
    "b1topicdiversity.effective"
  }

  list(
    name = "b1topicdiversity",
    coef.names = statistic_name,
    pkgname = "ergm.topicdiversity",

    # C expects actor-major rows: all K topics for EDM 1, then all K topics
    # for EDM 2, and so forth. R is column-major, so transpose first.
    inputs = as.double(t(theta)),

    # iinputs[1] = K; iinputs[2] = subtract-one indicator.
    iinputs = c(
      as.integer(ncol(theta)),
      as.integer(isTRUE(a$subtract.one))
    ),

    dependence = TRUE,
    emptynwstats = 0
  )
}
