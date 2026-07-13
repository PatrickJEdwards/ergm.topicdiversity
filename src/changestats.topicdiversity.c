/* File src/changestats.topicdiversity.c
 *
 * This package extends the Statnet ergm term API. It is distributed under
 * GPL-3 and retains the Statnet attribution described at
 * https://statnet.org/attribution .
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"

#include <stddef.h>

/*
 * Private storage layout:
 *
 *   state[0, ..., n1 - 1]
 *       Current degree of each mode-1 actor.
 *
 *   state[n1, ..., n1 + n1*K - 1]
 *       Current accumulated topic mass for each mode-1 actor, stored
 *       actor-major: actor 1 topics 1..K, actor 2 topics 1..K, etc.
 */

/*
 * Calculate the inverse Simpson effective number of topics:
 *
 *   D = 1 / sum_k p_k^2,
 *
 * where p_k is the actor's normalized attention to topic k.
 */
static double inverse_simpson_topic_count(
  const double *topic_mass,
  double degree,
  int K,
  int subtract_one
) {
  if (degree <= 0.0) {
    return 0.0;
  }

  double sum_squared_proportions = 0.0;

  for (int k = 0; k < K; ++k) {
    const double mass = topic_mass[k];

    /* Ignore exact zeros and tiny roundoff artifacts. */
    if (mass > 1e-12) {
      const double p = mass / degree;
      sum_squared_proportions += p * p;
    }
  }

  /* A nonisolated actor must have positive total topic mass. */
  if (sum_squared_proportions <= 0.0) {
    return 0.0;
  }

  const double effective = 1.0 / sum_squared_proportions;
  return subtract_one ? effective - 1.0 : effective;
}

static double inverse_simpson_topic_count_after_toggle(
  const double *topic_mass,
  double degree,
  const double *edm_topics,
  int K,
  int sign,
  int subtract_one
) {
  const double new_degree = degree + (double) sign;

  if (new_degree <= 0.0) {
    return 0.0;
  }

  double sum_squared_proportions = 0.0;

  for (int k = 0; k < K; ++k) {
    double new_mass = topic_mass[k] + (double) sign * edm_topics[k];

    /* Clamp only very small negative roundoff after removing a tie. */
    if (new_mass < 0.0 && new_mass > -1e-10) {
      new_mass = 0.0;
    }

    if (new_mass > 1e-12) {
      const double p = new_mass / new_degree;
      sum_squared_proportions += p * p;
    }
  }

  if (sum_squared_proportions <= 0.0) {
    return 0.0;
  }

  const double effective = 1.0 / sum_squared_proportions;
  return subtract_one ? effective - 1.0 : effective;
}

/* Initialize degree and topic-mass storage from the current network. */
I_CHANGESTAT_FN(i_b1topicdiversity) {
  const int K = IINPUT_PARAM[0];
  const Vertex n1 = BIPARTITE;
  const size_t storage_length = (size_t) n1 * (size_t) (K + 1);

  ALLOC_STORAGE(storage_length, double, state);

  /* Do not depend on allocator-specific zero initialization. */
  for (size_t z = 0; z < storage_length; ++z) {
    state[z] = 0.0;
  }

  double *degrees = state;
  double *topic_mass = state + n1;

  /* For bipartite ergm networks, tail is mode 1 and head is mode 2. */
  EXEC_THROUGH_NET_EDGES(tail, head, edge_id, {
    const Vertex mp_index = tail - 1;
    const Vertex edm_index = head - BIPARTITE - 1;

    degrees[mp_index] += 1.0;

    for (int k = 0; k < K; ++k) {
      topic_mass[(size_t) mp_index * (size_t) K + (size_t) k] +=
        INPUT_PARAM[(size_t) edm_index * (size_t) K + (size_t) k];
    }
  });
}

/* Keep private storage synchronized whenever a toggle is applied. */
U_CHANGESTAT_FN(u_b1topicdiversity) {
  const int K = IINPUT_PARAM[0];
  const Vertex n1 = BIPARTITE;
  const Vertex mp_index = tail - 1;
  const Vertex edm_index = head - BIPARTITE - 1;
  const int sign = edgestate ? -1 : 1;

  GET_STORAGE(double, state);

  double *degrees = state;
  double *topic_mass = state + n1;

  degrees[mp_index] += (double) sign;

  for (int k = 0; k < K; ++k) {
    topic_mass[(size_t) mp_index * (size_t) K + (size_t) k] +=
      (double) sign *
      INPUT_PARAM[(size_t) edm_index * (size_t) K + (size_t) k];
  }
}

/* Return the change in the sum of mode-1 inverse Simpson topic counts. */
C_CHANGESTAT_FN(c_b1topicdiversity) {
  const int K = IINPUT_PARAM[0];
  const int subtract_one = IINPUT_PARAM[1];
  const Vertex n1 = BIPARTITE;
  const Vertex mp_index = tail - 1;
  const Vertex edm_index = head - BIPARTITE - 1;
  const int sign = edgestate ? -1 : 1;

  GET_STORAGE(double, state);

  const double *degrees = state;
  const double *all_topic_mass = state + n1;
  const double *mp_topic_mass =
    all_topic_mass + (size_t) mp_index * (size_t) K;
  const double *edm_topics =
    INPUT_PARAM + (size_t) edm_index * (size_t) K;

  const double old_value = inverse_simpson_topic_count(
    mp_topic_mass,
    degrees[mp_index],
    K,
    subtract_one
  );

  const double new_value = inverse_simpson_topic_count_after_toggle(
    mp_topic_mass,
    degrees[mp_index],
    edm_topics,
    K,
    sign,
    subtract_one
  );

  CHANGE_STAT[0] += new_value - old_value;
}
