/**
 * Author: Isabel Rosa, isrosa@mit.edu
 **/

#include <assert.h>
#include <cilk/cilk.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "../../common/simulate.h"
#include "../include/misc_utils.h"

typedef struct simulator_state {
  simulator_spec_t s_spec;
  sphere_t *spheres;
} simulator_state_t;

simulator_state_t* init_simulator(const simulator_spec_t *spec) {
  simulator_state_t *state = (simulator_state_t*)malloc(sizeof(simulator_state_t));
  state->s_spec = *spec;
  state->spheres = malloc(2 * spec->n_spheres * sizeof(sphere_t));
  assert(state->spheres != NULL);
  memcpy(state->spheres, spec->spheres, sizeof(sphere_t) * spec->n_spheres);
  memcpy(state->spheres + state->s_spec.n_spheres, spec->spheres, sizeof(sphere_t) * state->s_spec.n_spheres);
  return state;
}

void destroy_simulator(simulator_state_t* state) {
  free(state->spheres);
  free(state);
}

void update_accel_sphere(sphere_t *spheres, int n_spheres, double g, int i) {
  double rx = 0;
  double ry = 0;
  double rz = 0;

  const vector_t zero_vec = {0, 0, 0};
  spheres[i + n_spheres].accel = zero_vec;
  for (int j = 0; j < n_spheres; j++) {
    if (i != j) {
      vector_t i_minus_j = qsubtract(spheres[i].pos, spheres[j].pos);
      vector_t j_minus_i = scale(-1, i_minus_j);
      vector_t force =
          scale(g * spheres[j].mass / pow(qsize(i_minus_j), 3), j_minus_i);
      rx += (double)force.x;
      ry += (double)force.y;
      rz += (double)force.z;
    }
  }
  const vector_t v = {.x = rx, .y = ry, .z = rz};
  spheres[i + n_spheres].accel = v;
}

static inline void add_interaction(sphere_t *spheres, int n_spheres, double g, int i, int j) {
  // r_ij = pos_j - pos_i
  vector_t r_ij = qsubtract(spheres[j].pos, spheres[i].pos);
  double dist_sq = (double)qdot(r_ij, r_ij);
  if (dist_sq == 0.0) return; 
  double inv_dist_cubed = 1.0 / (sqrt(dist_sq) * dist_sq); // 1/|r|^3

  float coeff_i = (float)(g * spheres[j].mass * inv_dist_cubed);
  float coeff_j = (float)(-g * spheres[i].mass * inv_dist_cubed);

  spheres[i + n_spheres].accel =
      qadd(spheres[i + n_spheres].accel, scale(coeff_i, r_ij));
  spheres[j + n_spheres].accel =
      qadd(spheres[j + n_spheres].accel, scale(coeff_j, r_ij));
}

static void serial_triangle(sphere_t *spheres, int n_spheres, double g, int lo, int hi) {
  for (int i = lo; i < hi; ++i) {
    for (int j = i + 1; j < hi; ++j) {
      add_interaction(spheres, n_spheres, g, i, j);
    }
  }
}

static void serial_rect(sphere_t *spheres, int n_spheres, double g,
                        int i0, int i1, int j0, int j1) {
  for (int i = i0; i < i1; ++i) {
    for (int j = j0; j < j1; ++j) {
      add_interaction(spheres, n_spheres, g, i, j);
    }
  }
}

void update_accelerations(sphere_t *spheres, int n_spheres, double g) {

  cilk_for (int i = 0; i < n_spheres; ++i) {
    update_accel_sphere(spheres, n_spheres, g, i);
  }
}

void update_velocities(sphere_t *spheres, int n_spheres, float t) {
  cilk_for (int i = 0; i < n_spheres; i++) {
    spheres[i + n_spheres].vel = qadd(spheres[i].vel, scale(t, spheres[i].accel));
  }
}

void update_positions(sphere_t *spheres, int n_spheres, float t) {
  cilk_for (int i = 0; i < n_spheres; i++) {
    spheres[i + n_spheres].pos = qadd(spheres[i].pos, scale(t, spheres[i].vel));
  }
}

// runs simulation for minCollisionTime timesteps
// perform collision between spheres at indices i and j
void do_ministep(sphere_t *spheres, int n_spheres, double g, float minCollisionTime, int i, int j) {
  update_accelerations(spheres, n_spheres, g);
  update_velocities(spheres, n_spheres, minCollisionTime);
  update_positions(spheres, n_spheres, minCollisionTime);

  cilk_for (int k = 0; k < n_spheres; k++) {
    spheres[k] = spheres[k + n_spheres];
  }

  if (i == -1 || j == -1) {
    return;
  }

  vector_t distVec = qsubtract(spheres[i].pos, spheres[j].pos);
  float scale1 = 2 * spheres[j].mass /
                 (float)((double)spheres[i].mass + (double)spheres[j].mass);
  float scale2 = 2 * spheres[i].mass /
                 (float)((double)spheres[i].mass + (double)spheres[j].mass);
  float distNorm = qdot(distVec, distVec);
  vector_t velDiff = qsubtract(spheres[i].vel, spheres[j].vel);
  vector_t scaledDist = scale(qdot(velDiff, distVec) / distNorm, distVec);
  spheres[i].vel = qsubtract(spheres[i].vel, scale(scale1, scaledDist));
  spheres[j].vel = qsubtract(spheres[j].vel, scale(-1 * scale2, scaledDist));
}

// Check if the spheres at indices i and j collide in the next
// timeToCollision timesteps
// 
// If so, modifies timeToCollision to be the time until spheres i and j collide.
int check_for_collision(sphere_t *spheres, int i, int j, float *timeToCollision) {
  vector_t distVec = qsubtract(spheres[i].pos, spheres[j].pos);
  float dist = qsize(distVec);
  float sumRadii = (float)((double)spheres[i].r + (double)spheres[j].r);

  // Shift frame of reference to act like sphere i is stationary
  // Not adjusting for acceleration because our simulation does not adjust for acceleration
  vector_t movevec = qsubtract(spheres[j].vel, spheres[i].vel);

  // Distance that sphere j moves in timeToCollision time
  float moveDist = (float)((double)qsize(movevec) * (double)*timeToCollision);

  // Break if the length the sphere moves in timeToCollision time is less than
  // distance between the centers of these spheres minus their radii
  if ((double)moveDist < (double)dist - (double)sumRadii ||
      (movevec.x == 0 && movevec.y == 0 && movevec.z == 0)) {
    return 0;
  }

  vector_t unitMovevec = scale(1 / qsize(movevec), movevec);

  // distAlongMovevec = ||distVec|| * cos(angle between unitMovevec and distVec)
  float distAlongMovevec = qdot(unitMovevec, distVec);

  // Check that sphere j is moving towards sphere i
  if (distAlongMovevec <= 0) {
    return 0;
  }

  float jToMovevecDistSq =
      (float)((double)(dist * dist) -
              (double)(distAlongMovevec * distAlongMovevec));

  // Break if the closest that sphere j will get to sphere i is more than
  // the sum of their radii
  float sumRadiiSquared = sumRadii * sumRadii;
  if (jToMovevecDistSq >= sumRadiiSquared) {
    return 0;
  }

  // We now have jToMovevecDistSq and sumRadii, two sides of a right triangle.
  // Use these to find the third side, sqrt(T)
  float extraDist = (float)((double)sumRadiiSquared - (double)jToMovevecDistSq);

  if (extraDist < 0) {
    return 0;
  }

  // Draw out the spheres to check why this is the distance sphere j moves
  // before hitting sphere i;)
  float distance = (float)((double)distAlongMovevec - (double)sqrt(extraDist));

  // Break if the distance sphere j has to move to touch sphere i is too big
  if (distance < 0 || moveDist < distance) {
    return 0;
  }

  *timeToCollision = distance / qsize(movevec);
  return 1;
}

void do_timestep(simulator_state_t* state, float timeStep) {
  float timeLeft = timeStep;

  // If collisions are getting too frequent, we cut time step early
  // This allows for smoother rendering without losing accuracy
  while (timeLeft > 0.000001) {
    float minCollisionTime = timeLeft;
    int indexCollider1 = -1;
    int indexCollider2 = -1;

    cilk_for (int i = 0; i < state->s_spec.n_spheres; i++) {
      cilk_for (int j = i + 1; j < state->s_spec.n_spheres; j++) {
        if (check_for_collision(state->spheres, i, j, &minCollisionTime)) {
          indexCollider1 = i;
          indexCollider2 = j;
        }
      }
    }

    do_ministep(state->spheres, state->s_spec.n_spheres, state->s_spec.g, minCollisionTime, indexCollider1, indexCollider2);

    timeLeft = timeLeft - minCollisionTime;
  }
}

sphere_t* simulate(simulator_state_t* state) {
  int n_spheres = state->s_spec.n_spheres;
  float timeStep = n_spheres > 1 ? (1 / log(n_spheres)) : 1;
  do_timestep(state, timeStep);
  return state->spheres;
}
