/**
 * Author: Isabel Rosa, isarosa@mit.edu, 
 * Jay Hilton, jhilton@mit.edu, 
 * Krit Boonsiriseth, talkon@mit.edu
 **/

#include <assert.h>
#include <cilk/cilk.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../common/render.h"
#include "../include/misc_utils.h"

typedef struct renderer_state {
  renderer_spec_t r_spec;
  float *img;
} renderer_state_t;

struct renderer_state* init_renderer(const renderer_spec_t *spec) {
  renderer_state_t *state = (renderer_state_t*)malloc(sizeof(renderer_state_t));
  state->r_spec = *spec;
  int n_pixels = state->r_spec.resolution * state->r_spec.resolution;
  state->img = calloc(3ull * (size_t)n_pixels, sizeof(float));
  assert(state->img != NULL);
  return (struct renderer_state*)state;
}

void destroy_renderer(struct renderer_state *state) {
  renderer_state_t *rs = (renderer_state_t*)state;
  free(rs->img);
  free(rs);
}

// Computes the ray from a given origin (usually the eye location) to the pixel (x, y)
// in image coordinates.
static ray_t origin_to_pixel(renderer_state_t *state, int x, int y) {
  ray_t viewingRay;

  const float pixel_size = state->r_spec.viewport_size / state->r_spec.resolution;
  
  // Center image frame
  float us = -state->r_spec.resolution / 2.0f + x;
  float vs = -state->r_spec.resolution / 2.0f + y;

  viewingRay.origin = state->r_spec.eye;
  viewingRay.dir =
    qsubtract(
      qadd(
        scale(us * pixel_size, state->r_spec.proj_plane_u), 
        scale(vs * pixel_size, state->r_spec.proj_plane_v)
      ),
      viewingRay.origin
    );
  viewingRay.dir = scale(1 / qsize(viewingRay.dir), viewingRay.dir);

  return viewingRay;
}

// Determines whether the ray r and the sphere s intersect.
// 
// If the ray and the sphere intersect, writes the distance to the closer intersection
// to `out`, and returns 1. Otherwise, returns 0.
int ray_sphere_intersection(ray_t *r, const sphere_t* s, float *out) {
  vector_t dist = qsubtract(r->origin, s->pos);
  // Uses quadratic formula to compute intersection
  float a = qdot(r->dir, r->dir);
  float b = 2 * qdot(r->dir, dist);
  float c = (float)((double)qdot(dist, dist) - (double)(s->r * s->r));
  float discr = (float)((double)(b * b) - (double)(4 * a * c));

  if (discr >= 0) {
    // Ray hits sphere
    float sqrtdiscr = sqrtf(discr);

    float min_dist;
    if (b >= 0) {
      float sol1 = (float)((double)-b - (double)sqrtdiscr) / (2 * a);
      float sol2 = (float)(2 * (double)c) / ((double)-b - (double)sqrtdiscr);
      min_dist = min(sol1, sol2);
    } else {
      float sol1 = (float)(2 * (double)c) / ((double)-b + (double)sqrtdiscr);
      float sol2 = (float)((double)-b + (double)sqrtdiscr) / (2 * a);
      min_dist = min(sol1, sol2);
    }

    // If new_t > 0 and smaller than original t, we
    // found a new, closer ray-sphere intersection
    if (min_dist > 0) {
      *out = min_dist;
      return 1;
    }
  }

  return 0;
}

// Sorts given spheres by length of the tangent (NOT in-place). 
// Returns a pointer to the sorted spheres. Returned pointer must be freed.
//
// Since the spheres are non-intersecting, this ensures that 
// if sphere S comes before sphere T in this ordering, then 
// sphere S is in front of sphere T in the rendering.


void set_pixel(renderer_state_t *state, int x, int y, float red, float green, float blue) {
    renderer_state_t *rs = (renderer_state_t*)state;
    int index = (x + y * rs->r_spec.resolution) * 3;
    state->img[index + 0] = min((float)red, 1.0);
    state->img[index + 1] = min((float)green, 1.0);
    state->img[index + 2] = min((float)blue, 1.0);
}

static inline vector_t vcross(vector_t a, vector_t b) {
  vector_t c;
  c.x = a.y * b.z - a.z * b.y;
  c.y = a.z * b.x - a.x * b.z;
  c.z = a.x * b.y - a.y * b.x;
  return c;
}

typedef struct {
  vector_t eye;
  vector_t u, v;
  vector_t u_hat, v_hat, n;
  float Lu, Lv;
  float pixel_size;
  int resolution;
} proj_ctx_t;

static inline int project_point_to_pixel(const proj_ctx_t *ctx, vector_t P, float *out_x, float *out_y) {
  vector_t EP = qsubtract(P, ctx->eye);
  float denom = qdot(ctx->n, EP);
  const float eps = 1e-8f;
  if (fabsf(denom) < eps) return 0; 
  float t = -qdot(ctx->n, ctx->eye) / denom;
  if (t <= 0.0f) return 0; 
  vector_t Q = qadd(ctx->eye, scale(t, EP));
  float alpha = qdot(Q, ctx->u_hat);
  float beta  = qdot(Q, ctx->v_hat);
  float us = alpha / (ctx->pixel_size * ctx->Lu);
  float vs = beta  / (ctx->pixel_size * ctx->Lv);
  float x = us + ctx->resolution / 2.0f;
  float y = vs + ctx->resolution / 2.0f;
  *out_x = x;
  *out_y = y;
  return 1;
}

static vector_t g_eye_for_sort;

static int sphere_depth_cmp(const void *pa, const void *pb) {
  const sphere_t *a = (const sphere_t*)pa;
  const sphere_t *b = (const sphere_t*)pb;
  vector_t da = qsubtract(a->pos, g_eye_for_sort);
  vector_t db = qsubtract(b->pos, g_eye_for_sort);
  float da2 = qdot(da, da);
  float db2 = qdot(db, db);
  if (da2 < db2) return -1;
  if (da2 > db2) return 1;
  return 0;
}



const float* render(renderer_state_t *state, const sphere_t *spheres, int n_spheres) {
    renderer_state_t *rs = (renderer_state_t*)state;
    const renderer_spec_t *spec = &rs->r_spec;
    const int resolution = spec->resolution;
    const int n_pixels = resolution * resolution;

    memset(rs->img, 0, sizeof(float) * 3 * (size_t)n_pixels);
    float *z_buffer = (float*)malloc(sizeof(float) * (size_t)n_pixels);
    if (!z_buffer) return rs->img;
    for (int i = 0; i < n_pixels; i++) {
        z_buffer[i] = INFINITY;
    }

    sphere_t *sorted_spheres = (sphere_t*)malloc(sizeof(sphere_t) * (size_t)n_spheres);
    if (!sorted_spheres) { free(z_buffer); return rs->img; }
    memcpy(sorted_spheres, spheres, sizeof(sphere_t) * (size_t)n_spheres);
    g_eye_for_sort = spec->eye;
    qsort(sorted_spheres, (size_t)n_spheres, sizeof(sphere_t), sphere_depth_cmp);
  


      double red = 0;
      double green = 0;
      double blue = 0;

      for (int j = 0; j < state->r_spec.n_lights; j++) {
        light_t currentLight = state->r_spec.lights[j];
        vector_t intersection_to_light = qsubtract(currentLight.pos, intersection);
        if (qdot(normal, intersection_to_light) <= 0)
          continue;

        ray_t lightRay;
        lightRay.origin = intersection;
        lightRay.dir = scale(1 / qsize(intersection_to_light), intersection_to_light);

        // Calculate Lambert diffusion
        float lambert = qdot(lightRay.dir, normal);
        red += (double)(currentLight.intensity.red * currentMat.diffuse.red *
                        lambert);
        green += (double)(currentLight.intensity.green *
                          currentMat.diffuse.green * lambert);
        blue += (double)(currentLight.intensity.blue * currentMat.diffuse.blue *
                         lambert);
      }

      set_pixel(state, x, y, red, green, blue);
      
  free(sorted_spheres);
  return state->img;
}
