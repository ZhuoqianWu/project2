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

// int ray_sphere_intersection(ray_t *r, const sphere_t* s, float *out) {
//     vector_t dist = qsubtract(r->origin, s->pos);
//     float a = qdot(r->dir, r->dir);
//     float b = 2.0f * qdot(r->dir, dist);
//     float c = (float)((double)qdot(dist, dist) - (double)(s->r * s->r));
    
//     // 计算判别式
//     float discr = (float)((double)(b * b) - double(4 * a * c));
//     if (discr < 0) {
//         return 0; // 无相交
//     }
    
//     float sqrtdiscr = sqrtf(discr);
//     float q;
//     if (b < 0) {
//         q = -0.5f * (b - sqrtdiscr);
//     } else {
//         q = -0.5f * (b + sqrtdiscr);
//     }
    
//     // 计算两个可能的解
//     float t1 = q / a;
//     float t2 = c / q;
    
//     // 确保t1 <= t2
//     if (t1 > t2) {
//         float temp = t1;
//         t1 = t2;
//         t2 = temp;
//     }
    
//     // 寻找最小正解
//     float min_dist = t1 > 0 ? t1 : t2;
//     if (min_dist <= 0) {
//         return 0; // 解无效（负值或零）
//     }
    
//     *out = min_dist;
//     return 1;
// }

// Sorts given spheres by length of the tangent (NOT in-place). 
// Returns a pointer to the sorted spheres. Returned pointer must be freed.
//
// Since the spheres are non-intersecting, this ensures that 
// if sphere S comes before sphere T in this ordering, then 
// sphere S is in front of sphere T in the rendering.


void set_pixel(struct renderer_state *state, int x, int y, float red, float green, float blue) {
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



const float* render(struct renderer_state *state, const sphere_t *spheres, int n_spheres) {
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

    vector_t temp_u = spec->proj_plane_u;
    vector_t temp_v = spec->proj_plane_v;
    if (qsize(temp_u) == 0.0f || qsize(temp_v) == 0.0f) {
        free(z_buffer);
        return rs->img;
    }
    vector_t cam_dir = vcross(temp_u, temp_v);
    cam_dir = scale(1.0f/qsize(cam_dir), cam_dir);
    vector_t cam_u = temp_u;
    cam_u = scale(1.0f/qsize(cam_u), cam_u);
    vector_t cam_v = vcross(cam_dir, cam_u);
    cam_v = scale(1.0f/qsize(cam_v), cam_v);

    const float viewport_height = spec->viewport_size;
    const float fov_factor = viewport_height / (2.0f * (float)resolution);

    sphere_t *sorted_spheres = (sphere_t*)malloc(sizeof(sphere_t) * (size_t)n_spheres);
    if (!sorted_spheres) { free(z_buffer); return rs->img; }
    memcpy(sorted_spheres, spheres, sizeof(sphere_t) * (size_t)n_spheres);
    g_eye_for_sort = spec->eye;
    qsort(sorted_spheres, (size_t)n_spheres, sizeof(sphere_t), sphere_depth_cmp);

    proj_ctx_t ctx;
    ctx.eye = spec->eye;
    ctx.u = spec->proj_plane_u;
    ctx.v = spec->proj_plane_v;
    ctx.Lu = qsize(ctx.u);
    ctx.Lv = qsize(ctx.v);
    if (ctx.Lu == 0.0f || ctx.Lv == 0.0f) {
        free(z_buffer); free(sorted_spheres); return rs->img;
    }
    ctx.u_hat = scale(1.0f / ctx.Lu, ctx.u);
    ctx.v_hat = scale(1.0f / ctx.Lv, ctx.v);
    ctx.n = vcross(ctx.u, ctx.v);
    float Ln = qsize(ctx.n);
    if (Ln == 0.0f) {
        free(z_buffer); free(sorted_spheres); return rs->img;
    }
    ctx.n = scale(1.0f / Ln, ctx.n);
    ctx.pixel_size = spec->viewport_size / (float)resolution;
    ctx.resolution = resolution;

    for (int i = 0; i < n_spheres; i++) {
        const sphere_t *s = &sorted_spheres[i];
        float xs[13], ys[13]; 
        int count = 0;
        vector_t c = s->pos;
        vector_t p[13];
        p[0] = c;
        p[1] = qadd(c, scale(s->r, ctx.u_hat));
        p[2] = qsubtract(c, scale(s->r, ctx.u_hat));
        p[3] = qadd(c, scale(s->r, ctx.v_hat));
        p[4] = qsubtract(c, scale(s->r, ctx.v_hat));
        vector_t diag1_dir = qadd(ctx.u_hat, ctx.v_hat);
        float diag1_len = qsize(diag1_dir);
        if (diag1_len > 0) {
            vector_t diag1_norm = scale(1.0f / diag1_len, diag1_dir);
            p[5] = qadd(c, scale(s->r, diag1_norm));
            p[6] = qsubtract(c, scale(s->r, diag1_norm));
        } else {
            p[5] = p[1]; p[6] = p[2];
        }
        
        vector_t diag2_dir = qsubtract(ctx.u_hat, ctx.v_hat);
        float diag2_len = qsize(diag2_dir);
        if (diag2_len > 0) {
            vector_t diag2_norm = scale(1.0f / diag2_len, diag2_dir);
            p[7] = qadd(c, scale(s->r, diag2_norm));
            p[8] = qsubtract(c, scale(s->r, diag2_norm));
        } else {
            p[7] = p[1]; p[8] = p[2];
        }
            p[9] = qadd(c, scale(s->r, ctx.n));
    p[10] = qsubtract(c, scale(s->r, ctx.n));
    
    vector_t extra_dir = qadd(ctx.u_hat, ctx.n);
    float extra_len = qsize(extra_dir);
    if (extra_len > 0) {
        vector_t extra_norm = scale(1.0f / extra_len, extra_dir);
        p[11] = qadd(c, scale(s->r, extra_norm));
        p[12] = qsubtract(c, scale(s->r, extra_norm));
    } else {
        p[11] = p[1]; p[12] = p[2];
    }
        for (int k = 0; k < 5; k++) {
            float xx, yy;
            if (project_point_to_pixel(&ctx, p[k], &xx, &yy)) {
                xs[count] = xx;
                ys[count] = yy;
                count++;
            }
        }
        if (count == 0) continue;

        float min_x = xs[0], max_x = xs[0];
        float min_y = ys[0], max_y = ys[0];
        for (int k = 1; k < count; k++) {
            if (xs[k] < min_x) min_x = xs[k];
            if (xs[k] > max_x) max_x = xs[k];
            if (ys[k] < min_y) min_y = ys[k];
            if (ys[k] > max_y) max_y = ys[k];
        }
        
        float distance_to_camera = qsize(qsubtract(s->pos, ctx.eye));
        float projected_radius = (s->r / fmaxf(distance_to_camera, 1e-8f)) * resolution;
      
        int scene_density_factor = (n_spheres > 200) ? 2 : 1;
        const int base_padding = (int)fmaxf(5 * scene_density_factor, 
                                        fminf(10 * scene_density_factor, 
                                                projected_radius * 0.15f * scene_density_factor));
      
        float width = max_x - min_x;
        float height = max_y - min_y;
        min_x -= width * 0.04f;
        max_x += width * 0.04f;
        min_y -= height * 0.04f;
        max_y += height * 0.04f;
        int x0 = (int)floorf(min_x) - padding;
        int x1 = (int)ceilf(max_x) + padding;
        int y0 = (int)floorf(min_y) - padding;
        int y1 = (int)ceilf(max_y) + padding;
        if (x1 < 0 || x0 >= resolution || y1 < 0 || y0 >= resolution) {
            bboxes[i].valid = 0;
            continue;
        }

        x0 = max(x0, 0);
        y0 = max(y0, 0);
        x1 = min(x1, resolution - 1);
        y1 = min(y1, resolution - 1);

        bboxes[i].x0 = x0;
        bboxes[i].y0 = y0;
        bboxes[i].x1 = x1;
        bboxes[i].y1 = y1;
        bboxes[i].valid = 1;
    }

        for (int y = y0; y <= y1; y++) {
            for (int x = x0; x <= x1; x++) {
                ray_t r = origin_to_pixel(rs, x, y);
                float t;
                if (!ray_sphere_intersection(&r, s, &t)) continue;
                int idx = x + y * resolution;
                if (t >= z_buffer[idx]) continue;
                z_buffer[idx] = t;

                material_t currentMat = s->mat;
                vector_t intersection = qadd(r.origin, scale(t, r.dir));
                vector_t normal = qsubtract(intersection, s->pos);
                float n_size = qsize(normal);
                if (n_size == 0.0f) continue;
                normal = scale(1.0f / n_size, normal);

                double red = 0.0, green = 0.0, blue = 0.0;
                for (int j = 0; j < spec->n_lights; j++) {
                    light_t currentLight = spec->lights[j];
                    vector_t intersection_to_light = qsubtract(currentLight.pos, intersection);
                    if (qdot(normal, intersection_to_light) <= 0.0f) continue;
                    ray_t lightRay;
                    lightRay.origin = intersection;
                    float itl_size = qsize(intersection_to_light);
                    if (itl_size == 0.0f) continue;
                    lightRay.dir = scale(1.0f / itl_size, intersection_to_light);
                    float lambert = qdot(lightRay.dir, normal);
                    red   += (double)(currentLight.intensity.red   * currentMat.diffuse.red   * lambert);
                    green += (double)(currentLight.intensity.green * currentMat.diffuse.green * lambert);
                    blue  += (double)(currentLight.intensity.blue  * currentMat.diffuse.blue  * lambert);
                }
                set_pixel(state, x, y, (float)red, (float)green, (float)blue);
            }
        }
    }
  free(z_buffer);
  free(sorted_spheres);
  return rs->img;
}

