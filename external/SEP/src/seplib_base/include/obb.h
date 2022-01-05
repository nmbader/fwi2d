#ifndef _oriented_bounding_box
#define _oriented_bounding_box
extern
#ifdef __cplusplus
"C"
#endif
/* takes X, Y coordinates and returns minimal oriented bounding box */
void oriented_bounding_box_f(int npts, const float *x, const float *y,
                           float xbb[4], float ybb[4]);
extern
#ifdef __cplusplus
"C"
#endif
void oriented_bounding_box_d(int npts, const double *x, const double *y,
                           double xbb[4], double ybb[4]);
#endif/*_oriented_bounding_box*/
