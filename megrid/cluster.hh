/******************************************************************************/
/* The C Clustering Library.
 * Copyright (C) 2002 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Laboratory of DNA Information Analysis,
 * Human Genome Center, Institute of Medical Science, University of Tokyo,
 * 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
 * Contact: mdehoon 'AT' gsc.riken.jp
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

#ifdef WINDOWS
#  include <windows.h>
#endif

#include "types.hh"

#define CLUSTERVERSION "1.50"

/* Chapter 2 */
real_t clusterdistance (int nrows, int ncolumns, real_t** data, int** mask,
  real_t weight[], int n1, int n2, int index1[], int index2[], char dist,
  char method, int transpose);
real_t** distancematrix (int ngenes, int ndata, real_t** data,
  int** mask, real_t* weight, char dist, int transpose);

/* Chapter 3 */
int getclustercentroids(int nclusters, int nrows, int ncolumns,
  real_t** data, int** mask, int clusterid[], real_t** cdata, int** cmask,
  int transpose, char method);
void getclustermedoids(int nclusters, int nelements, real_t** distance,
  int clusterid[], int centroids[], real_t errors[]);
void kcluster (int nclusters, int ngenes, int ndata, real_t** data,
  int** mask, real_t weight[], int transpose, int npass, char method, char dist,
  int clusterid[], real_t* error, int* ifound);
void kmedoids (int nclusters, int nelements, real_t** distance,
  int npass, int clusterid[], real_t* error, int* ifound);

/* Chapter 4 */
typedef struct {int left; int right; real_t distance;} Node;
/*
 * A Node struct describes a single node in a tree created by hierarchical
 * clustering. The tree can be represented by an array of n Node structs,
 * where n is the number of elements minus one. The integers left and right
 * in each Node struct refer to the two elements or subnodes that are joined
 * in this node. The original elements are numbered 0..nelements-1, and the
 * nodes -1..-(nelements-1). For each node, distance contains the distance
 * between the two subnodes that were joined.
 */

Node* treecluster (int nrows, int ncolumns, real_t** data, int** mask,
  real_t weight[], int transpose, char dist, char method, real_t** distmatrix);
void cuttree (int nelements, Node* tree, int nclusters, int clusterid[]);

/* Chapter 5 */
void somcluster (int nrows, int ncolumns, real_t** data, int** mask,
  const real_t weight[], int transpose, int nxnodes, int nynodes,
  real_t inittau, int niter, char dist, real_t*** celldata,
  int clusterid[][2]);

/* Chapter 6 */
int pca(int m, int n, real_t** u, real_t** v, real_t* w);

/* Utility routines, currently undocumented */
void sort(int n, const real_t data[], int index[]);
real_t mean(int n, real_t x[]);
real_t median (int n, real_t x[]);

real_t* calculate_weights(int nrows, int ncolumns, real_t** data, int** mask,
  real_t weights[], int transpose, char dist, real_t cutoff, real_t exponent);
