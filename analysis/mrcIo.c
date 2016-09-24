/** \file
 * MRC file writing utilities.
 */
#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mrcIo.h"

/** Initialize and fill an MRC file header. */
mrc_header_t *mrcIo_fill_header(mrc_header_t *mh, int32_t nx, int32_t ny, int32_t nz)
{
    if(sizeof(mrc_header_t) != 1024) {
        fprintf(stderr, "Error, mrc_header_t size is %zd instead of 1024.\n", sizeof(mrc_header_t));
        return NULL;
    }
    if(mh == NULL) {
        mh = (mrc_header_t*)malloc(sizeof(mrc_header_t));
        if(mh == NULL) {
            fprintf(stderr, "Error mallocing mrc_header_t.\n");
            return NULL;
        }
    }
    bzero(mh, sizeof(mrc_header_t));
    mh->nx = nx;
    mh->ny = ny;
    mh->nz = nz;
    mh->mode = 2;
    mh->nxstart = -nx/2;
    mh->nystart = -ny/2;
    mh->nzstart = -nz/2;
    mh->mx = nx;
    mh->my = ny;
    mh->mz = nz;
    mh->xlen = (float)nx;
    mh->ylen = (float)ny;
    mh->zlen = (float)nz;
    mh->mapc = 1;
    mh->mapr = 2;
    mh->maps = 3;
    mh->xorig = 0.0;
    mh->yorig = 0.0;
    mh->zorig = 0.0;
    mh->map[0] = 'M'; mh->map[1] = 'A'; mh->map[2] = 'P'; mh->map[3] = ' ';
    mh->machst[0] = 0x44; mh->machst[1] = 0x41;
    return mh;
}

/** Write data into an MRC file. */
int mrcIo_write_file(mrc_header_t *mh, const char *fname, const float *data)
{
    FILE *fp;
    size_t i, n;
    int32_t nx, ny, nz;
    float min, max, val;
    double mean, rms;

    if(mh == NULL) {
        fprintf(stderr, "Error mrc_header_t *mh == NULL.\n");
        return 1;
    }
    if((fp = fopen(fname, "wb")) == NULL) {
        fprintf(stderr, "Error opening file %s for writing.\n", fname);
        return 1;
    }
    
    /* calculate min, max, mean, rms */
    nx = mh->nx; ny = mh->ny; nz = mh->nz;
    min = FLT_MAX; max = FLT_MIN;
    mean = 0.0;
    rms = 0.0;
    n = (size_t)nx * (size_t)ny * (size_t)nz;
    for(i=0; i<n; i++) {
        val = data[i];
        if(val < min) min = val;
        if(val > max) max = val;
        mean += val;
        rms += val * val;
    }
    rms = sqrt(fabs(rms - mean * mean / (double)n)/(n - 1.0));
    mean /= (double)n;
    
    mh->dmin = min;
    mh->dmax = max;
    mh->dmean = mean;
    mh->rms = rms;

    fwrite(mh, sizeof(mrc_header_t), 1, fp);
    fwrite(data, sizeof(float) * n, 1, fp);
    fclose(fp);
    return 0;
}

