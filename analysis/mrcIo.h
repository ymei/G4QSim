#ifndef __MRCIO_H__
#define __MRCIO_H__

/** MRC file header format */
typedef struct
{
    int32_t nx;   /* number of points per dimension for the whole dataset */
    int32_t ny;
    int32_t nz;
    int32_t mode; /* should =2, 32-bit float point format */
    int32_t nxstart; /* image idx (can be a sub region in the data), */
    int32_t nystart; /* runs from n?start to n?start + m? - 1 */
    int32_t nzstart;
    int32_t mx; /* number of intervals (sub region size) */
    int32_t my; /* better make n? = m? */
    int32_t mz;
    float xlen; /* cell size (Angstroms), pixel spacing = ?len/m? */
    float ylen;
    float zlen;
    float alpha; /* cell angles (Degrees), ignore */
    float beta;
    float gamma;
    int32_t mapc; /* map column 1=x, 2=y, 3=z, set to 1 */
    int32_t mapr; /* map row, =2 */
    int32_t maps; /* map section, =3 */
    float dmin; /* pixel value min/max and mean, for proper scaling */
    float dmax;
    float dmean;
    int32_t ispg; /* space group number, ignore */
    int32_t nsymbt; /* number of bytes in extended header (for symmetry operators), 0 or 80 */
    int32_t extra[25]; /* all set to 0 */
    float xorig; /* origin of image */
    float yorig;
    float zorig;
    char map[4]; /* should contain "MAP " */
    char machst[4]; /* machine stamp, [0] = 0x11, [1] = 0x11 for big-endian */
                    /* [0] = 0x44, [1] = 0x41 for little-endian */
    float rms; /* rms deviation of map density from mean */
    int32_t nlabl; /* number of labels with useful data */
    char labels[10][80];
} mrc_header_t;

mrc_header_t *mrcIo_fill_header(mrc_header_t *mh, int32_t nx, int32_t ny, int32_t nz);
int mrcIo_write_file(mrc_header_t *mh, const char *fname, const float *data);

#endif /* __MRCIO_H__ */
