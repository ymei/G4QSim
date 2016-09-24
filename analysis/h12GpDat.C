/** \file
 * Save a TH1 histogram's contents into a text .dat file for plotting in gnuplot.
 *
 * In root, use .L or .x to load or run this.
 */
#include <stdio.h>
#include <stdlib.h>
#include <TH1.h>

int h12GpDat(TH1 *h1, const char *fname)
{
    FILE *fp;
    Int_t i, NX;
    TAxis *XAxis;
    Double_t x,v,err;

    if((fp = fopen(fname, "w"))==NULL) {
        perror(fname);
        return 0;
    }
    NX = h1->GetNbinsX();
    XAxis = h1->GetXaxis();
    for(i=1; i<=NX; i++) {
        x   = XAxis->GetBinLowEdge(i);
        v   = h1->GetBinContent(i);
        err = h1->GetBinError(i);
        fprintf(fp, "%g %24.16e %d %24.16e\n", x, v, i, err);
    }
    fclose(fp);
    return 1;
}
