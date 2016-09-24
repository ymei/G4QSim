/** \file
 * Make images out of processed data, optionally generate an MRC file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <TROOT.h>
#include <TFile.h>
#include <TParameter.h>
#include <TTree.h>
#include <TH1.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TColor.h>
#include <cfloat>
#include <string>
#include <sstream>
#include <vector>
#if defined(__MAKECINT__) || defined(__CINT__)
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<double> >+;
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<long> >+;
#endif

extern "C" {
#include "mrcIo.h"
}

#ifndef MIN
#define MIN(a,b) ((a)>(b)?(b):(a))
#endif

#define NX 72
size_t nCh = (72*72);
double dt = 0.663552; // [ms]

typedef std::vector< std::vector<float> > vecvec;

TFile *tfp;
TTree *tp1;
vecvec *wavf=0;
TH1D *coinH=0;
std::vector<double> *img=0, *imgMdF2x2=0, *imgMdF3x3=0;
TCanvas *mainCanvas;
TProfile2D *imgH;
TGaxis *x2Axis;

/** Setup a TCanvas for plotting */
TCanvas *setup_plots()
{
    gROOT->ProcessLine(".L ~/.rootstyle.C");
    gROOT->ProcessLine("RootStyle()");
    gROOT->SetStyle("Root-Style");
    gStyle->SetTitleAlign(13);
    gStyle->SetTitleW(0.5);
    gStyle->SetTitleX(0.0);
    
    mainCanvas = new TCanvas("c1", "", 1, 1, 460, 480);
    imgH = new TProfile2D("imgH", "", 72, 0.0, 72.0, 72, 0.0, 72.0);
    imgH->GetXaxis()->SetTitle("x [pix]");
    imgH->GetYaxis()->SetTitle("y [pix]");
    imgH->GetZaxis()->SetTitle("[mV]");
    imgH->GetZaxis()->SetTitleOffset(0.9);
    imgH->GetZaxis()->SetRangeUser(-5,15);
    imgH->Sumw2(true);
    imgH->Draw();
    x2Axis = new TGaxis(0,72,72,72,0,6,10,"-C");
    x2Axis->SetLabelFont(42);
    x2Axis->SetLabelSize(0.035);
    x2Axis->SetTitleFont(42);
    x2Axis->SetTitleSize(0.04);
    x2Axis->SetTitle("[mm]");
    x2Axis->Draw();
    return mainCanvas;
}

/** Make images or an MRC file.
 *  This function can be called directly in a root session
 *  after loading with .L
 */
int ImagePlots(const char *fname, ssize_t iStart=-1, ssize_t iStop=-1, int omrc=0)
{
    double val;
    ssize_t n, i, j, ii, jj, ix, iy, xx, yy;
    std::vector<double> vlist(9);
    std::ostringstream ss0;

    tfp = new TFile(fname, "READ");
    tfp->ls();
    coinH = (TH1D*)tfp->Get("coinH");
    tp1 = (TTree*)(tfp->Get("p1"));
    tp1->ls();
    std::cout << "p1 has " << tp1->GetEntries() << " entries." << std::endl;
    tp1->SetBranchAddress("wavf", &wavf);
    img       = new std::vector<double>(nCh);
    imgMdF2x2 = new std::vector<double>(nCh);
    imgMdF3x3 = new std::vector<double>(nCh);

    tp1->GetEntry(0);
    std::cout << (*wavf)[0].size() << " frames total" << std::endl;
    if(iStop < 0 || iStop > (*wavf)[0].size()) iStop = (*wavf)[0].size();
    if(iStart < 0) iStart = 0;
    if(iStart >= iStop) iStart = iStop-1;
    std::cout << "iStart = " << iStart << ", iStop = " << iStop << std::endl;

    mrc_header_t *mh;
    float *data;
    if(omrc==1) {
        n = iStop - iStart;        
        mh = mrcIo_fill_header(NULL, NX, nCh/NX, n);
        data = new float[nCh * n];        
    }

    setup_plots();
    for(i=iStart; i<iStop; i++) {
        printf("\n\n# i = %zd\n", i);
        for(j=0; j<nCh; j++) {
            (*img)[j] = (*wavf)[j][i];
        }
        imgH->Reset();
        for(iy=0; iy<nCh/NX; iy++) {
            for(ix=0; ix<NX; ix++) {
                /* 2x2 median filter the image */
                n = 0;
                for(jj=0; jj<=1; jj++) {
                    for(ii=0; ii<=1; ii++) {
                        yy = iy + jj; if(yy>=nCh/NX) yy -= 2;
                        xx = ix + ii; if(xx>=NX) xx -=2;
                        val = (*img)[yy * NX + xx];
                        vlist[n] = val;
                        n++;
                    }
                }
                std::nth_element(vlist.begin(), vlist.begin() + n/2, vlist.begin() + n);
                (*imgMdF2x2)[iy * NX + ix] = vlist[n/2];
                
                /* 3x3 median filter the image */                
                n = 0;
                for(jj=-1; jj<=1; jj++) {
                    for(ii=-1; ii<=1; ii++) {
                        yy = iy + jj; if(yy>=nCh/NX) yy -= 2; if(yy<0) yy += 2;
                        xx = ix + ii; if(xx>=NX) xx -= 2; if(xx<0) xx += 2;
                        val = (*img)[yy * NX + xx];
                        vlist[n] = val;
                        n++;
                    }
                }
                std::nth_element(vlist.begin(), vlist.begin() + n/2, vlist.begin() + n);
                (*imgMdF3x3)[iy * NX + ix] = vlist[n/2];
                // Fill
                imgH->Fill(ix, iy, 1000.0 * (*img)[iy * NX + ix]);
                if(omrc==1) {
                    data[(i-iStart) * nCh + iy * NX + ix] = (*imgMdF3x3)[iy * NX + ix];
                } else if(omrc>1) {
                    printf("%zd %zd %g\n", ix, iy, 1000.0 * (*img)[iy * NX + ix]);
                }
            }
        }
        ss0.clear(); ss0.str("");
        ss0 << "t = " << i * dt << " [ms]";
        imgH->SetTitle(ss0.str().c_str());
        imgH->Draw("colz");
        x2Axis->Draw();
        mainCanvas->Modified();
        mainCanvas->Update();
        ss0.clear(); ss0.str("");
        ss0 << "out" << std::setfill('0') << std::setw(7) << i << ".png";
        if(!omrc) {
            // mainCanvas->Print(ss0.str().c_str(), "png");
            mainCanvas->Print("out.gif+1", "gif");
        }
        // break;
    }
    printf("\n");
        
    if(omrc==1) {
        mh->zlen *= 1.0;
        mrcIo_write_file(mh, "out.mrc", data);
        delete []data;
        free(mh);
    }
    delete img;
    delete imgMdF2x2;
    delete imgMdF3x3;
    
    return 0;
}

#ifndef __CINT__
/** Main entry if this file is compiled outside of root */
int main(int argc, char **argv)
{
    ssize_t iStart=-1, iStop=-1;
    int omrc=0;
    if(argc<2) {
        fprintf(stderr, "%s rootfname [iStart=-1] [iStop=-1] [omrc=0]\n", argv[0]);
        return EXIT_FAILURE;
    }
    if(argc>2) iStart = atol(argv[2]);
    if(argc>3) iStop  = atol(argv[3]);
    if(argc>4) omrc   = atoi(argv[4]);
    
    ImagePlots(argv[1], iStart, iStop, omrc);

    tfp->CurrentFile()->Close();
    return EXIT_SUCCESS;
}
#endif
