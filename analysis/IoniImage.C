/** \file
 * Convert simulated energy deposition data into ionization images.
 *
 * Generated images are written into _p1.root
 *
 * Optionally an MRC file containing charge cloud can be written.
 */

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TH3I.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TUnuran.h>
#include <TUnuranMultiContDist.h>
#include <TMath.h>

extern "C" {
#include "mrcIo.h"
}

#if defined(__MAKECINT__) || defined(__CINT__)
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<double> >+;
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<long> >+;
#endif

#ifndef MIN
#define MIN(a,b) ((a)>(b)?(b):(a))
#endif

/** Cholesky decomposition.  Adapted from NR3.
 * @param[in] a input array, must be square, symmetric and positive-definite, with n*n elements.
 * @param[in] n dimension of array a and L.
 * @param[out] L decomposition result, a lower triangular matrix.
 *               If (*L)=NULL is supplied, *L is allocated.
 * @return 1 when success, 0 if a is not positive-definite.
 */
int cholesky_decomp(const double *a, ssize_t n, double **L)
{
    ssize_t i, j, k;
    double sum, *el;
    if(*L == NULL) {
        (*L) = (double*)calloc(n*n, sizeof(double));
        if(*L == NULL) {
            fprintf(stderr, "%s(): *L allocation failure.\n", __FUNCTION__);
            return 0;
        }
    }
    el = *L;
    for(i=0; i<n*n; i++) {el[i] = a[i];}
    for(i=0; i<n; i++) {
        for(j=i; j<n; j++) {
            for(sum=el[i*n+j],k=i-1; k>=0; k--) sum -= el[i*n+k] * el[j*n+k];
            if(i == j) {
                if(sum <= 0.0) {
                    fprintf(stderr, "%s(): a is not positive-definite.\n", __FUNCTION__);
                    return 0;
                }
                el[i*n+i] = sqrt(sum);
            } else {
                el[j*n+i] = sum/el[i*n+i];
            }
        }
    }
    for(i=0; i<n; i++)
        for(j=0; j<i; j++)
            el[j*n+i] = 0.0;
    return 1;
}
/** n-dimmensional Gaussian (multivariate normal) distribution.
 * @param[in] n dimension.
 * @param[in] L Cholesky decomposed covariance matrix.
 * @param[in] mean vector of mean values.
 * @param[out] pt n-dimmensional Gaussian deviate output.  pt has to be pre-allocated.
 */
double *rand_gaussnd(TRandom *tr, ssize_t n, const double *L, const double *mean, double *pt)
{
#ifndef RAND_GAUSSND_NMAX
#define RAND_GAUSSND_NMAX 1024
#endif /* for efficiency */
    double spt[RAND_GAUSSND_NMAX];

    ssize_t i, j;
    double u, v, x, y, q;
    for(i=0; i<n; i++) { /* fill vector spt of independent normal deviates. */
        do{
            u = tr->Uniform();
            v = 1.7156*(tr->Uniform()-0.5);
            x = u - 0.449871;
            y = fabs(v) + 0.386595;
            q = x*x + y*(0.19600*y-0.25472*x);
        } while (q > 0.27597 && (q > 0.27846 || v*v > -4.*log(u)*u*u));
        spt[i] = v/u;
    }
    for(i=0; i<n; i++) {
        pt[i] = 0.0;
        for(j=0; j<=i; j++) pt[i] += L[i*n+j] * spt[j];
        pt[i] += mean[i];
    }
    return pt;
#undef RAND_GAUSSND_NMAX
}

typedef std::vector< std::vector<float> > vecvec;
// for xenon at 10bar
double Fano=0.14;
double Wi=24.8; // eV
double Dt=1.0; // transverse diffusion (sigma) [mm]
double Dl=2.0; // longitudinal diffusion (sigma) [mm]
int nbin=400, bmin=-200, bmax=200; // binning

int IoniImage(TTree *t1, const char *ofname, int pProj=0, int omrc=-1)
{
    ssize_t i, j, k;

    mrc_header_t *mh=0;
    float *mdata=0;
    
    TFile *tfp=0;
    TTree *tp1=0;

    std::vector<int> *parentId=0;
    std::vector<double> *xp=0, *yp=0, *zp=0, *ed=0; // must be initialized to NULL
    int nIonTot;
    TH3I *cloudH;
    int nIon;
    double mIon, sIon;
    double cov[9] = {Dt*Dt,   0.0,   0.0,
                       0.0, Dt*Dt,   0.0,
                       0.0,   0.0, Dl*Dl};
    double L[9], *Lp, sxyz[3], mean[3]={0.0, 0.0, 0.0};
    TRandom3 *tr = new TRandom3;
    /*
    TUnuran *tur = new TUnuran(tr);
    TF3 *gs3 = new TF3("gs3", "(2.0*pi*[0]*[1]*[2])**(-3.0/2.0)*exp(-0.5*(x**2/[0]**2+y**2/[1]**2+z**2/[2]**2))",
                       -DBL_MAX,DBL_MAX, -DBL_MAX,DBL_MAX, -DBL_MAX,DBL_MAX);
    gs3->SetParameters(Dt, Dt, Dl);
    TUnuranMultiContDist dist(gs3);
    tur->Init(dist);
    */
    Lp = L;
    cholesky_decomp(cov, 3, &Lp);

    cloudH = new TH3I("cloudH", "Charge cloud",
                      nbin, bmin, bmax,
                      nbin, bmin, bmax,
                      nbin, bmin, bmax);

    
    t1->SetBranchAddress("parentid", &parentId);
    t1->SetBranchAddress("ed", &ed);
    t1->SetBranchAddress("xp", &xp);
    t1->SetBranchAddress("yp", &yp);
    t1->SetBranchAddress("zp", &zp);    

    if(omrc<0) {
        tfp = new TFile(ofname, "RECREATE");
        tfp->cd();
        tp1 = new TTree("p1", "IoniImage");
        tp1->Branch("nIonTot", &nIonTot, "nIonTot/I");
        tp1->Branch("cloudH", "TH3I", &cloudH);
    } else {
        mh = mrcIo_fill_header(NULL, nbin, nbin, nbin);
        mdata = (float*)calloc(nbin*nbin*nbin, sizeof(float));
    }

    for(i=0; i<t1->GetEntries(); i++) {
        if(omrc>=0) {
            i = MIN(omrc, t1->GetEntries());
        }
        
        t1->GetEntry(i);
        cloudH->Reset();
        nIonTot = 0;
        for(j=0; j<(ssize_t)ed->size(); j++) {
            mIon = (*ed)[j] * 1000.0 / Wi;
            sIon = std::sqrt(Fano * mIon);
            nIon = (int)tr->Gaus(mIon, sIon);
            if(nIon<0) nIon = 0;
            nIonTot += nIon;
            if(nIon == 0) continue;
            //handle diffusion
            for(k=0; k<nIon; k++) {
                // tur->SampleMulti(sxyz);
                rand_gaussnd(tr, 3, L, mean, sxyz);
                cloudH->Fill((*xp)[j] + sxyz[0],
                             (*yp)[j] + sxyz[1],
                             (*zp)[j] + sxyz[2]);
            }
        }
        if(omrc<0) {
            tp1->Fill();
        } else {
            break;
        }
    }

    if(omrc<0) {
        tp1->Write();
        tfp->CurrentFile()->Close();
        delete tfp;
    } else {
        for(i=0; i<nbin; i++) {
            for(j=0; j<nbin; j++) {
                for(k=0; k<nbin; k++) {
                    *(mdata + k*nbin*nbin + j*nbin + i)
                        = cloudH->GetBinContent(cloudH->GetBin(i+1, j+1, k+1));
                }
            }
        }
        mrcIo_write_file(mh, ofname, mdata);
        
        free(mdata);
        free(mh);
    }
    
    delete cloudH;
    //delete gs3;
    //delete tur;
    delete tr;

    return i;
}

#ifndef __CINT__
/** Main entry if this file is compiled outside of root */
int main(int argc, char **argv)
{
    char *sRootFname;
    std::string ofname;
    std::stringstream ss;
    int pProj=0, omrc=-1;
    ssize_t sz;

    TFile *tfs;
    TTree *t1s;
    
    if(argc > 1) {
        sRootFname = argv[1];
        if(argc > 2) pProj = atoi(argv[2]);
        if(argc > 3) omrc = atoi(argv[3]);
        if(argc > 4) Fano = atof(argv[4]);
        if(argc > 5) Wi = atof(argv[5]);
        
        ofname.assign(sRootFname);
        sz = ofname.size();
        if(omrc<0) {
            ofname.resize(sz+3);
            ofname.replace(sz-5, 8, "_p1.root");
        } else {
            ofname.resize(sz-5);
            ss << omrc;
            ofname += ss.str();
            ofname += ".mrc";
        }
        std::cout << "Output file: " << ofname << std::endl;
    } else {
        std::cerr << argv[0] << " SimulatedEdRootFname [pProj] [omrc] [Fano] [Wi]\n"
                  << "            [pProj] : plane of projection: 0-xy; 1-yz; 2-xz\n"
                  << "            [omrc]  : output to MRC: (<0)-no output; 0,1,2...-eventId\n"
                  << "            [Fano]  : Fano factor, default for Xe: 0.14\n"
                  << "            [Wi]    : W value, default for Xe: 24.8eV"
                  << std::endl;
        return EXIT_FAILURE;
    }

    // signal
    tfs = new TFile(sRootFname, "READ");
    std::cout << "******************** Signal ************************" << std::endl;
    tfs->ls();
    t1s = (TTree*)(tfs->Get("t1"));
    t1s->ls();
    std::cout << "t1 has " << t1s->GetEntries() << " entries." << std::endl;

    IoniImage(t1s, ofname.c_str(), pProj, omrc);

    tfs->Close();
    delete tfs;

    return EXIT_SUCCESS;
}
#endif
