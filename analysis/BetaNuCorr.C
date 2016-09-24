/** \file
 * Analyze beta-nu correlation simulation output.
 *
 * This file processes the .root file with simulated trajectory data
 * and computes the additional information useful for Beta-Neutrino
 * correlation experiments.  Additional information is saved in `t2'.
 *
 * run with
 * root -l -b -q '../analysis/BetaNuCorr.C("events.root")'
 *
 * Beta-Neutrino momentum correlation can be drawn by
 * t1->AddFriend("t2");
 * t1->Draw("(xm[p1idx[0]]*xm[p1idx[1]]+ym[p1idx[0]]*ym[p1idx[1]]+zm[p1idx[0]]*zm[p1idx[1]])");
 */

#include <stdio.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>

using namespace std;

void BetaNuCorr(string &fileName)
{
    TFile *tf;
    TTree *t1, *t2;

    vector<double> *time;
    vector<int> *parentId;
    vector<double> *xp, *yp, *zp;
    vector<int> *p1idx;
    vector<double> *xm, *ym, *zm;

    double v[3];
    int i;

    gROOT->Reset();

    gROOT->SetStyle("Plain");
    gROOT->ProcessLine("#include <vector>");

    tf = new TFile(fileName.c_str(), "UPDATE");
    tf->ls();

    t2 = new TTree("t2", "Derived quantities");
    xm = new vector<double>;
    ym = new vector<double>;
    zm = new vector<double>;
    p1idx = new vector<int>;

    t2->Branch("xm", "vector<double>", &xm);
    t2->Branch("ym", "vector<double>", &ym);
    t2->Branch("zm", "vector<double>", &zm);
    t2->Branch("p1idx", "vector<int>", &p1idx);

    t1 = (TTree*)tf->Get("t1");
    t1->ls();

    t1->SetBranchAddress("time", &time);
    t1->SetBranchAddress("parentid", &parentId);
    t1->SetBranchAddress("xp", &xp);
    t1->SetBranchAddress("yp", &yp);
    t1->SetBranchAddress("zp", &zp);    

    for(int j=0; j<t1->GetEntries(); j++) {
        t1->GetEntry(j);

        xm->clear(); ym->clear(); zm->clear();
        p1idx->clear();
        for(int k=0; k<time->size(); k++) {
            i = (int)(*time)[k];
            // cout << i << ": " << (*xp)[i] << ", " << (*yp)[i] << ", " << (*zp)[i] << endl;
            v[0] = (*xp)[i+1] - (*xp)[i];
            v[1] = (*yp)[i+1] - (*yp)[i];
            v[2] = (*zp)[i+1] - (*zp)[i];
            TMath::Normalize(v);
            // cout << i << ": " << v[0] << ", " << v[1] << ", " << v[2] << endl;

            xm->push_back(v[0]);
            ym->push_back(v[1]);
            zm->push_back(v[2]);

            if((*parentId)[k] == 1) {
                p1idx->push_back(k);
            }
        }
        t2->Fill();
    }

    tf->Write();
    tf->Close();
    delete tf;
    // delete t2; // already freed by tf->Close() ?
    delete xm; delete ym; delete zm;
    delete p1idx;
//    return EXIT_SUCCESS;
}
