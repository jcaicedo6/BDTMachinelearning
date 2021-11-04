#include <iostream>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>
#include <TLeaf.h>
#include <TMath.h>
#include "TMVA/Reader.h"
#include <TStopwatch.h>
#include <vector>


using namespace TMVA;



using namespace RooFit;


void Data_filter(TChain *chain, RooDataSet *data, RooRealVar x);
//void GetSandB(DataMC *t, Double_t &S, Double_t &B, Double_t bdt);
void PlotMass();
//void PlotSignificance();
void SetReader();

TGraph *gr;

TMVA::Reader *reader;

Int_t bins = 50;
Double_t Ymin_min = 1.701;
Double_t Ymin_max = 1.798;
void TMVABDTApplyAnalysisMmax()
{
    
    //PlotSignificance();
    PlotMass();
    
    return;
}




void PlotMass()
{
    //Adding data in a TChain
    TChain *chData_signal = new TChain("tf");
    chData_signal->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/signal_Mmin_c.root");
    TTree *taudata = (TTree*) chData_signal;

    TChain *chData_taupair = new TChain("tf");
    chData_taupair->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/bkg_Mmin_c.root");
    TTree *taupair = (TTree*) chData_taupair;
    
    TChain *chData_bkg = new TChain("tf");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/uub_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/ddb_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/ccb_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/ssb_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/charged_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/mixed_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/ee_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/eeee1_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/eeee2_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/eeee3_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/eemumu1_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/eemumu1_1_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/eemumu2_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/mumu1_Mmin_c.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_min/mumu2_Mmin_c.root");
    TTree *treeB = (TTree*) chData_bkg;

    //TChain *chData_bkg1 = new TChain("tf");
    //chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/ee_c.root");
    //chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eeee1_c.root");
    //chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eeee2_c.root");
    //chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eeee3_c.root");
    //chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eemumu1_c.root");
    //chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eemumu1_1_c.root");
    //chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eemumu2_c.root");
    ///chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/mumu1_c.root");
    //chData_bkg1->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/mumu2_c.root");
    //TTree *treeB1 = (TTree*) chData_bkg1;
    

    //TChain *chData_all = new TChain("tf");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/signal_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/bkg_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/uub_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/ddb_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/ccb_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/ssb_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/charged_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/mixed_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/ee_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eeee1_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eeee2_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eeee3_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eemumu1_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eemumu1_1_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/eemumu2_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/mumu1_c.root");
    //chData_all->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/Tau1prong/data_cuts/m_edge/mumu2_c.root");
    //TTree *treeA = (TTree*) chData_all;
    

    RooRealVar Mtau("Mtau", "M_{tau}", Ymin_min, Ymin_max, "GeV/c^{2}");
    RooDataSet *dataMtau = new RooDataSet("dataMtau","dataMtau",RooArgSet(Mtau));
    RooDataSet *dataMtauP = new RooDataSet("dataMtauP","dataMtauP",RooArgSet(Mtau));
    RooDataSet *dataMtauB = new RooDataSet("dataMtauB","dataMtauB",RooArgSet(Mtau));
    RooDataSet *dataMtauB1 = new RooDataSet("dataMtauB1","dataMtauB1",RooArgSet(Mtau));
    RooDataSet *dataMtauA = new RooDataSet("dataMtauA","dataMtauA",RooArgSet(Mtau));
    
    
    RooCategory tagCat("tagCat","tagging Category");
    tagCat.defineType("signal");
    tagCat.defineType("taupair");
    tagCat.defineType("bkg");
    //tagCat.defineType("bkg1");
    //tagCat.defineType("all");

    
    
    tagCat.setLabel("signal");
    Data_filter(chData_signal,dataMtau,Mtau);
    dataMtau->addColumn(tagCat);
    
    
    
    tagCat.setLabel("taupair");
    Data_filter(chData_taupair,dataMtauP,Mtau);
    dataMtauP->addColumn(tagCat);
    
    
    tagCat.setLabel("bkg");
    Data_filter(chData_bkg,dataMtauB,Mtau);
    dataMtauB->addColumn(tagCat);

    //tagCat.setLabel("bkg1");
    //Data_filter(chData_bkg1,dataMtauB1,Mtau);
    //dataMtauB1->addColumn(tagCat);
    
    
    dataMtau->append(*dataMtauP);
    dataMtau->append(*dataMtauB);
    //dataMtau->append(*dataMtauB1);

    //RooRealVar P11("P11","P11", 1.778, 1.76, 1.8);
    //RooRealVar P21("P21","P21", 0.01, -1, 1);
    //RooRealVar P31("P31","P31", -0.42/3.15, -100, 100);
    //RooRealVar P41("P41","P41", -0.10/3.15, -50.0, 50.0);
    //RooRealVar P51("P51","P51", -1.0/3.15, -50, 50);
    //RooRealVar P61("P61","P61", 0, -100, 100);
    //RooRealVar P71("P71","P71", 0, -100, 100);
    //RooRealVar P81("P81","P81", 0, -100, 100);
    //RooRealVar P91("P91","P91", 0, -100, 100);
    
    RooRealVar P11("P11","P11", 1.778, 0, 1000.0);
    RooRealVar P21("P21","P21", 0.01, -1000, 1000);
    RooRealVar P31("P31","P31", -0.42/3.15, -100, 100);
    RooRealVar P41("P41","P41", -0.10/3.15, -50.0, 50.0);
    RooRealVar P51("P51","P51", -1.0/3.15, -50, 50);
    RooRealVar P61("P61","P61", 0, -100, 100);
    RooRealVar P71("P71","P71", 0, -100, 100);
    RooRealVar P81("P81","P81", 0, -100, 100);



    //RooGenericPdf MassModelY("PsMassPDFY", "(@3 + @4*@0 + @6*@0*@0 + @7*@0*@0*@0 + @8*@0*@0*@0*@0) * atan((@0 - @1)/@2) + @5 * @0 + 1",
                              //RooArgSet(Mtau,P11,P21,P31,P41,P51,P61,P71,P81));

    //RooGenericPdf MassModelY("PsMassPDFY2", "(@1 + @6*@0 + @7*@0*@0 + @8*@0*@0*@0 + @9@0*@0*@0*@0) * atan((@0 - @3)/@2) + @4*@0 + @5",
                              //RooArgSet(Mtau,P11,P21,P31,P41,P51,P61,P71,P81));
    //RooGenericPdf MassModelY("PsMassPDFY", "(@3 + @4*@0 + @6*@0*@0 + @7*@0*@0*@0) * TMath::Erfc((@0 - @1)/@2) + @8*@0*@0 + @5 * @0 + 1",
                            // RooArgSet(Mtau,P11,P21,P31,P41,P51,P61,P71,P81));

    //RooFitResult *fitmasY = MassModelY.fitTo(*dataMtau, Save(kTRUE), Strategy(1), NumCPU(2));
    //fitmasY->Print("v");

    TCanvas *cYmin = new TCanvas("cMtau","cMtau",800,800);
    RooPlot* thrFrameY = Mtau.frame(Ymin_min,Ymin_max,bins);
    dataMtau->plotOn(thrFrameY, Name("Mtau"), DataError(RooAbsData::SumW2));
    dataMtau->plotOn(thrFrameY, Name("Mtau"), DataError(RooAbsData::SumW2),MarkerColor(kBlue),Cut("tagCat==tagCat::signal"));
    dataMtau->plotOn(thrFrameY, Name("Mtau"), DataError(RooAbsData::SumW2),MarkerColor(kGreen),Cut("tagCat==tagCat::bkg"));
    //dataMtau->plotOn(thrFrameY, Name("Mtau"), DataError(RooAbsData::SumW2),MarkerColor(kYellow),Cut("tagCat==tagCat::bkg1"));
    dataMtau->plotOn(thrFrameY, Name("Mtau"), DataError(RooAbsData::SumW2),MarkerColor(kRed),Cut("tagCat==tagCat::taupair"));
    //dataMtau->plotOn(thrFrameY, Name("Mtau"), DataError(RooAbsData::SumW2),MarkerColor(kBlack),Cut("tagCat==tagCat::all"));
    //MassModelY.plotOn(thrFrameY, Name("model"));
    //MassModelY.paramOn(thrFrameY,Layout(0.6,0.9,0.9));
    thrFrameY->GetXaxis()->SetTitleSize(0);
    thrFrameY->Draw();
    
}

void Data_filter(TChain *chain, RooDataSet *data, RooRealVar x)
{
    Double_t xmin,xmax;
    xmin = Ymin_min;
    xmax = Ymin_max;
    Long64_t nentries = chain->GetEntries();
    cout<<"No. events in data: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    
    //Bool_t IsSignal = false;
    
    for (Long64_t jentry=0; jentry!=nentries;jentry++)
    {
        Long64_t ientry = chain->LoadTree(jentry);
        if (ientry < 0) break;
        nb = chain->GetEntry(jentry);
        nbytes += nb;
        if(jentry%(nentries/10)==0) cout<<(100.0*jentry)/nentries<<" %"<<endl;

        if(chain->GetLeaf("iddata")->GetValue() == 2) continue;// IsSignal=true;
        //else IsSignal=false;

        Double_t val = chain->GetLeaf("BDT")->GetValue();

        if(val<0.24) continue;
        //if(IsSignal==false) continue;

        Double_t xval = 0;
        xval = TMath::Sqrt(chain->GetLeaf("Mmax")->GetValue());

        if(xval > xmin && xval < xmax)
        {
            x = xval;
            data->add(x);
        }

       
    }


    return;
}


