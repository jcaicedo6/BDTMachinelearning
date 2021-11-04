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

using namespace std;

//Here we declare global variables


Bool_t calculateFOMS(TLorentzVector mom_plus, TLorentzVector mom_minus,  Double_t Ecms,Double_t &Ymin, Double_t &Ymax)
{
  Bool_t accepted = true;
  Double_t CMS_E = Ecms;
  
  TLorentzVector pa = mom_plus;
  TLorentzVector pb = mom_minus;
  
  TVector3 n_a = (1.0/CMS_E)*pa.Vect();
  TVector3 n_b = (1.0/CMS_E)*pb.Vect();
  Double_t ab = n_a.Dot(n_b);
  
  Double_t za = pa.E()/CMS_E;
  Double_t zb = pb.E()/CMS_E;
  
  Double_t a2 = n_a.Mag2();
  Double_t b2 = n_b.Mag2();
  
  Double_t wa = (zb*zb - zb - b2 - 2*ab);
  Double_t wb = (za*za - za + a2);
  
  TVector3 H;
  H = wa*n_a + wb*n_b; 
  
  TVector3 acrossb = n_a.Cross(n_b);
  
  Double_t A1 = b2;
  Double_t A2 = a2;
  Double_t A3 = 2.0*ab;
  Double_t B1 = 2.0*(n_b.Dot(H));
  Double_t B2 = 2.0*(n_a.Dot(H));
  Double_t C1 = 4.0*acrossb.Mag2();
  Double_t D1 = H.Dot(H) -C1*(0.5 - za)*(0.5 - za);
  
  Double_t A = A1 + A2 + A3;
  Double_t B = B1 + B2;
  Double_t C = D1;
  
  
  if( (B*B - 4.0*A*C)<0 )
    {
      Ymin=-1000;
      Ymax=-1000;
      accepted = false;
    }
  else
    {
      Ymin = CMS_E*CMS_E*(-B - TMath::Sqrt(B*B - 4*A*C))/(2*A);
      Ymax = CMS_E*CMS_E*(-B + TMath::Sqrt(B*B - 4*A*C))/(2*A);
      accepted = true;
    }
    
  return accepted;
}



void produceTree(TChain *chain,vector<TString> leafs, TString filename, Int_t ident)
{

  Int_t nleafs = (Int_t)leafs.size();
  //Int_t nleafs = 16;
  TFile *rootfile = new TFile(filename,"RECREATE");
  TTree *treef = new TTree("tf","Friend tree");
  
  Float_t bdtval,Mmin,Mmax;
  Int_t idp;    
  
  treef->Branch("BDT", &bdtval, "BDT/F");
  treef->Branch("Mmin", &Mmin, "Mmin/F");
  treef->Branch("Mmax", &Mmax, "Mmax/F");
  treef->Branch("id", &idp, "id/I");
  
  //create the Reader object
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  Float_t var[15];
  for(Int_t i=0;i<nleafs;i++)
    {
      reader->AddVariable(leafs[i],&var[i]);
    }
  
  // Book methods
  reader->BookMVA("_BDT", TString("dataset/weights/_BDT.weights.xml"));
  
  //print(branchesTMVA)
  
  for(Int_t k=0;k<chain->GetEntries();k++)
    {
      chain->GetEntry(k);
      
      //Here we assign the values to the BDT variables
      //Note: I tried to use the same variables in the chain but somehow the direction get lost sometimes
      for(Int_t i=0;i<nleafs;i++)
        {
	  Double_t fval = chain->GetLeaf(leafs[i])->GetValue();
	  var[i] = fval;
        }
      
      Double_t fbdt = reader->EvaluateMVA("_BDT");
      bdtval=fbdt;
      
      //Here we calculate mass edge variables
      TLorentzVector pa = TLorentzVector();
      TLorentzVector pb = TLorentzVector();
      
      pa.SetPxPyPzE(chain->GetLeaf("track1__px_CMS")->GetValue(),
		    chain->GetLeaf("track1__py_CMS")->GetValue(),
		    chain->GetLeaf("track1__pz_CMS")->GetValue(),
		    chain->GetLeaf("track1__E_CMS")->GetValue());
      pb.SetPxPyPzE(chain->GetLeaf("track2__px_CMS")->GetValue(),
		    chain->GetLeaf("track2__py_CMS")->GetValue(),
		    chain->GetLeaf("track2__pz_CMS")->GetValue(),
		    chain->GetLeaf("track2__E_CMS")->GetValue());
      Double_t Ecms = chain->GetLeaf("Ecms")->GetValue();
      if(Ecms==0) Ecms = 10.58;
      Double_t Ymin,Ymax;
      Bool_t flagFOMS = calculateFOMS(pa,pb,Ecms,Ymin,Ymax);
      Mmin = Ymin;
      Mmax = Ymax;
      //By default is flagMOMS is false, Ymin and Ymax are equal to -1000
      
      idp = ident;
      
      treef->Fill();
    }        
  rootfile->cd();
  treef->Write();
  rootfile->Close();
}

int main(int argc, char **argv)
{
  char* inputfile = argv[1];
  char* outputfile = argv[2];
  Int_t id = (Int_t)atoi(argv[3]);
  
  TString infile = TString(inputfile);
  TString outfile = TString(outputfile);


  cout<<" Input file: "<<infile<<"  "<<endl;
  cout<<" Ouput file: "<<outfile<<"  "<<endl;
  cout<<" Type : "<<id<<endl;
  
  //TString outfile = "signal_c.root";

  TChain *chdata = new TChain("tau1x1");
  chdata->AddFile(inputfile);
  
  TString leafsall[16] ={"track1__cosToThrustOfEvent","track2__cosToThrustOfEvent","track1__clusterE",
    "track2__clusterE","visibleEnergyOfEventCMS","missingMomentumOfEventCMS",
    "missingMomentumOfEventCMS_theta","missingMass2OfEvent","track1__EoverP",
    "track2__EoverP","track1__pt","track2__pt","track1__pionID","track2__pionID",
    "track1__p","track2__p"};
  
  vector<TString> leafs;
  for(Int_t i=0;i<16;i++)
    {
      TString lf = leafsall[i];
      //leafs[i]=lf;
      leafs.push_back(lf);
    }
  
  TStopwatch t;
  t.Start();
  produceTree(chdata,leafs,outfile,id);
  t.Stop();
  t.Print(); 
  return 0;
}
