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


Bool_t calculateFOMS(TLorentzVector mom_plus, TLorentzVector mom_minus, Double_t Ecms,Double_t &M_tau_Edge)
{
  Bool_t accepted = true;
  Double_t CMS_E = Ecms;
  
  TLorentzVector pa = mom_plus;
  TLorentzVector pb = mom_minus;

  //tay mass compute
  Double_t px1,py1,pz1,px2,py2,pz2,E1,E2,Ecm;
  Double_t z_pi1, z_pi2, a_pi1x, a_pi1y, a_pi1z, a_pi2x, a_pi2y, a_pi2z;
  Double_t punto_a_b,A_0,B_0,C_0,D_0,M_nu_Edge,M_tau_Max;

  //Normalized energies (pion neutrino)
  z_pi1 = pa.E()/CMS_E;
  z_pi2 = pb.E()/CMS_E;
  //normalized moments in components
  a_pi1x = pa.Px()/CMS_E;
  a_pi1y = pa.Py()/CMS_E;
  a_pi1z = pa.Pz()/CMS_E;
  a_pi2x = pb.Px()/CMS_E;
  a_pi2y = pb.Py()/CMS_E;
  a_pi2z = pb.Pz()/CMS_E;
  // create the moments vectors
  TVector3 a_pi1, a_pi2, sum_a_b,cruz_a_b,comp1,comp2;
  a_pi1.SetXYZ(a_pi1x,a_pi1y,a_pi1z);
  a_pi2.SetXYZ(a_pi2x,a_pi2y,a_pi2z);

  // define some variables previous for add vectors, scalar and vectorial products
       
  sum_a_b = a_pi1+a_pi2;
  punto_a_b = a_pi1.Dot(a_pi2);
  cruz_a_b = a_pi1.Cross(a_pi2);
  comp1 = ((z_pi2*z_pi2) - z_pi2)*a_pi1 + ((z_pi1*z_pi1) - z_pi1)*a_pi2;
        
  // Define the parameters A0, B0, C0 y D0
  A_0 = sum_a_b.Mag2();//magnitud al cuadrado
  B_0 = (2.0*a_pi1.Mag2())*((z_pi2*z_pi2) - z_pi2) + (2.0*a_pi2.Mag2())*((z_pi1*z_pi1) - z_pi1) + (2.0*punto_a_b)*((z_pi1*z_pi1) + (z_pi2*z_pi2) - z_pi1 - z_pi2 - (sum_a_b.Mag2()));
  C_0 = 4.0*cruz_a_b.Mag2();
  D_0 = a_pi1.Mag2()*a_pi2.Mag2()*sum_a_b.Mag2() - cruz_a_b.Mag2() +(2.0*a_pi1.Mag2()*a_pi2.Mag2())*(z_pi1 + z_pi2 - (z_pi1*z_pi1) - (z_pi2*z_pi2)) + comp1.Mag2() - (2.0*punto_a_b)*(a_pi1.Mag2()*((z_pi2*z_pi2) - z_pi2) + a_pi2.Mag2()*((z_pi1*z_pi1) - z_pi1));
  // limit mass for the tau and the neutrino
  Double_t x = (4.0*(B_0*B_0) + 3.0*(C_0*C_0) - 16.0*(A_0*D_0) - 8.0*(B_0*C_0))/(16.0*(A_0*C_0));
  Double_t y = (4.0*(B_0*B_0) - (C_0*C_0) - 16.0*(A_0*D_0))/(16*(A_0*C_0));
  //Double_t mu_tau2_Edge2 = ((TMath::Sqrt((B_0*B_0) - 4.0*(A_0*D_0))) - B_0)/(2.0*(A_0));
        
  //Double_t mu_tau_Max2 = ((B_0 - C_0)*(B_0 - C_0))/(4.0*(A_0*C_0)) - ((D_0)/(C_0));
  Double_t mu_tau2_Edge2 = ((TMath::Sqrt((B_0*B_0) - 4.0*(A_0*D_0))) - B_0)/(2.0*(A_0));

  // Now we check if the vertex is in quadrant IV (y<0)
  if (y >= 0)
  {
     M_tau_Edge = TMath::Sqrt(x)*CMS_E; 
  }

  else 
  {
     M_tau_Edge = M_tau_Edge = TMath::Sqrt(mu_tau2_Edge2)*CMS_E;
     M_nu_Edge = 0;
  }

  //if (y < 0)
  //{
    //  Double_t rval = B_0*B_0 - 4*A_0*D_0;
      //if (rval < 0)
      //{
        //accepted = false;
        //x = +1000;
        //y = +1000;

     // }
      //else
      //{
        //x = (TMath::Sqrt(rval) - B_0)/(2*A_0);
        //y = 0;
      //}
  //}

  // Now lets check that x>0
  //if(x < 0)
  //{
    //M_tau_Edge = +1000;
    //accepted = false;
  //}
  //else 
  //{
    //M_tau_Edge = TMath::Sqrt(CMS_E*CMS_E*x);
  //}
  return accepted;
  
}


void produceTree(TChain *chain,vector<TString> leafs, TString filename, Int_t ident)
{

  Int_t nleafs = (Int_t)leafs.size();
  //Int_t nleafs = 16;
  TFile *rootfile = new TFile(filename,"RECREATE");
  TTree *treef = new TTree("tf","Friend tree");
  
  Float_t bdtval,Mtau;
  Int_t idp,iddata;
  
  //let's create a new id where id=1 is bdt, id=2 is bias and id=3 is data    
  
  treef->Branch("BDT", &bdtval, "BDT/F");
  treef->Branch("Mtau", &Mtau, "Mtau/F");
  treef->Branch("id", &idp, "id/I");
  treef->Branch("iddata", &iddata, "iddata/I");
  
 
  
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
  Long64_t nentries, nentries_20, nentries_30, nentries_50;
  nentries = chain->GetEntries();
  nentries_20 = 0.20*nentries;
  nentries_30 = 0.30*nentries;
  nentries_50 = 0.50*nentries;
  
  for(Int_t k=0;k<nentries_20;k++)
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
      Double_t M_tau_Edge;
      Bool_t flagFOMS = calculateFOMS(pa,pb,Ecms,M_tau_Edge);
      Mtau = M_tau_Edge;
      //By default is flagMOMS is false, Ymin and Ymax are equal to -1000
      
      idp = ident;
      iddata = 1;

      
      treef->Fill();
    }        
  

  for(Int_t k=nentries_20;k<nentries_50;k++)
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
      Double_t M_tau_Edge;
      Bool_t flagFOMS = calculateFOMS(pa,pb,Ecms,M_tau_Edge);
      Mtau = M_tau_Edge;
      //By default is flagMOMS is false, Ymin and Ymax are equal to -1000
      
      idp = ident;
      iddata = 2 ;

      
      treef->Fill();
    }  

    for(Int_t k=nentries_50;k<nentries;k++)
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
      Double_t M_tau_Edge;
      Bool_t flagFOMS = calculateFOMS(pa,pb,Ecms,M_tau_Edge);
      Mtau = M_tau_Edge;
      //By default is flagMOMS is false, Ymin and Ymax are equal to -1000
      
      idp = ident;
      iddata = 3;

      
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
