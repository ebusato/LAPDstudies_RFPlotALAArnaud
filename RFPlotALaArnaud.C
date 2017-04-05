// Example: root -l 'RFPlotALaArnaud.C("~/godaq_rootfiles/analysis_v2.12-calibG2/run111.root")'
// h2 can be partly obtained quickly by doing:
//    .L RFPlotALaArnaud.C 
//    RFPlotALaArnaudDirectly("T30[LORIdx1] > 20 && T30[LORIdx1] < 50 && T30[LORIdx2] > 20 && T30[LORIdx2] < 50", "analysis_v2.16-calibG2/run110LOR.root", "analysis_v2.16-calibG2/run111LOR.root")
//    RFPlotALaArnaudDirectly("NoLORs == 1", "analysis_v2.14-calibG2/run110.root", "analysis_v2.14-calibG2/run111.root")

#include "TFile.h"
#include "TTree.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using namespace RooFit;

TF1* fitSine(TGraph* g) {
	TF1* f = new TF1("f","[0] + [1]*TMath::Sin([2]*x + [3])", 0, 200);
	f->SetParameter(2, 2*TMath::Pi()*24.85e6/1e9);
	f->SetParLimits(2, 2*TMath::Pi()/20., 2*TMath::Pi()/60.);
	f->SetLineColor(kRed);
	g->Fit("f", "Q", "", 0, 200);
	return f;
}

void RFPlotALaArnaud(int nEvents, TString fileName0, TString fileName1="", TString fileName2="", TString fileName3="", TString fileName4="", 
		       TString fileName5="", TString fileName6="", TString fileName7="", 
		       TString fileName8="", TString fileName9="") {
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(1);
	
	TChain ch("tree");
	ch.Add(fileName0);
	if(fileName1 != "") {
		ch.Add(fileName1);
	}
	if(fileName2 != "") {
		ch.Add(fileName2);
	}
	if(fileName3 != "") {
		ch.Add(fileName3);
	}
	if(fileName4 != "") {
		ch.Add(fileName4);
	}
	if(fileName5 != "") {
		ch.Add(fileName5);
	}
	if(fileName6 != "") {
		ch.Add(fileName6);
	}
	if(fileName7 != "") {
		ch.Add(fileName7);
	}
	if(fileName8 != "") {
		ch.Add(fileName8);
	}
	if(fileName9 != "") {
		ch.Add(fileName9);
	}
	
        TTreeReader reader(&ch);
        TTreeReaderValue<UInt_t> Run(reader, "Run");
        TTreeReaderValue<UInt_t> Evt(reader, "Evt");
	TTreeReaderArray<Double_t> E(reader, "E");
	TTreeReaderArray<Double_t> PulseRF(reader, "PulseRF");
	TTreeReaderArray<Double_t> SampleTimes(reader, "SampleTimes");
	TTreeReaderValue<Int_t> NoLORs(reader, "NoLORs");
        TTreeReaderArray<Int_t> LORIdx1(reader, "LORIdx1");
        TTreeReaderArray<Int_t> LORIdx2(reader, "LORIdx2");
	TTreeReaderArray<Double_t> LORTMean(reader, "LORTMean");
	TTreeReaderArray<Double_t> LORTRF(reader, "LORTRF");

	TH2F* h1 = new TH2F("h1", "h1", 100, 0, 40, 100, 0, 600);
	TH2F* h2 = new TH2F("h2", "h2", 100, 0, 40, 100, 0, 600);
	
	TH1F* hDeltaTRF = new TH1F("hDeltaTRF", "hDeltaTRF", 100, -100, 100);
	TH1F* hPeriod = new TH1F("hPeriod", "hPeriod", 5000, 0, 100);
	
	//TCanvas* cc = new TCanvas("cc", "cc");
	TGraph* gRF = new TGraph(999);
        while (reader.Next() && *Evt < nEvents) {
		if(*Evt%500 == 0) {
			cout << "Event " << *Evt << endl;
		}
// 		if(*Evt != 233)
// 			continue;
		if(*NoLORs >= 10) {
			continue;
		}
		
		for(int i = 0; i < 999; i++) {
			gRF->SetPoint(i, SampleTimes[i], PulseRF[i]);
		}
// 		gRF->Draw("ap");
// 		gPad->SetGridx();
// 		gPad->SetGridy();
		TF1* f = fitSine(gRF);
		double omega = f->GetParameter(2);
		double phase = f->GetParameter(3);
		double period = 2*TMath::Pi()/omega;
		hPeriod->Fill(period);
		std::vector<double> roots;
		int p = 0;
		while(roots.size() < 5) {
			double root = (TMath::Pi()*p-phase)/omega;
			if(root > 0 && f->Derivative(root) > 0) { 
				roots.push_back(root);
			}
			p++;
		}
		
		//cout << "roots size: " << roots.size() << endl;
		
// 		cout << (-1*phase)/omega << " " << f->Derivative((-1*phase)/omega) << " " 
// 			(2*TMath::Pi()*1-phase)/omega << " " << f->Derivative((2*TMath::Pi()*1-1*phase)/omega) << " "
// 			(2*TMath::Pi()*2-phase)/omega << " " << endl;
		
// 		for(int i=0; i<5; i++) {
// 			cout << "root " << i << " = " << roots[i] << " " << f->Derivative(roots[i]) << endl;
// 		}
		
		//cout << "NoLORs = " << *NoLORs << endl;
		for(int j = 0; j < *NoLORs; j++) {
// 			cout << "IDs: " << LORIdx1[j] << " "<< LORIdx2[j] << endl;
// 			cout << " Energies: " << E[LORIdx1[j]] << " " << E[LORIdx2[j]] << endl;
// 			cout << " Times: " << LORTRF[j] << " " << LORTMean[j] << endl;
			h1 -> Fill(LORTMean[j] - LORTRF[j], E[LORIdx1[j]]);
			h1 -> Fill(LORTMean[j] - LORTRF[j], E[LORIdx2[j]]);
			
			double TRFfit = 0;
			if(LORTMean[j] <= roots[0]) {
				TRFfit = roots[0] - period;
			} else if(LORTMean[j] >= roots[roots.size()-1]) {
				TRFfit = roots[roots.size()-1];
			} else {
				for(int k = 0; k < roots.size(); k++) {
					if(k < roots.size() -1) {
						if(LORTMean[j] > roots[k] && LORTMean[j] < roots[k+1]) {
							TRFfit = roots[k];
							break;
						}
					} else {
						cout << "LORTMean[j] = " << LORTMean[j] << endl;
						cout << "roots[k] = " << roots[k] << endl;
						cout << "ERROR, this should never happen" << endl;
						exit(0);
					}
					
				}
			}
			if(LORTMean[j] - TRFfit > period) {
				cout << "LORTMean[j] = " << LORTMean[j] << endl;
				cout << "TRFfit = " << TRFfit << endl;
				cout << "period = " << period << endl;
				cout << "ERROR, this should not happen (Evt = " <<  *Evt << ")" << endl;
				//exit(0);
			}
			
			h2 -> Fill(LORTMean[j] - TRFfit, E[LORIdx1[j]]);
			h2 -> Fill(LORTMean[j] - TRFfit, E[LORIdx2[j]]);
			
			hDeltaTRF->Fill(TRFfit - LORTRF[j]);
		}
		delete f;
        }
	
	TCanvas* c1 = new TCanvas("c1", "c1");
	h1->Draw("colz");
	TCanvas* c2 = new TCanvas("c2", "c2");
	h2->Draw("colz");
	TCanvas* c3 = new TCanvas("c3", "c3");
	hDeltaTRF->Draw();
	TCanvas* c4 = new TCanvas("c4", "c4");
	hPeriod->Draw();
	
	TFile* fOut = new TFile("RFPlotALaArnaud.root", "recreate");
	h1->Write();
	h2->Write();
	hDeltaTRF->Write();
	hPeriod->Write();
	fOut->Write();
}

void RFPlotALaArnaudDirectly(TCut cut, TString fileName0, TString fileName1="", TString fileName2="", TString fileName3="", TString fileName4="", 
		       TString fileName5="", TString fileName6="", TString fileName7="", 
		       TString fileName8="", TString fileName9="") {
	TChain ch("tree");
	ch.Add(fileName0.Data());
	if(fileName1 != "") {
		ch.Add(fileName1);
	}
	if(fileName2 != "") {
		ch.Add(fileName2);
	}
	if(fileName3 != "") {
		ch.Add(fileName3);
	}
	if(fileName4 != "") {
		ch.Add(fileName4);
	}
	if(fileName5 != "") {
		ch.Add(fileName5);
	}
	if(fileName6 != "") {
		ch.Add(fileName6);
	}
	if(fileName7 != "") {
		ch.Add(fileName7);
	}
	if(fileName8 != "") {
		ch.Add(fileName8);
	}
	if(fileName9 != "") {
		ch.Add(fileName9);
	}
	TCanvas* c0 = new TCanvas("c0", "c0");
	ch.Draw("RateLvsL3 : Entry$");
	
	TCanvas* c1 = new TCanvas("c1", "c1");
	TH2F* hArnaud = new TH2F("hArnaud", "hArnaud", 100, 0, 40, 200, 0, 1000);
	ch.Draw("E[LORIdx1] : LORTMean - LORTRF>>hArnaud",  cut, "colz");
	ch.Draw("E[LORIdx2] : LORTMean - LORTRF>>+hArnaud",  cut, "colz");
	hArnaud->Draw("colz");
	
	TCanvas* c2 = new TCanvas("c2", "c2");
	TH1F* hESpillOut = new TH1F("hESpillOut", "hESpillOut", 100, 0, 1100);
	TH1F* hESpillIn = new TH1F("hESpillIn", "hESpillIn", 100, 0, 1100);
	hESpillIn->SetLineColor(kRed);
	ch.Draw("E[LORIdx1]>>hESpillOut", cut && "abs(LORTMean - LORTRF - 7) > 5");
	ch.Draw("E[LORIdx2]>>+hESpillOut", cut && "abs(LORTMean - LORTRF - 7) > 5");
	ch.Draw("E[LORIdx1]>>hESpillIn", cut && "abs(LORTMean - LORTRF - 7) < 5");
	ch.Draw("E[LORIdx2]>>+hESpillIn", cut && "abs(LORTMean - LORTRF - 7) < 5");
	hESpillOut->Draw();
	hESpillIn->Draw("same");
	
	TCanvas* c3 = new TCanvas("c3", "c3");
	TH1F* hZmarSpillOut = new TH1F("hZmarSpillOut", "hZmarSpillOut", 100, -100, 100);
	TH1F* hZmarSpillIn = new TH1F("hZmarSpillIn", "hZmarSpillIn", 100, -100, 100);
	hZmarSpillIn->SetLineColor(kRed);
	ch.Draw("LORZmar>>hZmarSpillOut", cut && "abs(LORTMean - LORTRF - 7) > 5");
	ch.Draw("LORZmar>>+hZmarSpillOut", cut && "abs(LORTMean - LORTRF - 7) > 5");
	ch.Draw("LORZmar>>hZmarSpillIn", cut && "abs(LORTMean - LORTRF - 7) < 5");
	ch.Draw("LORZmar>>+hZmarSpillIn", cut && "abs(LORTMean - LORTRF - 7) < 5");
	hZmarSpillIn->Scale(1/hZmarSpillIn->Integral());
	hZmarSpillOut->Scale(1/hZmarSpillOut->Integral());
	hZmarSpillOut->Draw();
	hZmarSpillIn->Draw("same");
	hZmarSpillOut->Draw("same");
}


TH1F* Draw(TTree* t, TString var, TCut cut, TString hName, int Nbins, double xmin, double xmax, int color, int linesize, bool normalize=false)
{
	TH1F* h = new TH1F(hName.Data(), hName.Data(), Nbins, xmin, xmax);
	TString temp(var);
	temp+=">>";
	temp+=hName;
	t->Draw(temp.Data(), cut);
	h->SetLineColor(color);
	h->SetLineWidth(linesize);
	if(normalize) {
		h->Scale(1/h->Integral());
	}
	return h;
}

// RooDataSet* MakeDataSetFromTH1(TH1F* h)
TH1F* MakeKernelPDFFromTH1(TH1F* h)
{
  RooRealVar* z = new RooRealVar("z", "z", -100, 100);
  z->setBins(1000);
  RooDataSet* ds = new RooDataSet("ds","ds",RooArgSet(*z)) ;
  for(int i=0; i<h->GetNbinsX(); i++) {
	double binContent = h->GetBinContent(i);
	if(binContent!=0) {
	  double binCenter = h->GetBinCenter(i);
	  z->setVal(binCenter);
	  for(int j=0; j<binContent; j++) {
	  	ds->add(RooArgSet(*z));
	  }
	}
  }
  RooDataHist* dh = new RooDataHist("dh", "dh", *z, Import(*h));
  
  RooKeysPdf kest1("kest1","kest1",*z,*ds,RooKeysPdf::MirrorBoth, 1) ;
//   RooKeysPdf kest1("kest1","kest1",*z,*ds,RooKeysPdf::NoMirror) ;
//   RooPlot* frame = z->frame() ;
// //   ds->plotOn(frame);
//   dh->plotOn(frame);
//   kest1.plotOn(frame);
//   frame->Draw();
//   
  TH1F* hKeys = (TH1F*) kest1.createHistogram("hKeys", *z);
//   hKeys->Draw("histsame");
  return hKeys;
}
  //RooDataSet* data = new RooDataSet("data","data",RooArgSet(x));
  //x=4;
  //data->add(RooArgSet(x));


void MakeSpillOutPlots()
{
	TFile* f0 = new TFile("analysis_v2.18-calibG2/run91LOR.root", "read");
	TFile* f1 = new TFile("analysis_v2.18-calibG2/run110LOR.root", "read");
// 	TFile* f0 = new TFile("analysis_v2.18-calibG2/run98LOR.root", "read");
// 	TFile* f1 = new TFile("analysis_v2.18-calibG2/run99LOR.root", "read");
	
	TTree* t0 = (TTree*) f0->Get("tree");
	TTree* t1 = (TTree*) f1->Get("tree");
	
	TCut cutTimes("T30[LORIdx1] > 20 && T30[LORIdx1] < 50 && T30[LORIdx2] > 20 && T30[LORIdx2] < 50");
	TCut cutSpillOut("abs(LORTMean - LORTRF - 7) > 5");
	TCut cut = cutTimes && cutSpillOut;
	
// 	TCut cut0 = "Evt > 2000 && Evt < 60000";
// 	TCut cut1 = "Evt > 2000 && Evt < 60000";
	
	TCut cut0 = "Evt > 2000 && Evt < 60000";
	TCut cut1 = "Evt > 2000 && Evt < 60000";
	
	TCanvas* c0 = new TCanvas("c0", "c0");
	c0->Divide(2,2);
	c0->cd(1);
	t0->Draw("RateLvsL3 : Evt");
	c0->cd(2);
	t1->Draw("RateLvsL3 : Evt");
	c0->cd(3);
	t0->Draw("RateLvsL3 : Evt", cut0);
	c0->cd(4);
	t1->Draw("RateLvsL3 : Evt", cut1);
	
	
	TCanvas* c1 = new TCanvas("c1", "c1");
	c1->Divide(2,1);
	c1->cd(1);
	TH2F* hArnaud_0 = new TH2F("hArnaud_0", "hArnaud_0", 100, 0, 40, 200, 0, 1000);
	t0->Draw("E[LORIdx1] : LORTMean - LORTRF>>hArnaud_0",  cut0 && cutTimes, "colz");
	t0->Draw("E[LORIdx2] : LORTMean - LORTRF>>+hArnaud_0",  cut0 && cutTimes, "colz");
	hArnaud_0->Draw("colz");
	c1->cd(2);
	TH2F* hArnaud_1 = new TH2F("hArnaud_1", "hArnaud_1", 100, 0, 40, 200, 0, 1000);
	t1->Draw("E[LORIdx1] : LORTMean - LORTRF>>hArnaud_1",  cut1 && cutTimes, "colz");
	t1->Draw("E[LORIdx2] : LORTMean - LORTRF>>+hArnaud_1",  cut1 && cutTimes, "colz");
	hArnaud_1->Draw("colz");
	
	TCanvas* c2 = new TCanvas("c2", "c2");
	TH1F* hETemp_0 = Draw(t0, "E[LORIdx1]", cut0 && cut, "hE_01", 100, 0, 1000, kRed, 1);
	TH1F* hE_0 = Draw(t0, "E[LORIdx2]", cut0 && cut, "hE_02", 100, 0, 1000, kRed, 1);
	hE_0->Add(hETemp_0);
	TH1F* hETemp_1 = Draw(t1, "E[LORIdx1]", cut1 && cut, "hE_11", 100, 0, 1000, kGreen+2, 1);
	TH1F* hE_1 = Draw(t1, "E[LORIdx2]", cut1 && cut, "hE_12", 100, 0, 1000, kGreen+2, 1);
	hE_1->Add(hETemp_1);
	hE_0->Draw();
	hE_1->Draw("same");
	
	TCanvas* c3 = new TCanvas("c3", "c3");
	TH1F* hZmar_0 = Draw(t0, "LORZmar", cut0 && cut, "hZmar", 2000, -100, 100, kRed, 2);
	TH1F* hZmar_1 = Draw(t1, "LORZmar", cut1 && cut, "hZmar", 2000, -100, 100, kGreen+2, 1);
	TH1F* hKeys_0 = MakeKernelPDFFromTH1(hZmar_0);
	TH1F* hKeys_1 = MakeKernelPDFFromTH1(hZmar_1);
        hKeys_0->SetLineColor(kMagenta);
	hKeys_1->SetLineColor(kBlue);
	hZmar_0->Scale(1/hZmar_0->Integral());
	hZmar_0->Draw();
	hZmar_1->Scale(1/hZmar_1->Integral());
	hZmar_1->Draw("same");
	
// 	TCanvas* c4 = new TCanvas("c4", "c4");
// 	RooDataSet* data_0 = MakeDataSetFromTH1(hZmar_0);
	hKeys_0->Scale(16/hKeys_0->Integral());
	hKeys_0->Draw("same");
	hKeys_1->Scale(16/hKeys_1->Integral());
	hKeys_1->Draw("same");
// 	h->Draw();
	
}