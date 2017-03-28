// Example: root -l 'RFPlotALaArnaud.C("~/godaq_rootfiles/analysis_v2.12-calibG2/run111.root")'

#include "TFile.h"
#include "TTree.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

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
	TH1F* hPeriod = new TH1F("hPeriod", "hPeriod", 100, 0, 100);
	
	TCanvas* cc = new TCanvas("cc", "cc");
        while (reader.Next() && *Evt < nEvents) {
		if(*Evt%500 == 0) {
			cout << "Event " << *Evt << endl;
		}
// 		if(*Evt != 233)
// 			continue;
		if(*NoLORs >= 10) {
			continue;
		}
		TGraph* gRF = new TGraph(999);
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
			
			if(LORTMean[j] < 20 || LORTMean[j] > 50) {
				continue;
			}
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
		delete gRF;
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
}
