#include "TFile.h"
#include "MyPlanePack.h"

using namespace std;


void MyPlanePack::ResetEvent() {
	Weight = 1.0;

	Full.ResetEvent();
	East.ResetEvent();
	West.ResetEvent();

	IsGood = false;
	for(int i=0; i<_NMAX+1; i++) {
		RawReso[i] = 0;
		CorReso[i] = 0;
	}
}


void MyPlanePack::ClaimHistogram() {
	hCosFull = new TProfile2D(Name+"CosFull", Name+"CosFull", _NMAX,0.5,_NMAX+0.5, _KMAX,0.5,_KMAX+0.5);
	hSinFull = new TProfile2D(Name+"SinFull", Name+"SinFull", _NMAX,0.5,_NMAX+0.5, _KMAX,0.5,_KMAX+0.5);
	hCosEast = new TProfile2D(Name+"CosEast", Name+"CosEast", _NMAX,0.5,_NMAX+0.5, _KMAX,0.5,_KMAX+0.5);
	hSinEast = new TProfile2D(Name+"SinEast", Name+"SinEast", _NMAX,0.5,_NMAX+0.5, _KMAX,0.5,_KMAX+0.5);
	hCosWest = new TProfile2D(Name+"CosWest", Name+"CosWest", _NMAX,0.5,_NMAX+0.5, _KMAX,0.5,_KMAX+0.5);
	hSinWest = new TProfile2D(Name+"SinWest", Name+"SinWest", _NMAX,0.5,_NMAX+0.5, _KMAX,0.5,_KMAX+0.5);
	hCosFull->Sumw2();
	hSinFull->Sumw2();
	hCosEast->Sumw2();
	hSinEast->Sumw2();
	hCosWest->Sumw2();
	hSinWest->Sumw2();

	for(int i=1; i<_NMAX+1; i++) {
		hRawFullDist[i] = new TH1D(Name+Form("RawFullDistPsi%d",i), Name+Form("RawFullDistPsi%d",i), 120,0.0,2*M_PI/i);
		hRawEastDist[i] = new TH1D(Name+Form("RawEastDistPsi%d",i), Name+Form("RawEastDistPsi%d",i), 120,0.0,2*M_PI/i);
		hRawWestDist[i] = new TH1D(Name+Form("RawWestDistPsi%d",i), Name+Form("RawWestDistPsi%d",i), 120,0.0,2*M_PI/i);
		hRawSubsDist[i] = new TH2D(Name+Form("RawSubsDistPsi%d",i), Name+Form("RawSubsDistPsi%d",i), 32,0.0,2*M_PI/i, 32,0.0,2*M_PI/i);
		hCorFullDist[i] = new TH1D(Name+Form("CorFullDistPsi%d",i), Name+Form("CorFullDistPsi%d",i), 120,0.0,2*M_PI/i);
		hCorEastDist[i] = new TH1D(Name+Form("CorEastDistPsi%d",i), Name+Form("CorEastDistPsi%d",i), 120,0.0,2*M_PI/i);
		hCorWestDist[i] = new TH1D(Name+Form("CorWestDistPsi%d",i), Name+Form("CorWestDistPsi%d",i), 120,0.0,2*M_PI/i);
		hCorSubsDist[i] = new TH2D(Name+Form("CorSubsDistPsi%d",i), Name+Form("CorSubsDistPsi%d",i), 32,0.0,2*M_PI/i, 32,0.0,2*M_PI/i);
		hRawFullDist[i]->Sumw2();
		hRawEastDist[i]->Sumw2();
		hRawWestDist[i]->Sumw2();
		hRawSubsDist[i]->Sumw2();
		hCorFullDist[i]->Sumw2();
		hCorEastDist[i]->Sumw2();
		hCorWestDist[i]->Sumw2();
		hCorSubsDist[i]->Sumw2();
	}
}


bool MyPlanePack::InputHistogram(TString inputpath) {
	TFile *fInput = new TFile(inputpath, "READ");
	if(fInput->IsZombie()) {
		cout << "MyPlanePack: Invalid fInput from " << inputpath << "!!!" << endl;
		cout << "Only RawPsi this time." << endl;
		return false;
	}

	TProfile2D *hInputCosFull = (TProfile2D*)fInput->Get(Name+"CosFull");
	TProfile2D *hInputSinFull = (TProfile2D*)fInput->Get(Name+"SinFull");
	TProfile2D *hInputCosEast = (TProfile2D*)fInput->Get(Name+"CosEast");
	TProfile2D *hInputSinEast = (TProfile2D*)fInput->Get(Name+"SinEast");
	TProfile2D *hInputCosWest = (TProfile2D*)fInput->Get(Name+"CosWest");
	TProfile2D *hInputSinWest = (TProfile2D*)fInput->Get(Name+"SinWest");

	IsGoodInput = true;
	if(!hInputCosFull) IsGoodInput = false;
	if(!hInputSinFull) IsGoodInput = false;
	if(!hInputCosEast) IsGoodInput = false;
	if(!hInputSinEast) IsGoodInput = false;
	if(!hInputCosWest) IsGoodInput = false;
	if(!hInputSinWest) IsGoodInput = false;
	if(!IsGoodInput) {
		cout << "MyPlanePack: Invalid hInput from " << inputpath << "!!!" << endl;
		cout << "Only RawPsi this time." << endl;
		return false;
	}
	hInputCosFull->SetName("Input"+Name+"CosFull");
	hInputSinFull->SetName("Input"+Name+"SinFull");
	hInputCosEast->SetName("Input"+Name+"CosEast");
	hInputSinEast->SetName("Input"+Name+"SinEast");
	hInputCosWest->SetName("Input"+Name+"CosWest");
	hInputSinWest->SetName("Input"+Name+"SinWest");

	Full.InputShift(hInputCosFull, hInputSinFull);
	West.InputShift(hInputCosWest, hInputSinWest);
	East.InputShift(hInputCosEast, hInputSinEast);

	fInput->Close("R");

	cout << "MyPlanePack: " << Name << " shift file loaded." << endl;

	return true;
}


bool MyPlanePack::InputHistogram(TDirectoryFile *fInput) {
	if(!fInput || fInput->IsZombie()) {
		cout << "MyPlanePack: Invalid fInput!!!" << endl;
		cout << "Only RawPsi this time." << endl;
		return false;
	}

	TProfile2D *hInputCosFull = (TProfile2D*)fInput->Get(Name+"CosFull");
	TProfile2D *hInputSinFull = (TProfile2D*)fInput->Get(Name+"SinFull");
	TProfile2D *hInputCosEast = (TProfile2D*)fInput->Get(Name+"CosEast");
	TProfile2D *hInputSinEast = (TProfile2D*)fInput->Get(Name+"SinEast");
	TProfile2D *hInputCosWest = (TProfile2D*)fInput->Get(Name+"CosWest");
	TProfile2D *hInputSinWest = (TProfile2D*)fInput->Get(Name+"SinWest");

	IsGoodInput = true;
	if(!hInputCosFull) IsGoodInput = false;
	if(!hInputSinFull) IsGoodInput = false;
	if(!hInputCosEast) IsGoodInput = false;
	if(!hInputSinEast) IsGoodInput = false;
	if(!hInputCosWest) IsGoodInput = false;
	if(!hInputSinWest) IsGoodInput = false;
	if(!IsGoodInput) {
		cout << "MyPlanePack: Invalid hInput from!!!" << endl;
		cout << "Only RawPsi this time." << endl;
		return false;
	}
	hInputCosFull->SetName("Input"+Name+"CosFull");
	hInputSinFull->SetName("Input"+Name+"SinFull");
	hInputCosEast->SetName("Input"+Name+"CosEast");
	hInputSinEast->SetName("Input"+Name+"SinEast");
	hInputCosWest->SetName("Input"+Name+"CosWest");
	hInputSinWest->SetName("Input"+Name+"SinWest");

	Full.InputShift(hInputCosFull, hInputSinFull);
	West.InputShift(hInputCosWest, hInputSinWest);
	East.InputShift(hInputCosEast, hInputSinEast);

	cout << "MyPlanePack: " << Name << " shift file loaded." << endl;

	return true;
}


MyPlanePack::MyPlanePack() {
	Name = "";
	IsGoodInput = false;
	IsAutoFill = true;

	Full = MyPlaneInfo(Name+"Full");
	East = MyPlaneInfo(Name+"East");
	West = MyPlaneInfo(Name+"West");

	ResetEvent();
	ClaimHistogram();
}


MyPlanePack::MyPlanePack(TString name) {
	Name = name;
	IsGoodInput = false;
	IsAutoFill = true;

	Full = MyPlaneInfo(Name+"Full");
	East = MyPlaneInfo(Name+"East");
	West = MyPlaneInfo(Name+"West");

	ResetEvent();
	ClaimHistogram();
}


MyPlanePack::~MyPlanePack() {
}


void MyPlanePack::Add(double phi, double eta, double wts) {
	int isub = eta<0?0:1;
	double mphi = phi + isub*M_PI; // eta sign effect

	Full.Add(mphi, wts);
	if(isub==0) East.Add(mphi, wts);
	if(isub==1) West.Add(mphi, wts);
}


void MyPlanePack::Fill() {
	if(!IsGood) return;

	if(Full.GetIsGood()) Full.OutputShift(hCosFull, hSinFull, Weight);
	if(West.GetIsGood()) West.OutputShift(hCosWest, hSinWest, Weight);
	if(East.GetIsGood()) East.OutputShift(hCosEast, hSinEast, Weight);

	for(int i=1; i<_NMAX+1; i++) {
		hRawFullDist[i]->Fill(Full.GetRawPsi(i), Weight);
		hRawEastDist[i]->Fill(East.GetRawPsi(i), Weight);
		hRawWestDist[i]->Fill(West.GetRawPsi(i), Weight);
		hRawSubsDist[i]->Fill(East.GetRawPsi(i), West.GetRawPsi(i), Weight);
		hCorFullDist[i]->Fill(Full.GetCorPsi(i), Weight);
		hCorEastDist[i]->Fill(East.GetCorPsi(i), Weight);
		hCorWestDist[i]->Fill(West.GetCorPsi(i), Weight);
		hCorSubsDist[i]->Fill(East.GetCorPsi(i), West.GetCorPsi(i), Weight);
	}
}


void MyPlanePack::Calc() {
	Full.Calc();
	East.Calc();
	West.Calc();

	IsGood = East.GetIsGood() && West.GetIsGood();
	for(int i=1; i<_NMAX+1; i++) {
		RawReso[i] = cos(i*(East.GetRawPsi(i)-West.GetRawPsi(i)));
		CorReso[i] = cos(i*(East.GetCorPsi(i)-West.GetCorPsi(i)));
	}

	if(IsAutoFill) Fill();
}


void MyPlanePack::Finish() {
	hCosFull->Write();
	hSinFull->Write();
	hCosEast->Write();
	hSinEast->Write();
	hCosWest->Write();
	hSinWest->Write();

	for(int i=1; i<_NMAX+1; i++) {
		hRawFullDist[i]->Write();
		hRawEastDist[i]->Write();
		hRawWestDist[i]->Write();
		hRawSubsDist[i]->Write();
		hCorFullDist[i]->Write();
		hCorEastDist[i]->Write();
		hCorWestDist[i]->Write();
		hCorSubsDist[i]->Write();
	}
}
