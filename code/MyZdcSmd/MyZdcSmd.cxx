#include "MyZdcSmd.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;


void MyZdcSmd::ResetEvent() {
	mCent = 0;
	mWght = 1.0;
	mIsGoodEvent = true;

	mPv.SetXYZ(0.0, 0.0, 0.0);
	mPvx = mPv.X();
	mPvy = mPv.Y();

	// [East,West][Vert,Hori][Slats]
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				mZdcSmdRawAdc[iEW][iVH][iST] = 0;
				mZdcSmdCorAdc[iEW][iVH][iST] = 0;
			}
		}
	}

	// Raw: Raw
	// Cor: Raw + PedeGain
	// Rct: Raw + PedeGain + Recenter
	// Shf: Raw + PedeGain + Recenter + Shift
	mEastRawQ.Set(0.0, 0.0);
	mEastCorQ.Set(0.0, 0.0);
	mEastRctQ.Set(0.0, 0.0);
	mEastShfQ.Set(0.0, 0.0);
	mWestRawQ.Set(0.0, 0.0);
	mWestCorQ.Set(0.0, 0.0);
	mWestRctQ.Set(0.0, 0.0);
	mWestShfQ.Set(0.0, 0.0);
	mFullRawQ.Set(0.0, 0.0);
	mFullCorQ.Set(0.0, 0.0);
	mFullRctQ.Set(0.0, 0.0);
	mFullShfQ.Set(0.0, 0.0);
	mEastRawPsi = -999.;
	mEastCorPsi = -999.;
	mEastRctPsi = -999.;
	mEastShfPsi = -999.;
	mWestRawPsi = -999.;
	mWestCorPsi = -999.;
	mWestRctPsi = -999.;
	mWestShfPsi = -999.;
	mFullRawPsi = -999.;
	mFullCorPsi = -999.;
	mFullRctPsi = -999.;
	mFullShfPsi = -999.;

	// PvRct: Raw + PedeGain + PvRecenter
	// PvShf: Raw + PedeGain + PvRecenter + Shift
	mEastPvRctQ.Set(0.0, 0.0);
	mEastPvShfQ.Set(0.0, 0.0);
	mWestPvRctQ.Set(0.0, 0.0);
	mWestPvShfQ.Set(0.0, 0.0);
	mFullPvRctQ.Set(0.0, 0.0);
	mFullPvShfQ.Set(0.0, 0.0);
	mEastPvRctPsi = -999.;
	mEastPvShfPsi = -999.;
	mWestPvRctPsi = -999.;
	mWestPvShfPsi = -999.;
	mFullPvRctPsi = -999.;
	mFullPvShfPsi = -999.;
}


void MyZdcSmd::Init() {
	ResetEvent();

	mVerbose = 2;
	mMerge = 0;
	mGainMod = 0;
	mCenterMod = 0;
	mIsGoodInput = true;
	mIsAutoFill = true; // default true; if false, must call Fill() later

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				Position[iEW][iVH][iST] = GetPosition(iEW, iVH, iST);
				Poserror[iEW][iVH][iST] = 0;
			}
		}
	}

	// record the parameters from the first calculation, keep consistency
	hZdcSmdPede = new TProfile(mName+"ZdcSmdPede", mName+"ZdcSmdPede", nAL,0,nAL);
	hZdcSmdGain = new TProfile(mName+"ZdcSmdGain", mName+"ZdcSmdGain", nAL,0,nAL);
	hZdcSmdPeak = new TProfile(mName+"ZdcSmdPeak", mName+"ZdcSmdPeak", nAL,0,nAL);
	hZdcRawCenter = new TProfile(mName+"ZdcRawCenter", mName+"ZdcRawCenter", 4,0,4);
	hZdcCorCenter = new TProfile(mName+"ZdcCorCenter", mName+"ZdcCorCenter", 4,0,4);
	hZdcSmdPede->Sumw2();
	hZdcSmdGain->Sumw2();
	hZdcSmdPeak->Sumw2();
	hZdcRawCenter->Sumw2();
	hZdcCorCenter->Sumw2();

	// [East,West][Vert,Hori][Slats]
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				mZdcSmdPede[iEW][iVH][iST] = 0;
				mZdcSmdGain[iEW][iVH][iST] = 1.0;
				mZdcSmdPeak[iEW][iVH][iST] = 0;
			}
		}
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				TString StrTmp = Form("EW%dVH%dST%d", iEW, iVH, iST);
				hZdcSmdRawAdc[iEW][iVH][iST] = new TH1D(mName+"ZdcSmdRawAdc"+StrTmp, mName+"ZdcSmdRawAdc"+StrTmp, 3000,-5,2995);
				//hZdcSmdCorAdc[iEW][iVH][iST] = new TH1D(mName+"ZdcSmdCorAdc"+StrTmp, mName+"ZdcSmdCorAdc"+StrTmp, 3000,-5,2995);
				hZdcSmdRawAdc[iEW][iVH][iST]->Sumw2();
				//hZdcSmdCorAdc[iEW][iVH][iST]->Sumw2();
			}
		}
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			TString StrTmp = Form("EW%dVH%d", iEW, iVH);
			hRawCenter[iEW][iVH] = new TProfile(mName+"RawCenter"+StrTmp, mName+"RawCenter"+StrTmp, 1,0,1);
			hCorCenter[iEW][iVH] = new TProfile(mName+"CorCenter"+StrTmp, mName+"CorCenter"+StrTmp, 1,0,1);
			hRawCenter[iEW][iVH]->Sumw2();
			hCorCenter[iEW][iVH]->Sumw2();
			RawCenter[iEW][iVH] = 0;
			CorCenter[iEW][iVH] = 0;
			RefCenter[iEW][iVH] = 0;
		}
	}

	hAvePvx = new TProfile(mName+"AvePvx", mName+"AvePvx", 1,0,1);
	hAvePvy = new TProfile(mName+"AvePvy", mName+"AvePvy", 1,0,1);
	hAvePvx->Sumw2();
	hAvePvy->Sumw2();
	AvePvx = 0;
	AvePvy = 0;

	for(int iCent=0; iCent<nCent; iCent++) {
		TString StrCent = Form("Cent%d", iCent);
		hEastRctPsiCos[iCent] = new TProfile(mName+"EastRctPsiCos"+StrCent, mName+"EastRctPsiCos"+StrCent, nShf,0.5,nShf+0.5);
		hEastRctPsiSin[iCent] = new TProfile(mName+"EastRctPsiSin"+StrCent, mName+"EastRctPsiSin"+StrCent, nShf,0.5,nShf+0.5);
		hWestRctPsiCos[iCent] = new TProfile(mName+"WestRctPsiCos"+StrCent, mName+"WestRctPsiCos"+StrCent, nShf,0.5,nShf+0.5);
		hWestRctPsiSin[iCent] = new TProfile(mName+"WestRctPsiSin"+StrCent, mName+"WestRctPsiSin"+StrCent, nShf,0.5,nShf+0.5);
		hFullRctPsiCos[iCent] = new TProfile(mName+"FullRctPsiCos"+StrCent, mName+"FullRctPsiCos"+StrCent, nShf,0.5,nShf+0.5);
		hFullRctPsiSin[iCent] = new TProfile(mName+"FullRctPsiSin"+StrCent, mName+"FullRctPsiSin"+StrCent, nShf,0.5,nShf+0.5);
		hEastRctPsiCos[iCent]->Sumw2();
		hEastRctPsiSin[iCent]->Sumw2();
		hWestRctPsiCos[iCent]->Sumw2();
		hWestRctPsiSin[iCent]->Sumw2();
		hFullRctPsiCos[iCent]->Sumw2();
		hFullRctPsiSin[iCent]->Sumw2();
		for(int iShf=0; iShf<nShf; iShf++) {
			EastRctPsiCos[iCent][iShf] = 0;
			EastRctPsiSin[iCent][iShf] = 0;
			WestRctPsiCos[iCent][iShf] = 0;
			WestRctPsiSin[iCent][iShf] = 0;
			FullRctPsiCos[iCent][iShf] = 0;
			FullRctPsiSin[iCent][iShf] = 0;
		}
		hEastPvRctPsiCos[iCent] = new TProfile(mName+"EastPvRctPsiCos"+StrCent, mName+"EastPvRctPsiCos"+StrCent, nShf,0.5,nShf+0.5);
		hEastPvRctPsiSin[iCent] = new TProfile(mName+"EastPvRctPsiSin"+StrCent, mName+"EastPvRctPsiSin"+StrCent, nShf,0.5,nShf+0.5);
		hWestPvRctPsiCos[iCent] = new TProfile(mName+"WestPvRctPsiCos"+StrCent, mName+"WestPvRctPsiCos"+StrCent, nShf,0.5,nShf+0.5);
		hWestPvRctPsiSin[iCent] = new TProfile(mName+"WestPvRctPsiSin"+StrCent, mName+"WestPvRctPsiSin"+StrCent, nShf,0.5,nShf+0.5);
		hFullPvRctPsiCos[iCent] = new TProfile(mName+"FullPvRctPsiCos"+StrCent, mName+"FullPvRctPsiCos"+StrCent, nShf,0.5,nShf+0.5);
		hFullPvRctPsiSin[iCent] = new TProfile(mName+"FullPvRctPsiSin"+StrCent, mName+"FullPvRctPsiSin"+StrCent, nShf,0.5,nShf+0.5);
		hEastPvRctPsiCos[iCent]->Sumw2();
		hEastPvRctPsiSin[iCent]->Sumw2();
		hWestPvRctPsiCos[iCent]->Sumw2();
		hWestPvRctPsiSin[iCent]->Sumw2();
		hFullPvRctPsiCos[iCent]->Sumw2();
		hFullPvRctPsiSin[iCent]->Sumw2();
		for(int iShf=0; iShf<nShf; iShf++) {
			EastPvRctPsiCos[iCent][iShf] = 0;
			EastPvRctPsiSin[iCent][iShf] = 0;
			WestPvRctPsiCos[iCent][iShf] = 0;
			WestPvRctPsiSin[iCent][iShf] = 0;
			FullPvRctPsiCos[iCent][iShf] = 0;
			FullPvRctPsiSin[iCent][iShf] = 0;
		}
	}

	for(int iCent=0; iCent<nCent; iCent++) {
		TString StrCent = Form("Cent%d", iCent);
		hEastRawPsi[iCent] = new TH1D(mName+"EastRawPsi"+StrCent, mName+"EastRawPsi"+StrCent, 360,0,2*M_PI);
		hEastCorPsi[iCent] = new TH1D(mName+"EastCorPsi"+StrCent, mName+"EastCorPsi"+StrCent, 360,0,2*M_PI);
		hEastRctPsi[iCent] = new TH1D(mName+"EastRctPsi"+StrCent, mName+"EastRctPsi"+StrCent, 360,0,2*M_PI);
		hEastShfPsi[iCent] = new TH1D(mName+"EastShfPsi"+StrCent, mName+"EastShfPsi"+StrCent, 360,0,2*M_PI);
		hWestRawPsi[iCent] = new TH1D(mName+"WestRawPsi"+StrCent, mName+"WestRawPsi"+StrCent, 360,0,2*M_PI);
		hWestCorPsi[iCent] = new TH1D(mName+"WestCorPsi"+StrCent, mName+"WestCorPsi"+StrCent, 360,0,2*M_PI);
		hWestRctPsi[iCent] = new TH1D(mName+"WestRctPsi"+StrCent, mName+"WestRctPsi"+StrCent, 360,0,2*M_PI);
		hWestShfPsi[iCent] = new TH1D(mName+"WestShfPsi"+StrCent, mName+"WestShfPsi"+StrCent, 360,0,2*M_PI);
		hFullRawPsi[iCent] = new TH1D(mName+"FullRawPsi"+StrCent, mName+"FullRawPsi"+StrCent, 360,0,2*M_PI);
		hFullCorPsi[iCent] = new TH1D(mName+"FullCorPsi"+StrCent, mName+"FullCorPsi"+StrCent, 360,0,2*M_PI);
		hFullRctPsi[iCent] = new TH1D(mName+"FullRctPsi"+StrCent, mName+"FullRctPsi"+StrCent, 360,0,2*M_PI);
		hFullShfPsi[iCent] = new TH1D(mName+"FullShfPsi"+StrCent, mName+"FullShfPsi"+StrCent, 360,0,2*M_PI);
		h2SubRawPsi[iCent] = new TH2D(mName+"SubRawPsi"+StrCent, mName+"SubRawPsi"+StrCent, 32,0,2*M_PI, 32,0,2*M_PI);
		h2SubCorPsi[iCent] = new TH2D(mName+"SubCorPsi"+StrCent, mName+"SubCorPsi"+StrCent, 32,0,2*M_PI, 32,0,2*M_PI);
		h2SubRctPsi[iCent] = new TH2D(mName+"SubRctPsi"+StrCent, mName+"SubRctPsi"+StrCent, 32,0,2*M_PI, 32,0,2*M_PI);
		h2SubShfPsi[iCent] = new TH2D(mName+"SubShfPsi"+StrCent, mName+"SubShfPsi"+StrCent, 32,0,2*M_PI, 32,0,2*M_PI);
		hEastRawPsi[iCent]->Sumw2();
		hEastCorPsi[iCent]->Sumw2();
		hEastRctPsi[iCent]->Sumw2();
		hEastShfPsi[iCent]->Sumw2();
		hWestRawPsi[iCent]->Sumw2();
		hWestCorPsi[iCent]->Sumw2();
		hWestRctPsi[iCent]->Sumw2();
		hWestShfPsi[iCent]->Sumw2();
		hFullRawPsi[iCent]->Sumw2();
		hFullCorPsi[iCent]->Sumw2();
		hFullRctPsi[iCent]->Sumw2();
		hFullShfPsi[iCent]->Sumw2();
		h2SubRawPsi[iCent]->Sumw2();
		h2SubCorPsi[iCent]->Sumw2();
		h2SubRctPsi[iCent]->Sumw2();
		h2SubShfPsi[iCent]->Sumw2();
		hEastPvRctPsi[iCent] = new TH1D(mName+"EastPvRctPsi"+StrCent, mName+"EastPvRctPsi"+StrCent, 360,0,2*M_PI);
		hEastPvShfPsi[iCent] = new TH1D(mName+"EastPvShfPsi"+StrCent, mName+"EastPvShfPsi"+StrCent, 360,0,2*M_PI);
		hWestPvRctPsi[iCent] = new TH1D(mName+"WestPvRctPsi"+StrCent, mName+"WestPvRctPsi"+StrCent, 360,0,2*M_PI);
		hWestPvShfPsi[iCent] = new TH1D(mName+"WestPvShfPsi"+StrCent, mName+"WestPvShfPsi"+StrCent, 360,0,2*M_PI);
		hFullPvRctPsi[iCent] = new TH1D(mName+"FullPvRctPsi"+StrCent, mName+"FullPvRctPsi"+StrCent, 360,0,2*M_PI);
		hFullPvShfPsi[iCent] = new TH1D(mName+"FullPvShfPsi"+StrCent, mName+"FullPvShfPsi"+StrCent, 360,0,2*M_PI);
		h2SubPvRctPsi[iCent] = new TH2D(mName+"SubPvRctPsi"+StrCent, mName+"SubPvRctPsi"+StrCent, 32,0,2*M_PI, 32,0,2*M_PI);
		h2SubPvShfPsi[iCent] = new TH2D(mName+"SubPvShfPsi"+StrCent, mName+"SubPvShfPsi"+StrCent, 32,0,2*M_PI, 32,0,2*M_PI);
		hEastPvRctPsi[iCent]->Sumw2();
		hEastPvShfPsi[iCent]->Sumw2();
		hWestPvRctPsi[iCent]->Sumw2();
		hWestPvShfPsi[iCent]->Sumw2();
		hFullPvRctPsi[iCent]->Sumw2();
		hFullPvShfPsi[iCent]->Sumw2();
		h2SubPvRctPsi[iCent]->Sumw2();
		h2SubPvShfPsi[iCent]->Sumw2();
	}

	hResoRawPsi = new TProfile(mName+"ResoRawPsi", mName+"ResoRawPsi", nCent,0,nCent);
	hResoCorPsi = new TProfile(mName+"ResoCorPsi", mName+"ResoCorPsi", nCent,0,nCent);
	hResoRctPsi = new TProfile(mName+"ResoRctPsi", mName+"ResoRctPsi", nCent,0,nCent);
	hResoShfPsi = new TProfile(mName+"ResoShfPsi", mName+"ResoShfPsi", nCent,0,nCent);
	hResoRawPsi->Sumw2();
	hResoCorPsi->Sumw2();
	hResoRctPsi->Sumw2();
	hResoShfPsi->Sumw2();
	hResoPvRctPsi = new TProfile(mName+"ResoPvRctPsi", mName+"ResoPvRctPsi", nCent,0,nCent);
	hResoPvShfPsi = new TProfile(mName+"ResoPvShfPsi", mName+"ResoPvShfPsi", nCent,0,nCent);
	hResoPvRctPsi->Sumw2();
	hResoPvShfPsi->Sumw2();

	h2EastRawQ = new TH2D(mName+"EastRawQ", mName+"EastRawQ", 50,-0.25,10.25, 50,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2EastCorQ = new TH2D(mName+"EastCorQ", mName+"EastCorQ", 50,-0.25,10.25, 50,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2WestRawQ = new TH2D(mName+"WestRawQ", mName+"WestRawQ", 50,-10.25,0.25, 50,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2WestCorQ = new TH2D(mName+"WestCorQ", mName+"WestCorQ", 50,-10.25,0.25, 50,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2FullRawQ = new TH2D(mName+"FullRawQ", mName+"FullRawQ", 50,-0.25,10.25, 50,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2FullCorQ = new TH2D(mName+"FullCorQ", mName+"FullCorQ", 50,-0.25,10.25, 50,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2EastRawQ->Sumw2();
	h2EastCorQ->Sumw2();
	h2WestRawQ->Sumw2();
	h2WestCorQ->Sumw2();
	h2FullRawQ->Sumw2();
	h2FullCorQ->Sumw2();
}


MyZdcSmd::MyZdcSmd() {
	mDir = "";
	mName = "";
	mStep = 1;
	Init();
}


MyZdcSmd::MyZdcSmd(TString name) {
	mDir = "";
	mName = name;
	mStep = 1;
	Init();
}


MyZdcSmd::MyZdcSmd(TString dir, TString name) {
	mDir = dir;
	mName = name;
	mStep = 1;
	Init();
}


MyZdcSmd::~MyZdcSmd() {
}


bool MyZdcSmd::IsGoodSlat(int iEW, int iVH, int iST) {
	if(iEW<0 || iEW>=nEW) return false;

	if(iVH==1) {
		return ((iST>=0) && (iST<nST));
	} else if(iVH==0) {
		return ((iST>=0) && (iST<nST-1));
	} else {
		return false;
	}

	return false;
}


double MyZdcSmd::GetPosition(int iEW, int iVH, int iST) {
	if(!IsGoodSlat(iEW, iVH, iST)) return -999;

	double zdcsmd_x[7] = {0.5,2,3.5,5,6.5,8,9.5};
	double zdcsmd_y[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

	if(iEW==0 && iVH==0) return  zdcsmd_x[iST];
	if(iEW==1 && iVH==0) return -zdcsmd_x[iST];
	if(iEW==0 && iVH==1) return  zdcsmd_y[iST]/sqrt(2.);
	if(iEW==1 && iVH==1) return  zdcsmd_y[iST]/sqrt(2.);
	
	return 0;
}


double MyZdcSmd::GetRctPosition(int iEW, int iVH, int iST) {
	return (GetPosition(iEW, iVH, iST) - CorCenter[iEW][iVH]);
}


double MyZdcSmd::GetPvRctPosition(int iEW, int iVH, int iST) {
	double d = ((iVH==0)?(mPvx-AvePvx):(mPvy-AvePvy));
	return (GetRctPosition(iEW, iVH, iST) - d);
}


double MyZdcSmd::ZdcSmdAdcRawToCor(int iEW, int iVH, int iST, double raw) {
	// Pedestal + Gain correction
	double Pede = mZdcSmdPede[iEW][iVH][iST];
	double Gain = mZdcSmdGain[iEW][iVH][iST];
	double cor = (raw - Pede) / Gain;
	if(cor<0) cor = 0;
	return cor;
}


int MyZdcSmd::GetLastBin(TH1D *h) {
	int n = h->GetXaxis()->GetNbins();
	int m = 0;
	for(int i=n; i>0; i--) {
		if(h->GetBinContent(i)>0) {
			m ++;
		} else {
			m = 0;
		}
		if(m>4) return i;
	}
	return 0;
}


vector<double> MyZdcSmd::GetPede() const {
	vector<double> tmp;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				tmp.push_back(mZdcSmdPede[iEW][iVH][iST]);
			}
		}
	}
	return tmp;
}
void MyZdcSmd::SetPede(vector<double> tmp) {
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				mZdcSmdPede[iEW][iVH][iST] = tmp[iST+nST*(iVH+nVH*iEW)];
			}
		}
	}
}

vector<double> MyZdcSmd::GetGain() const {
	vector<double> tmp;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				tmp.push_back(mZdcSmdGain[iEW][iVH][iST]);
			}
		}
	}
	return tmp;
}
void MyZdcSmd::SetGain(vector<double> tmp) {
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				mZdcSmdGain[iEW][iVH][iST] = tmp[iST+nST*(iVH+nVH*iEW)];
			}
		}
	}
}

vector<double> MyZdcSmd::GetPeak() const {
	vector<double> tmp;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				tmp.push_back(mZdcSmdPeak[iEW][iVH][iST]);
			}
		}
	}
	return tmp;
}
void MyZdcSmd::SetPeak(vector<double> tmp) {
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				mZdcSmdPeak[iEW][iVH][iST] = tmp[iST+nST*(iVH+nVH*iEW)];
			}
		}
	}
}


bool MyZdcSmd::InputFile(TString filepath) {
	mInputFilePath = filepath;
	TFile *fInput = new TFile(filepath, "READ");
	cout << "MyZdcSmd: " << mName << endl;
	if(fInput->IsZombie()) {
		cout << "MyZdcSmd: invalid input file! "  << filepath << endl;
		cout << "MyZdcSmd: go back to step 1!" << endl;
		mStep = 1;
		cout << "MyZdcSmd: mStep = " << mStep << "." << endl;
		return false;
	}

	cout << "MyZdcSmd: file loaded from " << filepath << endl;

	TString StrDir = "";
	if(mDir!="") StrDir = mDir+"/";
	TString StrDirName = StrDir + mName;

	TProfile *hInputAvePvx = (TProfile*)fInput->Get(StrDirName+"AvePvx");
	TProfile *hInputAvePvy = (TProfile*)fInput->Get(StrDirName+"AvePvy");
	if(!hInputAvePvx) {
		cout << "MyZdcSmd: " << StrDirName << " histograms not found!" << endl;
		if(mRefDir=="") StrDirName = mRefName;
		if(mRefDir!="") StrDirName = mRefDir + "/" + mRefName;
		cout << "MyZdcSmd: search backup histograms: " << StrDirName << endl;
		hInputAvePvx = (TProfile*)fInput->Get(StrDirName+"AvePvx");
		hInputAvePvy = (TProfile*)fInput->Get(StrDirName+"AvePvy");
	}
	if(!hInputAvePvx) {
		cout << "MyZdcSmd: invalid input histogram!" << endl;
		mStep = 1;
		cout << "MyZdcSmd: mStep = " << mStep << "." << endl;
		return false;
	}
	if(hInputAvePvx) hInputAvePvx->SetName("Input"+mName+"AvePvx");
	if(hInputAvePvy) hInputAvePvy->SetName("Input"+mName+"AvePvy");
	if(hInputAvePvx) AvePvx = hInputAvePvx->GetBinContent(1);
	if(hInputAvePvy) AvePvy = hInputAvePvy->GetBinContent(1);


	TProfile *hInputZdcSmdPede = (TProfile*)fInput->Get(StrDirName+"ZdcSmdPede");
	TProfile *hInputZdcSmdGain = (TProfile*)fInput->Get(StrDirName+"ZdcSmdGain");
	TProfile *hInputZdcSmdPeak = (TProfile*)fInput->Get(StrDirName+"ZdcSmdPeak");
	TProfile *hInputZdcRawCenter = (TProfile*)fInput->Get(StrDirName+"ZdcRawCenter");
	TProfile *hInputZdcCorCenter = (TProfile*)fInput->Get(StrDirName+"ZdcCorCenter");
	if(hInputZdcSmdPede) hInputZdcSmdPede->SetName("Input"+mName+"ZdcSmdPede");
	if(hInputZdcSmdGain) hInputZdcSmdGain->SetName("Input"+mName+"ZdcSmdGain");
	if(hInputZdcSmdPeak) hInputZdcSmdPeak->SetName("Input"+mName+"ZdcSmdPeak");
	if(hInputZdcRawCenter) hInputZdcRawCenter->SetName("Input"+mName+"ZdcRawCenter");
	if(hInputZdcCorCenter) hInputZdcCorCenter->SetName("Input"+mName+"ZdcCorCenter");
	bool mIsGoodInputGain = false;
	if(hInputZdcSmdGain) if(hInputZdcSmdGain->GetEntries()>0) mIsGoodInputGain = true;
	if(!mIsGoodInputGain) cout << "MyZdcSmd: input rootfile does NOT have gain parameter; mGainMod = " << mGainMod << endl;

	TH1D *hInputZdcSmdRawAdc[nEW][nVH][nST];
	TF1 *funcInputZdcSmdRawAdc[nEW][nVH][nST];
	double SumSlat = 0;
	double SumGain = 0;
	double MaxLast = 0;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				TString StrTmp = Form("EW%dVH%dST%d", iEW, iVH, iST);
				hInputZdcSmdRawAdc[iEW][iVH][iST] = (TH1D*)fInput->Get(StrDirName+"ZdcSmdRawAdc"+StrTmp);
				hInputZdcSmdRawAdc[iEW][iVH][iST]->SetName("Input"+mName+"ZdcSmdRawAdc"+StrTmp);
				if(!IsGoodSlat(iEW, iVH, iST)) continue;
				if(hInputZdcSmdRawAdc[iEW][iVH][iST]->GetEntries()==0) {
					mIsGoodInput = false;
					continue;
				}
				if((mGainMod==0 && !mIsGoodInputGain) || (mGainMod==3)) {
					int bPeak = hInputZdcSmdRawAdc[iEW][iVH][iST]->GetMaximumBin();
					double Peak = hInputZdcSmdRawAdc[iEW][iVH][iST]->GetBinLowEdge(bPeak);
					int bLast = GetLastBin(hInputZdcSmdRawAdc[iEW][iVH][iST]);
					double Last = hInputZdcSmdRawAdc[iEW][iVH][iST]->GetBinLowEdge(bLast);
					//double Tail = 2*Peak;
					//if(bLast>0 && 0.8*Last>2*Peak) Tail = 0.8*Last;
					//if(Tail >= Last) Tail = 0.8*Last;
					double Tail = Last - (Last-Peak)*0.30;
					//funcInputZdcSmdRawAdc[iEW][iVH][iST] = new TF1("funcInput"+mName+"ZdcSmdRawAdc"+StrTmp, "[0]*exp(-(x-[1])/[2])", Tail,0.8*Last>2*Peak?Last:800);
					funcInputZdcSmdRawAdc[iEW][iVH][iST] = new TF1("funcInput"+mName+"ZdcSmdRawAdc"+StrTmp, "[0]*exp(-(x-[1])/[2])", Tail,Last);
					double ValTail = hInputZdcSmdRawAdc[iEW][iVH][iST]->GetBinContent(hInputZdcSmdRawAdc[iEW][iVH][iST]->GetXaxis()->FindBin(Tail));
					double ValLast = hInputZdcSmdRawAdc[iEW][iVH][iST]->GetBinContent(hInputZdcSmdRawAdc[iEW][iVH][iST]->GetXaxis()->FindBin(Last));
					//cout << StrTmp << ": (" << Tail << ", " << ValTail << ")   (" << Last << ", " << ValLast << ")" << endl;
					funcInputZdcSmdRawAdc[iEW][iVH][iST]->SetParameter(0, ValTail);
					funcInputZdcSmdRawAdc[iEW][iVH][iST]->FixParameter(1, Tail);
					funcInputZdcSmdRawAdc[iEW][iVH][iST]->SetParameter(2, -(Tail-Last)/log(ValTail));
					funcInputZdcSmdRawAdc[iEW][iVH][iST]->SetLineColor(kRed);
					funcInputZdcSmdRawAdc[iEW][iVH][iST]->SetLineWidth(5);
					hInputZdcSmdRawAdc[iEW][iVH][iST]->Fit(funcInputZdcSmdRawAdc[iEW][iVH][iST], "R0Q");
					double Gain = funcInputZdcSmdRawAdc[iEW][iVH][iST]->GetParameter(2);
					mZdcSmdPede[iEW][iVH][iST] = 0; // TODO: we just use default pedestal here. I think it has be corrected in data production, so PicoDst does not need this.
					mZdcSmdGain[iEW][iVH][iST] = Gain;
					mZdcSmdPeak[iEW][iVH][iST] = Peak;
					SumSlat += 1.0;
					SumGain += Gain;
					if(MaxLast<Last) MaxLast = Last;
				}
			}
		}
	}
	double AveGain = SumGain / SumSlat;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				if(!IsGoodSlat(iEW, iVH, iST)) continue;
				if((mGainMod==0 && !mIsGoodInputGain) || (mGainMod==3)) mZdcSmdGain[iEW][iVH][iST] /= AveGain;
				if(mGainMod==1) mZdcSmdGain[iEW][iVH][iST] = 1.0;
				if(mGainMod==0 && mIsGoodInputGain) {
					int iAL = iST + nST*(iVH + nVH*iEW);
					mZdcSmdGain[iEW][iVH][iST] = hInputZdcSmdGain->GetBinContent(iAL+1);
				}
			}
		}
	}

	if(mVerbose>=2) {
	cout << "MyZdcSmd: Pedestal & Gain Table" << endl;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				//if(!IsGoodSlat(iEW, iVH, iST)) continue;
				cout << StrEW[iEW] << " " << StrVH[iVH] << " " << "Slat" << iST << " " << flush;
				cout << "Peak=" << mZdcSmdPeak[iEW][iVH][iST] << "   " << flush;
				cout << "Pede=" << mZdcSmdPede[iEW][iVH][iST] << "   " << flush;
				cout << "Gain=" << mZdcSmdGain[iEW][iVH][iST] << endl;
			}
		}
	}
	}

	cout << "MyZdcSmd: Pedestal & Gain done." << endl;

	double WtsRawAdc[nEW][nVH][nST] = {0};
	double WtsCorAdc[nEW][nVH][nST] = {0};
	double SumRawAdc[nEW][nVH][nST] = {0};
	double SumCorAdc[nEW][nVH][nST] = {0};
	double ErrRawAdc[nEW][nVH][nST] = {0};
	double ErrCorAdc[nEW][nVH][nST] = {0};
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				if(!IsGoodSlat(iEW, iVH, iST)) continue;
				for(int i=0; i<hInputZdcSmdRawAdc[iEW][iVH][iST]->GetXaxis()->GetNbins(); i++) {
					double Entry  = hInputZdcSmdRawAdc[iEW][iVH][iST]->GetBinContent(i+1);
					double Error  = hInputZdcSmdRawAdc[iEW][iVH][iST]->GetBinError(i+1);
					double RawAdc = hInputZdcSmdRawAdc[iEW][iVH][iST]->GetXaxis()->GetBinLowEdge(i+1);
					double CorAdc = ZdcSmdAdcRawToCor(iEW, iVH, iST, RawAdc);
					WtsRawAdc[iEW][iVH][iST] += Entry;
					WtsCorAdc[iEW][iVH][iST] += Entry;
					SumRawAdc[iEW][iVH][iST] += Entry*RawAdc;
					SumCorAdc[iEW][iVH][iST] += Entry*CorAdc;
					ErrRawAdc[iEW][iVH][iST] += pow(Error*RawAdc, 2.0);
					ErrCorAdc[iEW][iVH][iST] += pow(Error*CorAdc, 2.0);
				}
			}
		}
	}

	//TGraphErrors *ogZdcRawAdc[nEW][nVH];
	//TGraphErrors *ogZdcCorAdc[nEW][nVH];
	TGraph *ogZdcRawAdc[nEW][nVH];
	TGraph *ogZdcCorAdc[nEW][nVH];
	TF1 *funcZdcRawAdc[nEW][nVH];
	TF1 *funcZdcCorAdc[nEW][nVH];
	double RawGausCenter[nEW][nVH];
	double CorGausCenter[nEW][nVH];
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			TString StrTmp = Form("EW%dVH%d", iEW, iVH);
			int nTmp = iVH==0?nST-1:nST;
			//ogZdcRawAdc[iEW][iVH] = new TGraphErrors(nTmp, Position[iEW][iVH], SumRawAdc[iEW][iVH], Poserror[iEW][iVH], ErrRawAdc[iEW][iVH]);
			//ogZdcCorAdc[iEW][iVH] = new TGraphErrors(nTmp, Position[iEW][iVH], SumCorAdc[iEW][iVH], Poserror[iEW][iVH], ErrCorAdc[iEW][iVH]);
			ogZdcRawAdc[iEW][iVH] = new TGraph(nTmp, Position[iEW][iVH], SumRawAdc[iEW][iVH]);
			ogZdcCorAdc[iEW][iVH] = new TGraph(nTmp, Position[iEW][iVH], SumCorAdc[iEW][iVH]);
			ogZdcRawAdc[iEW][iVH]->Sort();
			ogZdcCorAdc[iEW][iVH]->Sort();
			funcZdcRawAdc[iEW][iVH] = new TF1("func"+mName+"ZdcRawAdc"+StrTmp, "gaus", ogZdcRawAdc[iEW][iVH]->GetX()[0]-0.5, ogZdcRawAdc[iEW][iVH]->GetX()[nTmp-1]+0.5);
			funcZdcCorAdc[iEW][iVH] = new TF1("func"+mName+"ZdcCorAdc"+StrTmp, "gaus", ogZdcCorAdc[iEW][iVH]->GetX()[0]-0.5, ogZdcCorAdc[iEW][iVH]->GetX()[nTmp-1]+0.5);
			funcZdcRawAdc[iEW][iVH]->SetParameter(1, RawCenter[iEW][iVH]);
			funcZdcCorAdc[iEW][iVH]->SetParameter(1, CorCenter[iEW][iVH]);
			ogZdcRawAdc[iEW][iVH]->Fit(funcZdcRawAdc[iEW][iVH], "Q0");
			ogZdcCorAdc[iEW][iVH]->Fit(funcZdcCorAdc[iEW][iVH], "Q0");
			RawGausCenter[iEW][iVH] = funcZdcRawAdc[iEW][iVH]->GetParameter(1);
			CorGausCenter[iEW][iVH] = funcZdcCorAdc[iEW][iVH]->GetParameter(1);
		}
	}

	double WtsRawCenter[nEW][nVH] = {0};
	double WtsCorCenter[nEW][nVH] = {0};
	double SumRawCenter[nEW][nVH] = {0};
	double SumCorCenter[nEW][nVH] = {0};
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				if(!IsGoodSlat(iEW, iVH, iST)) continue;
				double position = GetPosition(iEW, iVH, iST);
				WtsRawCenter[iEW][iVH] += SumRawAdc[iEW][iVH][iST];
				WtsCorCenter[iEW][iVH] += SumCorAdc[iEW][iVH][iST];
				SumRawCenter[iEW][iVH] += SumRawAdc[iEW][iVH][iST]*position;
				SumCorCenter[iEW][iVH] += SumCorAdc[iEW][iVH][iST]*position;
			}
			//RawCenter[iEW][iVH] = SumRawCenter[iEW][iVH] / WtsRawCenter[iEW][iVH];
			//CorCenter[iEW][iVH] = SumCorCenter[iEW][iVH] / WtsCorCenter[iEW][iVH];
			//if(WtsRawCenter[iEW][iVH]<=0) RawCenter[iEW][iVH] = RefCenter[iEW][iVH];
			//if(WtsCorCenter[iEW][iVH]<=0) CorCenter[iEW][iVH] = RefCenter[iEW][iVH];
			if(mCenterMod==0 && mIsGoodInputGain) RawCenter[iEW][iVH] = hInputZdcRawCenter->GetBinContent(iVH+nVH*iEW + 1);
			if(mCenterMod==0 && mIsGoodInputGain) CorCenter[iEW][iVH] = hInputZdcCorCenter->GetBinContent(iVH+nVH*iEW + 1);
			if(mCenterMod==0 && !mIsGoodInputGain && WtsRawCenter[iEW][iVH]> 0) RawCenter[iEW][iVH] = SumRawCenter[iEW][iVH] / WtsRawCenter[iEW][iVH];
			if(mCenterMod==0 && !mIsGoodInputGain && WtsCorCenter[iEW][iVH]> 0) CorCenter[iEW][iVH] = SumCorCenter[iEW][iVH] / WtsCorCenter[iEW][iVH];
			if(mCenterMod==0 && !mIsGoodInputGain && WtsRawCenter[iEW][iVH]<=0) RawCenter[iEW][iVH] = RefCenter[iEW][iVH];
			if(mCenterMod==0 && !mIsGoodInputGain && WtsCorCenter[iEW][iVH]<=0) CorCenter[iEW][iVH] = RefCenter[iEW][iVH];
			if(mCenterMod==1) RawCenter[iEW][iVH] = RefCenter[iEW][iVH];
			if(mCenterMod==1) CorCenter[iEW][iVH] = RefCenter[iEW][iVH];
			if(mCenterMod==2 && mIsGoodInputGain) RawCenter[iEW][iVH] = hInputZdcRawCenter->GetBinContent(iVH+nVH*iEW + 1);
			if(mCenterMod==2 && mIsGoodInputGain) CorCenter[iEW][iVH] = hInputZdcCorCenter->GetBinContent(iVH+nVH*iEW + 1);
			if(mCenterMod==2 && !mIsGoodInputGain) RawCenter[iEW][iVH] = RawGausCenter[iEW][iVH];
			if(mCenterMod==2 && !mIsGoodInputGain) CorCenter[iEW][iVH] = CorGausCenter[iEW][iVH];
		}
	}

	TProfile *hInputRawCenter[nEW][nVH];
	TProfile *hInputCorCenter[nEW][nVH];
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			TString StrTmp = Form("EW%dVH%d", iEW, iVH);
			hInputRawCenter[iEW][iVH] = (TProfile*)fInput->Get(StrDirName+"RawCenter"+StrTmp);
			hInputCorCenter[iEW][iVH] = (TProfile*)fInput->Get(StrDirName+"CorCenter"+StrTmp);
			hInputRawCenter[iEW][iVH]->SetName("Input"+mName+"RawCenter"+StrTmp);
			hInputCorCenter[iEW][iVH]->SetName("Input"+mName+"CorCenter"+StrTmp);
		}
	}

	if(mVerbose>=1) {
	cout << "MyZdcSmd: Center Table" << endl;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			cout << StrEW[iEW] << " " << StrXY[iVH] << " " << flush;
			cout << "Raw " << RawCenter[iEW][iVH] << "   " << hInputRawCenter[iEW][iVH]->GetBinContent(1) << "   " << RawGausCenter[iEW][iVH] << endl;
			cout << StrEW[iEW] << " " << StrXY[iVH] << " " << flush;
			cout << "Cor " << CorCenter[iEW][iVH] << "   " << hInputCorCenter[iEW][iVH]->GetBinContent(1) << "   " << CorGausCenter[iEW][iVH] << endl;
		}
	}
	}

	mStep = 2;
	cout << "MyZdcSmd: Recenter done." << endl;
	cout << "MyZdcSmd: mStep = " << mStep << "." << endl;

	TProfile *hInputEastRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hInputEastRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hInputWestRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hInputWestRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hInputFullRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hInputFullRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hInputEastPvRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hInputEastPvRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hInputWestPvRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hInputWestPvRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hInputFullPvRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hInputFullPvRctPsiSin[nCent]; // Psi=Psi1
	for(int iCent=0; iCent<nCent; iCent++) {
		TString StrCent = Form("Cent%d", iCent);
		hInputEastRctPsiCos[iCent] = (TProfile*)fInput->Get(StrDirName+"EastRctPsiCos"+StrCent);
		hInputEastRctPsiSin[iCent] = (TProfile*)fInput->Get(StrDirName+"EastRctPsiSin"+StrCent);
		hInputWestRctPsiCos[iCent] = (TProfile*)fInput->Get(StrDirName+"WestRctPsiCos"+StrCent);
		hInputWestRctPsiSin[iCent] = (TProfile*)fInput->Get(StrDirName+"WestRctPsiSin"+StrCent);
		hInputFullRctPsiCos[iCent] = (TProfile*)fInput->Get(StrDirName+"FullRctPsiCos"+StrCent);
		hInputFullRctPsiSin[iCent] = (TProfile*)fInput->Get(StrDirName+"FullRctPsiSin"+StrCent);
		if(hInputEastRctPsiCos[iCent]) hInputEastRctPsiCos[iCent]->SetName("Input"+mName+"EastRctPsiCos"+StrCent);
		if(hInputEastRctPsiSin[iCent]) hInputEastRctPsiSin[iCent]->SetName("Input"+mName+"EastRctPsiSin"+StrCent);
		if(hInputWestRctPsiCos[iCent]) hInputWestRctPsiCos[iCent]->SetName("Input"+mName+"WestRctPsiCos"+StrCent);
		if(hInputWestRctPsiSin[iCent]) hInputWestRctPsiSin[iCent]->SetName("Input"+mName+"WestRctPsiSin"+StrCent);
		if(hInputFullRctPsiCos[iCent]) hInputFullRctPsiCos[iCent]->SetName("Input"+mName+"FullRctPsiCos"+StrCent);
		if(hInputFullRctPsiSin[iCent]) hInputFullRctPsiSin[iCent]->SetName("Input"+mName+"FullRctPsiSin"+StrCent);
		hInputEastPvRctPsiCos[iCent] = (TProfile*)fInput->Get(StrDirName+"EastPvRctPsiCos"+StrCent);
		hInputEastPvRctPsiSin[iCent] = (TProfile*)fInput->Get(StrDirName+"EastPvRctPsiSin"+StrCent);
		hInputWestPvRctPsiCos[iCent] = (TProfile*)fInput->Get(StrDirName+"WestPvRctPsiCos"+StrCent);
		hInputWestPvRctPsiSin[iCent] = (TProfile*)fInput->Get(StrDirName+"WestPvRctPsiSin"+StrCent);
		hInputFullPvRctPsiCos[iCent] = (TProfile*)fInput->Get(StrDirName+"FullPvRctPsiCos"+StrCent);
		hInputFullPvRctPsiSin[iCent] = (TProfile*)fInput->Get(StrDirName+"FullPvRctPsiSin"+StrCent);
		if(hInputEastPvRctPsiCos[iCent]) hInputEastPvRctPsiCos[iCent]->SetName("Input"+mName+"EastPvRctPsiCos"+StrCent);
		if(hInputEastPvRctPsiSin[iCent]) hInputEastPvRctPsiSin[iCent]->SetName("Input"+mName+"EastPvRctPsiSin"+StrCent);
		if(hInputWestPvRctPsiCos[iCent]) hInputWestPvRctPsiCos[iCent]->SetName("Input"+mName+"WestPvRctPsiCos"+StrCent);
		if(hInputWestPvRctPsiSin[iCent]) hInputWestPvRctPsiSin[iCent]->SetName("Input"+mName+"WestPvRctPsiSin"+StrCent);
		if(hInputFullPvRctPsiCos[iCent]) hInputFullPvRctPsiCos[iCent]->SetName("Input"+mName+"FullPvRctPsiCos"+StrCent);
		if(hInputFullPvRctPsiSin[iCent]) hInputFullPvRctPsiSin[iCent]->SetName("Input"+mName+"FullPvRctPsiSin"+StrCent);
	}

#ifdef MyZdcSmdPlot
	gStyle->SetOptTitle(1);

	TCanvas *cInputZdcSmdRawAdc[nEW][nVH][nST];
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				if(!IsGoodSlat(iEW, iVH, iST)) continue;
				TString StrTmp = Form("EW%dVH%dST%d", iEW, iVH, iST);
				cInputZdcSmdRawAdc[iEW][iVH][iST] = new TCanvas(mName+"InputZdcSmdRawAdc"+StrTmp, mName+"InputZdcSmdRawAdc"+StrTmp);
				hInputZdcSmdRawAdc[iEW][iVH][iST]->GetXaxis()->SetRangeUser(0,1.1*MaxLast);
				hInputZdcSmdRawAdc[iEW][iVH][iST]->Draw("");
				funcInputZdcSmdRawAdc[iEW][iVH][iST]->Draw("same");
				gPad->SetLogy();
				cInputZdcSmdRawAdc[iEW][iVH][iST]->SaveAs("./plot/"+mName+"InputZdcSmdRawAdc"+StrTmp+".png");
			}
		}
	}

	TH2D *h2ZdcRawAdc[nEW];
	TH2D *h2ZdcCorAdc[nEW];
	h2ZdcRawAdc[0] = new TH2D(mName+"ZdcRawAdcEW0", "ZdcSmd ADCx*ADCy East", 7,-0.25,10.25, 8,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2ZdcRawAdc[1] = new TH2D(mName+"ZdcRawAdcEW1", "ZdcSmd ADCx*ADCy West", 7,-10.25,0.25, 8,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2ZdcCorAdc[0] = new TH2D(mName+"ZdcCorAdcEW0", "ZdcSmd ADCx*ADCy East", 7,-0.25,10.25, 8,0.25/sqrt(2.0),16.25/sqrt(2.0));
	h2ZdcCorAdc[1] = new TH2D(mName+"ZdcCorAdcEW1", "ZdcSmd ADCx*ADCy West", 7,-10.25,0.25, 8,0.25/sqrt(2.0),16.25/sqrt(2.0));

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iST=0; iST<nST; iST++) { // iVH=0; x
			if(!IsGoodSlat(iEW, 0, iST)) continue;
			double x = GetPosition(iEW, 0, iST);
			for(int jST=0; jST<nST; jST++) { // iVH=1; y
				if(!IsGoodSlat(iEW, 1, jST)) continue;
				double y = GetPosition(iEW, 1, jST);
				h2ZdcRawAdc[iEW]->Fill(x, y, SumRawAdc[iEW][0][iST]*SumRawAdc[iEW][1][jST]);
				h2ZdcCorAdc[iEW]->Fill(x, y, SumCorAdc[iEW][0][iST]*SumCorAdc[iEW][1][jST]);
			}
		}
		h2ZdcRawAdc[iEW]->Scale(1.0/h2ZdcRawAdc[iEW]->GetMaximum());
		h2ZdcCorAdc[iEW]->Scale(1.0/h2ZdcCorAdc[iEW]->GetMaximum());
	}

	gStyle->SetOptTitle(1);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadTopMargin(0.10);

	TGraph *ogRawCenter[nEW];
	TGraph *ogCorCenter[nEW];
	TGraph *ogRawGausCenter[nEW];
	TGraph *ogCorGausCenter[nEW];
	for(int iEW=0; iEW<nEW; iEW++) {
		ogRawCenter[iEW] = new TGraph(1, &RawCenter[iEW][0], &RawCenter[iEW][1]);
		ogCorCenter[iEW] = new TGraph(1, &CorCenter[iEW][0], &CorCenter[iEW][1]);
		ogRawCenter[iEW]->SetMarkerStyle(20);
		ogCorCenter[iEW]->SetMarkerStyle(20);
		ogRawGausCenter[iEW] = new TGraph(1, &RawGausCenter[iEW][0], &RawGausCenter[iEW][1]);
		ogCorGausCenter[iEW] = new TGraph(1, &CorGausCenter[iEW][0], &CorGausCenter[iEW][1]);
		ogRawGausCenter[iEW]->SetMarkerStyle(21);
		ogCorGausCenter[iEW]->SetMarkerStyle(21);
		ogRawGausCenter[iEW]->SetMarkerColor(2);
		ogCorGausCenter[iEW]->SetMarkerColor(2);
	}

	TCanvas *c2ZdcRawAdc[nEW];
	TLatex *texZdcRawAdc[nEW];
	for(int iEW=0; iEW<nEW; iEW++) {
		TString StrEW = Form("EW%d", iEW);
		c2ZdcRawAdc[iEW] = new TCanvas("ZdcRawAdc"+StrEW, "ZdcRawAdc"+StrEW, 500,511);
		h2ZdcRawAdc[iEW]->GetXaxis()->SetTitle("x (cm)");
		h2ZdcRawAdc[iEW]->GetYaxis()->SetTitle("y (cm)");
		h2ZdcRawAdc[iEW]->GetYaxis()->SetTitleOffset(0.9);
		h2ZdcRawAdc[iEW]->Draw("colz");
		ogRawCenter[iEW]->Draw("P same");
		//ogRawGausCenter[iEW]->Draw("P same");
		texZdcRawAdc[iEW] = new TLatex();
		texZdcRawAdc[iEW]->DrawLatexNDC(0.80,0.943, "#font[42]{raw}");
		texZdcRawAdc[iEW]->DrawLatex(ogRawCenter[iEW]->GetX()[0]+0.2, ogRawCenter[iEW]->GetY()[0]-0.4, Form("#font[42]{(%.3f, %.3f)}",ogRawCenter[iEW]->GetX()[0], ogRawCenter[iEW]->GetY()[0]));
		//texZdcRawAdc[iEW]->DrawLatex(ogRawGausCenter[iEW]->GetX()[0]-4.5, ogRawGausCenter[iEW]->GetY()[0]+0.2, Form("#color[2]{#font[42]{(%.3f, %.3f)}}",ogRawGausCenter[iEW]->GetX()[0], ogRawGausCenter[iEW]->GetY()[0]));
		c2ZdcRawAdc[iEW]->SaveAs("./plot/"+mName+"ZdcRawAdc"+StrEW+".pdf");
	}

	TCanvas *c2ZdcCorAdc[nEW];
	TLatex *texZdcCorAdc[nEW];
	for(int iEW=0; iEW<nEW; iEW++) {
		TString StrEW = Form("EW%d", iEW);
		c2ZdcCorAdc[iEW] = new TCanvas("ZdcCorAdc"+StrEW, "ZdcCorAdc"+StrEW, 500,511);
		h2ZdcCorAdc[iEW]->GetXaxis()->SetTitle("x (cm)");
		h2ZdcCorAdc[iEW]->GetYaxis()->SetTitle("y (cm)");
		h2ZdcCorAdc[iEW]->GetYaxis()->SetTitleOffset(0.9);
		h2ZdcCorAdc[iEW]->Draw("colz");
		ogCorCenter[iEW]->Draw("P same");
		//ogCorGausCenter[iEW]->Draw("P same");
		texZdcCorAdc[iEW] = new TLatex();
		texZdcCorAdc[iEW]->DrawLatexNDC(0.70,0.943, "#font[42]{#splitline{corrected}{for gain rate}}");
		texZdcCorAdc[iEW]->DrawLatex(ogCorCenter[iEW]->GetX()[0]+0.2, ogCorCenter[iEW]->GetY()[0]-0.4, Form("#font[42]{(%.3f, %.3f)}",ogCorCenter[iEW]->GetX()[0], ogCorCenter[iEW]->GetY()[0]));
		//texZdcCorAdc[iEW]->DrawLatex(ogCorGausCenter[iEW]->GetX()[0]-4.5, ogCorGausCenter[iEW]->GetY()[0]+0.2, Form("#color[2]{#font[42]{(%.3f, %.3f)}}",ogCorGausCenter[iEW]->GetX()[0], ogCorGausCenter[iEW]->GetY()[0]));
		c2ZdcCorAdc[iEW]->SaveAs("./plot/"+mName+"ZdcCorAdc"+StrEW+".pdf");
	}

	TCanvas *cZdcRawAdc[nEW][nVH];
	for(int iEW=0; iEW<nEW; iEW++) {
		TString StrEW = Form("EW%d", iEW);
		for(int iVH=0; iVH<nVH; iVH++) {
			TString StrVH = Form("VH%d", iVH);
			cZdcRawAdc[iEW][iVH] = new TCanvas("ZdcRawAdc"+StrEW+StrVH, "ZdcRawAdc"+StrEW+StrVH);
			ogZdcRawAdc[iEW][iVH]->SetTitle("ZdcRawAdc"+StrEW+StrVH);
			ogZdcRawAdc[iEW][iVH]->SetMarkerStyle(24);
			ogZdcRawAdc[iEW][iVH]->Draw("PLA");
			funcZdcRawAdc[iEW][iVH]->Draw("same");
			cZdcRawAdc[iEW][iVH]->SaveAs("./plot/"+mName+"ZdcRawAdc"+StrEW+StrVH+".pdf");
		}
	}

	TCanvas *cZdcCorAdc[nEW][nVH];
	for(int iEW=0; iEW<nEW; iEW++) {
		TString StrEW = Form("EW%d", iEW);
		for(int iVH=0; iVH<nVH; iVH++) {
			TString StrVH = Form("VH%d", iVH);
			cZdcCorAdc[iEW][iVH] = new TCanvas("ZdcCorAdc"+StrEW+StrVH, "ZdcCorAdc"+StrEW+StrVH);
			ogZdcCorAdc[iEW][iVH]->SetTitle("ZdcCorAdc"+StrEW+StrVH);
			ogZdcCorAdc[iEW][iVH]->SetMarkerStyle(24);
			ogZdcCorAdc[iEW][iVH]->Draw("PLA");
			funcZdcCorAdc[iEW][iVH]->Draw("same");
			cZdcCorAdc[iEW][iVH]->SaveAs("./plot/"+mName+"ZdcCorAdc"+StrEW+StrVH+".pdf");
		}
	}
#endif

	if(!hInputEastRctPsiCos[1]) {
		cout << "MyZdcSmd: Shift not ready!!!" << endl;
		cout << "MyZdcSmd: mStep = " << mStep << "." << endl;
		return true;
	}

	for(int iCent=0; iCent<nCent; iCent++) {
		for(int iShf=1; iShf<nShf; iShf++) {
			// Cos(iShf*Psi), Sin(iShf*Psi)
			EastRctPsiCos[iCent][iShf] = hInputEastRctPsiCos[iCent]->GetBinContent(iShf);
			EastRctPsiSin[iCent][iShf] = hInputEastRctPsiSin[iCent]->GetBinContent(iShf);
			WestRctPsiCos[iCent][iShf] = hInputWestRctPsiCos[iCent]->GetBinContent(iShf);
			WestRctPsiSin[iCent][iShf] = hInputWestRctPsiSin[iCent]->GetBinContent(iShf);
			FullRctPsiCos[iCent][iShf] = hInputFullRctPsiCos[iCent]->GetBinContent(iShf);
			FullRctPsiSin[iCent][iShf] = hInputFullRctPsiSin[iCent]->GetBinContent(iShf);
			EastPvRctPsiCos[iCent][iShf] = hInputEastPvRctPsiCos[iCent]->GetBinContent(iShf);
			EastPvRctPsiSin[iCent][iShf] = hInputEastPvRctPsiSin[iCent]->GetBinContent(iShf);
			WestPvRctPsiCos[iCent][iShf] = hInputWestPvRctPsiCos[iCent]->GetBinContent(iShf);
			WestPvRctPsiSin[iCent][iShf] = hInputWestPvRctPsiSin[iCent]->GetBinContent(iShf);
			FullPvRctPsiCos[iCent][iShf] = hInputFullPvRctPsiCos[iCent]->GetBinContent(iShf);
			FullPvRctPsiSin[iCent][iShf] = hInputFullPvRctPsiSin[iCent]->GetBinContent(iShf);
		}
	}

	mStep = 3;
	cout << "MyZdcSmd: Shift done." << endl;
	cout << "MyZdcSmd: mStep = " << mStep << "." << endl;
	cout << "MyZdcSmd: Ready!" << endl;

	return true;
}


double MyZdcSmd::Phi02Pi(double phi) {
	double tmpphi = phi;
	while(tmpphi<0) {
		tmpphi += 2*M_PI;
	}
	while(tmpphi>=2*M_PI) {
		tmpphi -= 2*M_PI;
	}
	return tmpphi;
}


void MyZdcSmd::CalcGau() {
	TGraph *ogOneZdcRawAdc[nEW][nVH];
	TGraph *ogOneZdcCorAdc[nEW][nVH];
	TF1 *funcOneZdcRawAdc[nEW][nVH];
	TF1 *funcOneZdcCorAdc[nEW][nVH];
	double RawWts[nEW][nVH] = {0};
	double CorWts[nEW][nVH] = {0};
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			TString StrTmp = Form("EW%dVH%d", iEW, iVH);
			int nTmp = iVH==0?nST-1:nST;
			//ogOneZdcRawAdc[iEW][iVH] = new TGraphErrors(nTmp, Position[iEW][iVH], SumRawAdc[iEW][iVH], Poserror[iEW][iVH], ErrRawAdc[iEW][iVH]);
			//ogOneZdcCorAdc[iEW][iVH] = new TGraphErrors(nTmp, Position[iEW][iVH], SumCorAdc[iEW][iVH], Poserror[iEW][iVH], ErrCorAdc[iEW][iVH]);
			ogOneZdcRawAdc[iEW][iVH] = new TGraph(nTmp, Position[iEW][iVH], mZdcSmdRawAdc[iEW][iVH]);
			ogOneZdcCorAdc[iEW][iVH] = new TGraph(nTmp, Position[iEW][iVH], mZdcSmdCorAdc[iEW][iVH]);
			ogOneZdcRawAdc[iEW][iVH]->Sort();
			ogOneZdcCorAdc[iEW][iVH]->Sort();
			funcOneZdcRawAdc[iEW][iVH] = new TF1("func"+mName+"OneZdcRawAdc"+StrTmp, "gaus", ogOneZdcRawAdc[iEW][iVH]->GetX()[0]-0.5, ogOneZdcRawAdc[iEW][iVH]->GetX()[nTmp-1]+0.5);
			funcOneZdcCorAdc[iEW][iVH] = new TF1("func"+mName+"OneZdcCorAdc"+StrTmp, "gaus", ogOneZdcCorAdc[iEW][iVH]->GetX()[0]-0.5, ogOneZdcCorAdc[iEW][iVH]->GetX()[nTmp-1]+0.5);
			funcOneZdcRawAdc[iEW][iVH]->SetParameter(1, RawCenter[iEW][iVH]);
			funcOneZdcCorAdc[iEW][iVH]->SetParameter(1, CorCenter[iEW][iVH]);
			ogOneZdcRawAdc[iEW][iVH]->Fit(funcOneZdcRawAdc[iEW][iVH], "Q0");
			ogOneZdcCorAdc[iEW][iVH]->Fit(funcOneZdcCorAdc[iEW][iVH], "Q0");
			for(int iST=0; iST<nST; iST++) {
				RawWts[iEW][iVH] += mZdcSmdRawAdc[iEW][iVH][iST];
				CorWts[iEW][iVH] += mZdcSmdCorAdc[iEW][iVH][iST];
			}
		}
	}

	mEastRawQ.Set(RawWts[0][0]>0?funcOneZdcRawAdc[0][0]->GetParameter(1):0, RawWts[0][1]>0?funcOneZdcRawAdc[0][1]->GetParameter(1):0);
	mWestRawQ.Set(RawWts[1][0]>0?funcOneZdcRawAdc[1][0]->GetParameter(1):0, RawWts[1][1]>0?funcOneZdcRawAdc[1][1]->GetParameter(1):0);
	mFullRawQ = 0.5*(mEastRawQ - mWestRawQ);
	mEastCorQ.Set(CorWts[0][0]>0?funcOneZdcCorAdc[0][0]->GetParameter(1):0, CorWts[0][1]>0?funcOneZdcCorAdc[0][1]->GetParameter(1):0);
	mWestCorQ.Set(CorWts[1][0]>0?funcOneZdcCorAdc[1][0]->GetParameter(1):0, CorWts[1][1]>0?funcOneZdcCorAdc[1][1]->GetParameter(1):0);
	mFullCorQ = 0.5*(mEastCorQ - mWestCorQ);
	mEastRctQ = mEastCorQ - TVector2(CorCenter[0][0], CorCenter[0][1]);
	mWestRctQ = mWestCorQ - TVector2(CorCenter[1][0], CorCenter[1][1]);
	mFullRctQ = 0.5*(mEastRctQ - mWestRctQ);
	mEastPvRctQ = mEastRctQ - TVector2(mPvx-AvePvx, mPvy-AvePvy);
	mWestPvRctQ = mWestRctQ - TVector2(mPvx-AvePvx, mPvy-AvePvy);
	mFullPvRctQ = 0.5*(mEastPvRctQ - mWestPvRctQ);

	mIsGoodEvent = mIsGoodInput;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			if(RawWts[iEW][iVH]<=0 || CorWts[iEW][iVH]<=0) mIsGoodEvent = false;
#ifndef MyZdcSmdPlot
			delete ogOneZdcRawAdc[iEW][iVH];
			delete ogOneZdcCorAdc[iEW][iVH];
			delete funcOneZdcRawAdc[iEW][iVH];
			delete funcOneZdcCorAdc[iEW][iVH];
#endif
		}
	}

#ifdef MyZdcSmdPlot
	TCanvas *cOneZdcCorAdc[nEW][nVH];
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			TString StrOneZdcCorAdc = Form("EW%dVH%d", iEW, iVH);
			cOneZdcCorAdc[iEW][iVH] = new TCanvas("OneZdcCorAdc"+StrOneZdcCorAdc, "OneZdcCorAdc"+StrOneZdcCorAdc);
			ogOneZdcCorAdc[iEW][iVH]->Draw("PLA");
			funcOneZdcCorAdc[iEW][iVH]->Draw("same");
		}
	}
#endif

	return;
}


void MyZdcSmd::CalcAve() {
	double WtsEastRawX = 0;
	double SumEastRawX = 0;
	double WtsEastCorX = 0;
	double SumEastCorX = 0;
	double WtsEastRctX = 0;
	double SumEastRctX = 0;
	double WtsEastPvRctX = 0;
	double SumEastPvRctX = 0;
	for(int iST=0; iST<nST; iST++) { // East, X
		int iEW = 0; // East
		int iVH = 0; // X
		if(!IsGoodSlat(iEW, iVH, iST)) continue;
		double position = GetPosition(iEW, iVH, iST);
		WtsEastRawX += mZdcSmdRawAdc[iEW][iVH][iST];
		SumEastRawX += mZdcSmdRawAdc[iEW][iVH][iST]*position;
		WtsEastCorX += mZdcSmdCorAdc[iEW][iVH][iST];
		SumEastCorX += mZdcSmdCorAdc[iEW][iVH][iST]*position;
		double rctposition = GetRctPosition(iEW, iVH, iST);
		WtsEastRctX += mZdcSmdCorAdc[iEW][iVH][iST];
		SumEastRctX += mZdcSmdCorAdc[iEW][iVH][iST]*rctposition;
		double pvrctposition = GetPvRctPosition(iEW, iVH, iST);
		WtsEastPvRctX += mZdcSmdCorAdc[iEW][iVH][iST];
		SumEastPvRctX += mZdcSmdCorAdc[iEW][iVH][iST]*pvrctposition;
	}

	double WtsEastRawY = 0;
	double SumEastRawY = 0;
	double WtsEastCorY = 0;
	double SumEastCorY = 0;
	double WtsEastRctY = 0;
	double SumEastRctY = 0;
	double WtsEastPvRctY = 0;
	double SumEastPvRctY = 0;
	for(int iST=0; iST<nST; iST++) { // East, Y
		int iEW = 0; // East
		int iVH = 1; // Y
		if(!IsGoodSlat(iEW, iVH, iST)) continue;
		double position = GetPosition(iEW, iVH, iST);
		WtsEastRawY += mZdcSmdRawAdc[iEW][iVH][iST];
		SumEastRawY += mZdcSmdRawAdc[iEW][iVH][iST]*position;
		WtsEastCorY += mZdcSmdCorAdc[iEW][iVH][iST];
		SumEastCorY += mZdcSmdCorAdc[iEW][iVH][iST]*position;
		double rctposition = GetRctPosition(iEW, iVH, iST);
		WtsEastRctY += mZdcSmdCorAdc[iEW][iVH][iST];
		SumEastRctY += mZdcSmdCorAdc[iEW][iVH][iST]*rctposition;
		double pvrctposition = GetPvRctPosition(iEW, iVH, iST);
		WtsEastPvRctY += mZdcSmdCorAdc[iEW][iVH][iST];
		SumEastPvRctY += mZdcSmdCorAdc[iEW][iVH][iST]*pvrctposition;
	}

	double WtsWestRawX = 0;
	double SumWestRawX = 0;
	double WtsWestCorX = 0;
	double SumWestCorX = 0;
	double WtsWestRctX = 0;
	double SumWestRctX = 0;
	double WtsWestPvRctX = 0;
	double SumWestPvRctX = 0;
	for(int iST=0; iST<nST; iST++) { // West, X
		int iEW = 1; // West
		int iVH = 0; // X
		if(!IsGoodSlat(iEW, iVH, iST)) continue;
		double position = GetPosition(iEW, iVH, iST);
		WtsWestRawX += mZdcSmdRawAdc[iEW][iVH][iST];
		SumWestRawX += mZdcSmdRawAdc[iEW][iVH][iST]*position;
		WtsWestCorX += mZdcSmdCorAdc[iEW][iVH][iST];
		SumWestCorX += mZdcSmdCorAdc[iEW][iVH][iST]*position;
		double rctposition = GetRctPosition(iEW, iVH, iST);
		WtsWestRctX += mZdcSmdCorAdc[iEW][iVH][iST];
		SumWestRctX += mZdcSmdCorAdc[iEW][iVH][iST]*rctposition;
		double pvrctposition = GetPvRctPosition(iEW, iVH, iST);
		WtsWestPvRctX += mZdcSmdCorAdc[iEW][iVH][iST];
		SumWestPvRctX += mZdcSmdCorAdc[iEW][iVH][iST]*pvrctposition;
	}

	double WtsWestRawY = 0;
	double SumWestRawY = 0;
	double WtsWestCorY = 0;
	double SumWestCorY = 0;
	double WtsWestRctY = 0;
	double SumWestRctY = 0;
	double WtsWestPvRctY = 0;
	double SumWestPvRctY = 0;
	for(int iST=0; iST<nST; iST++) { // West, Y
		int iEW = 1; // West
		int iVH = 1; // Y
		if(!IsGoodSlat(iEW, iVH, iST)) continue;
		double position = GetPosition(iEW, iVH, iST);
		WtsWestRawY += mZdcSmdRawAdc[iEW][iVH][iST];
		SumWestRawY += mZdcSmdRawAdc[iEW][iVH][iST]*position;
		WtsWestCorY += mZdcSmdCorAdc[iEW][iVH][iST];
		SumWestCorY += mZdcSmdCorAdc[iEW][iVH][iST]*position;
		double rctposition = GetRctPosition(iEW, iVH, iST);
		WtsWestRctY += mZdcSmdCorAdc[iEW][iVH][iST];
		SumWestRctY += mZdcSmdCorAdc[iEW][iVH][iST]*rctposition;
		double pvrctposition = GetPvRctPosition(iEW, iVH, iST);
		WtsWestPvRctY += mZdcSmdCorAdc[iEW][iVH][iST];
		SumWestPvRctY += mZdcSmdCorAdc[iEW][iVH][iST]*pvrctposition;
	}

	double WtsFullRawX = WtsEastRawX + WtsWestRawX;
	double SumFullRawX = SumEastRawX - SumWestRawX;
	double WtsFullRawY = WtsEastRawY + WtsWestRawY;
	double SumFullRawY = SumEastRawY - SumWestRawY;
	double WtsFullCorX = WtsEastCorX + WtsWestCorX;
	double SumFullCorX = SumEastCorX - SumWestCorX;
	double WtsFullCorY = WtsEastCorY + WtsWestCorY;
	double SumFullCorY = SumEastCorY - SumWestCorY;
	double WtsFullRctX = WtsEastRctX + WtsWestRctX;
	double SumFullRctX = SumEastRctX - SumWestRctX;
	double WtsFullRctY = WtsEastRctY + WtsWestRctY;
	double SumFullRctY = SumEastRctY - SumWestRctY;
	double WtsFullPvRctX = WtsEastPvRctX + WtsWestPvRctX;
	double SumFullPvRctX = SumEastPvRctX - SumWestPvRctX;
	double WtsFullPvRctY = WtsEastPvRctY + WtsWestPvRctY;
	double SumFullPvRctY = SumEastPvRctY - SumWestPvRctY;

	mEastRawQ.Set(WtsEastRawX>0?SumEastRawX/WtsEastRawX:0, WtsEastRawY>0?SumEastRawY/WtsEastRawY:0);
	mWestRawQ.Set(WtsWestRawX>0?SumWestRawX/WtsWestRawX:0, WtsWestRawY>0?SumWestRawY/WtsWestRawY:0);
	mFullRawQ.Set(WtsFullRawX>0?SumFullRawX/WtsFullRawX:0, WtsFullRawY>0?SumFullRawY/WtsFullRawY:0);
	if(mMerge==1) mFullRawQ = 0.5*(mEastRawQ - mWestRawQ);
	mEastCorQ.Set(WtsEastCorX>0?SumEastCorX/WtsEastCorX:0, WtsEastCorY>0?SumEastCorY/WtsEastCorY:0);
	mWestCorQ.Set(WtsWestCorX>0?SumWestCorX/WtsWestCorX:0, WtsWestCorY>0?SumWestCorY/WtsWestCorY:0);
	mFullCorQ.Set(WtsFullCorX>0?SumFullCorX/WtsFullCorX:0, WtsFullCorY>0?SumFullCorY/WtsFullCorY:0);
	if(mMerge==1) mFullCorQ = 0.5*(mEastCorQ - mWestCorQ);
	mEastRctQ.Set(WtsEastRctX>0?SumEastRctX/WtsEastRctX:0, WtsEastRctY>0?SumEastRctY/WtsEastRctY:0);
	mWestRctQ.Set(WtsWestRctX>0?SumWestRctX/WtsWestRctX:0, WtsWestRctY>0?SumWestRctY/WtsWestRctY:0);
	mFullRctQ.Set(WtsFullRctX>0?SumFullRctX/WtsFullRctX:0, WtsFullRctY>0?SumFullRctY/WtsFullRctY:0);
	if(mMerge==1) mFullRctQ = 0.5*(mEastRctQ - mWestRctQ);
	mEastPvRctQ.Set(WtsEastPvRctX>0?SumEastPvRctX/WtsEastPvRctX:0, WtsEastPvRctY>0?SumEastPvRctY/WtsEastPvRctY:0);
	mWestPvRctQ.Set(WtsWestPvRctX>0?SumWestPvRctX/WtsWestPvRctX:0, WtsWestPvRctY>0?SumWestPvRctY/WtsWestPvRctY:0);
	mFullPvRctQ.Set(WtsFullPvRctX>0?SumFullPvRctX/WtsFullPvRctX:0, WtsFullPvRctY>0?SumFullPvRctY/WtsFullPvRctY:0);
	if(mMerge==1) mFullPvRctQ = 0.5*(mEastPvRctQ - mWestPvRctQ);

	mIsGoodEvent = mIsGoodInput;
	if(WtsEastRawX<=0 || WtsEastRawY<=0) mIsGoodEvent = false;
	if(WtsWestRawX<=0 || WtsWestRawY<=0) mIsGoodEvent = false;
	if(WtsEastCorX<=0 || WtsEastCorY<=0) mIsGoodEvent = false;
	if(WtsWestCorX<=0 || WtsWestCorY<=0) mIsGoodEvent = false;
	if(WtsEastRctX<=0 || WtsEastRctY<=0) mIsGoodEvent = false;
	if(WtsWestRctX<=0 || WtsWestRctY<=0) mIsGoodEvent = false;
	if(WtsEastPvRctX<=0 || WtsEastPvRctY<=0) mIsGoodEvent = false;
	if(WtsWestPvRctX<=0 || WtsWestPvRctY<=0) mIsGoodEvent = false;

	return;
}


void MyZdcSmd::CalcPsi() {
	mEastRawPsi = Phi02Pi(mEastRawQ.Phi());
	mWestRawPsi = Phi02Pi(mWestRawQ.Phi());
	mFullRawPsi = Phi02Pi(mFullRawQ.Phi());
	mEastCorPsi = Phi02Pi(mEastCorQ.Phi());
	mWestCorPsi = Phi02Pi(mWestCorQ.Phi());
	mFullCorPsi = Phi02Pi(mFullCorQ.Phi());
	mEastRctPsi = Phi02Pi(mEastRctQ.Phi());
	mWestRctPsi = Phi02Pi(mWestRctQ.Phi());
	mFullRctPsi = Phi02Pi(mFullRctQ.Phi());
	mEastPvRctPsi = Phi02Pi(mEastPvRctQ.Phi());
	mWestPvRctPsi = Phi02Pi(mWestPvRctQ.Phi());
	mFullPvRctPsi = Phi02Pi(mFullPvRctQ.Phi());

	mEastShfPsi = mEastRctPsi;
	for(int iShf=1; iShf<nShf; iShf++) mEastShfPsi += 2.0/iShf*(-EastRctPsiSin[mCent][iShf]*cos(iShf*mEastRctPsi) + EastRctPsiCos[mCent][iShf]*sin(iShf*mEastRctPsi));
	mEastShfPsi = Phi02Pi(mEastShfPsi);
	mEastShfQ.Set(mEastRctQ.Mod()*cos(mEastShfPsi), mEastRctQ.Mod()*sin(mEastShfPsi));

	mWestShfPsi = mWestRctPsi;
	for(int iShf=1; iShf<nShf; iShf++) mWestShfPsi += 2.0/iShf*(-WestRctPsiSin[mCent][iShf]*cos(iShf*mWestRctPsi) + WestRctPsiCos[mCent][iShf]*sin(iShf*mWestRctPsi));
	mWestShfPsi = Phi02Pi(mWestShfPsi);
	mWestShfQ.Set(mWestRctQ.Mod()*cos(mWestShfPsi), mWestRctQ.Mod()*sin(mWestShfPsi));

	mFullShfPsi = mFullRctPsi;
	for(int iShf=1; iShf<nShf; iShf++) mFullShfPsi += 2.0/iShf*(-FullRctPsiSin[mCent][iShf]*cos(iShf*mFullRctPsi) + FullRctPsiCos[mCent][iShf]*sin(iShf*mFullRctPsi));
	mFullShfPsi = Phi02Pi(mFullShfPsi);
	mFullShfQ.Set(mFullRctQ.Mod()*cos(mFullShfPsi), mFullRctQ.Mod()*sin(mFullShfPsi));

	mEastPvShfPsi = mEastPvRctPsi;
	for(int iShf=1; iShf<nShf; iShf++) mEastPvShfPsi += 2.0/iShf*(-EastPvRctPsiSin[mCent][iShf]*cos(iShf*mEastPvRctPsi) + EastPvRctPsiCos[mCent][iShf]*sin(iShf*mEastPvRctPsi));
	mEastPvShfPsi = Phi02Pi(mEastPvShfPsi);
	mEastPvShfQ.Set(mEastPvRctQ.Mod()*cos(mEastPvShfPsi), mEastPvRctQ.Mod()*sin(mEastPvShfPsi));

	mWestPvShfPsi = mWestPvRctPsi;
	for(int iShf=1; iShf<nShf; iShf++) mWestPvShfPsi += 2.0/iShf*(-WestPvRctPsiSin[mCent][iShf]*cos(iShf*mWestPvRctPsi) + WestPvRctPsiCos[mCent][iShf]*sin(iShf*mWestPvRctPsi));
	mWestPvShfPsi = Phi02Pi(mWestPvShfPsi);
	mWestPvShfQ.Set(mWestPvRctQ.Mod()*cos(mWestPvShfPsi), mWestPvRctQ.Mod()*sin(mWestPvShfPsi));

	mFullPvShfPsi = mFullPvRctPsi;
	for(int iShf=1; iShf<nShf; iShf++) mFullPvShfPsi += 2.0/iShf*(-FullPvRctPsiSin[mCent][iShf]*cos(iShf*mFullPvRctPsi) + FullPvRctPsiCos[mCent][iShf]*sin(iShf*mFullPvRctPsi));
	mFullPvShfPsi = Phi02Pi(mFullPvShfPsi);
	mFullPvShfQ.Set(mFullPvRctQ.Mod()*cos(mFullPvShfPsi), mFullPvRctQ.Mod()*sin(mFullPvShfPsi));

	return;
}


void MyZdcSmd::Calc() {
	if(mCenterMod==2) {
		CalcGau();
	} else {
		CalcAve();
	}
	CalcPsi();

	return;
}


void MyZdcSmd::Fill() {
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				int iAL = iST + nST*(iVH + nVH*iEW); 
				hZdcSmdPede->Fill(iAL, mZdcSmdPede[iEW][iVH][iST], mWght);
				hZdcSmdGain->Fill(iAL, mZdcSmdGain[iEW][iVH][iST], mWght);
				hZdcSmdPeak->Fill(iAL, mZdcSmdPeak[iEW][iVH][iST], mWght);
			}
			hZdcRawCenter->Fill(iVH+nVH*iEW, RawCenter[iEW][iVH], mWght);
			hZdcCorCenter->Fill(iVH+nVH*iEW, CorCenter[iEW][iVH], mWght);
		}
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				hZdcSmdRawAdc[iEW][iVH][iST]->Fill(mZdcSmdRawAdc[iEW][iVH][iST], mWght);
				//hZdcSmdCorAdc[iEW][iVH][iST]->Fill(mZdcSmdCorAdc[iEW][iVH][iST], mWght);
			}
		}
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				if(!IsGoodSlat(iEW, iVH, iST)) continue;
				double position = GetPosition(iEW, iVH, iST);
				hRawCenter[iEW][iVH]->Fill(0.5, position, mZdcSmdRawAdc[iEW][iVH][iST]*mWght);
				hCorCenter[iEW][iVH]->Fill(0.5, position, mZdcSmdCorAdc[iEW][iVH][iST]*mWght);
			}
		}
	}

	hAvePvx->Fill(0.5, mPvx, mWght);
	hAvePvy->Fill(0.5, mPvy, mWght);

	for(int iShf=1; iShf<nShf; iShf++) {
		hEastRctPsiCos[mCent]->Fill(iShf, cos(iShf*mEastRctPsi), mWght);
		hEastRctPsiSin[mCent]->Fill(iShf, sin(iShf*mEastRctPsi), mWght);
		hWestRctPsiCos[mCent]->Fill(iShf, cos(iShf*mWestRctPsi), mWght);
		hWestRctPsiSin[mCent]->Fill(iShf, sin(iShf*mWestRctPsi), mWght);
		hFullRctPsiCos[mCent]->Fill(iShf, cos(iShf*mFullRctPsi), mWght);
		hFullRctPsiSin[mCent]->Fill(iShf, sin(iShf*mFullRctPsi), mWght);
		hEastPvRctPsiCos[mCent]->Fill(iShf, cos(iShf*mEastPvRctPsi), mWght);
		hEastPvRctPsiSin[mCent]->Fill(iShf, sin(iShf*mEastPvRctPsi), mWght);
		hWestPvRctPsiCos[mCent]->Fill(iShf, cos(iShf*mWestPvRctPsi), mWght);
		hWestPvRctPsiSin[mCent]->Fill(iShf, sin(iShf*mWestPvRctPsi), mWght);
		hFullPvRctPsiCos[mCent]->Fill(iShf, cos(iShf*mFullPvRctPsi), mWght);
		hFullPvRctPsiSin[mCent]->Fill(iShf, sin(iShf*mFullPvRctPsi), mWght);
	}

	hEastRawPsi[mCent]->Fill(mEastRawPsi, mWght);
	hEastCorPsi[mCent]->Fill(mEastCorPsi, mWght);
	hEastRctPsi[mCent]->Fill(mEastRctPsi, mWght);
	hEastShfPsi[mCent]->Fill(mEastShfPsi, mWght);
	hWestRawPsi[mCent]->Fill(mWestRawPsi, mWght);
	hWestCorPsi[mCent]->Fill(mWestCorPsi, mWght);
	hWestRctPsi[mCent]->Fill(mWestRctPsi, mWght);
	hWestShfPsi[mCent]->Fill(mWestShfPsi, mWght);
	hFullRawPsi[mCent]->Fill(mFullRawPsi, mWght);
	hFullCorPsi[mCent]->Fill(mFullCorPsi, mWght);
	hFullRctPsi[mCent]->Fill(mFullRctPsi, mWght);
	hFullShfPsi[mCent]->Fill(mFullShfPsi, mWght);
	h2SubRawPsi[mCent]->Fill(mEastRawPsi, mWestRawPsi, mWght);
	h2SubCorPsi[mCent]->Fill(mEastCorPsi, mWestCorPsi, mWght);
	h2SubRctPsi[mCent]->Fill(mEastRctPsi, mWestRctPsi, mWght);
	h2SubShfPsi[mCent]->Fill(mEastShfPsi, mWestShfPsi, mWght);
	hEastPvRctPsi[mCent]->Fill(mEastPvRctPsi, mWght);
	hEastPvShfPsi[mCent]->Fill(mEastPvShfPsi, mWght);
	hWestPvRctPsi[mCent]->Fill(mWestPvRctPsi, mWght);
	hWestPvShfPsi[mCent]->Fill(mWestPvShfPsi, mWght);
	hFullPvRctPsi[mCent]->Fill(mFullPvRctPsi, mWght);
	hFullPvShfPsi[mCent]->Fill(mFullPvShfPsi, mWght);
	h2SubPvRctPsi[mCent]->Fill(mEastPvRctPsi, mWestPvRctPsi, mWght);
	h2SubPvShfPsi[mCent]->Fill(mEastPvShfPsi, mWestPvShfPsi, mWght);

	hResoRawPsi->Fill(mCent, cos(mEastRawPsi-mWestRawPsi+M_PI), mWght);
	hResoCorPsi->Fill(mCent, cos(mEastCorPsi-mWestCorPsi+M_PI), mWght);
	hResoRctPsi->Fill(mCent, cos(mEastRctPsi-mWestRctPsi+M_PI), mWght);
	hResoShfPsi->Fill(mCent, cos(mEastShfPsi-mWestShfPsi+M_PI), mWght);
	hResoPvRctPsi->Fill(mCent, cos(mEastPvRctPsi-mWestPvRctPsi+M_PI), mWght);
	hResoPvShfPsi->Fill(mCent, cos(mEastPvShfPsi-mWestPvShfPsi+M_PI), mWght);

	h2EastRawQ->Fill(mEastRawQ.X(), mEastRawQ.Y(), mWght);
	h2EastCorQ->Fill(mEastCorQ.X(), mEastCorQ.Y(), mWght);
	h2WestRawQ->Fill(mWestRawQ.X(), mWestRawQ.Y(), mWght);
	h2WestCorQ->Fill(mWestCorQ.X(), mWestCorQ.Y(), mWght);
	h2FullRawQ->Fill(mFullRawQ.X(), mFullRawQ.Y(), mWght);
	h2FullCorQ->Fill(mFullCorQ.X(), mFullCorQ.Y(), mWght);
}


void MyZdcSmd::InitEvent(int cent, double *zdcsmdadc, TVector3 pv, double wght=1.0) {
	ResetEvent();
	mPv = pv;
	mPvx = mPv.X();
	mPvy = mPv.Y();

	mCent = cent;
	mWght = wght;

	// input EastHor 0-7; EastVer, 8-15; WestHor 16-23; WestVer, 24-31;
	// East 0, West 1; horizontal 1, vertical 0
	// iST = strip-1
	for(int iST=0; iST<nST; iST++) {
		mZdcSmdRawAdc[0][1][iST] = zdcsmdadc[iST];      //  EastHori
		mZdcSmdRawAdc[0][0][iST] = zdcsmdadc[iST+8];    //  EastVert
		mZdcSmdRawAdc[1][1][iST] = zdcsmdadc[iST+16];   //  WestHori
		mZdcSmdRawAdc[1][0][iST] = zdcsmdadc[iST+24];   //  WestVert
	}
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				if(!IsGoodSlat(iEW, iVH, iST)) continue;
				mZdcSmdCorAdc[iEW][iVH][iST] = ZdcSmdAdcRawToCor(iEW, iVH, iST, mZdcSmdRawAdc[iEW][iVH][iST]);
			}
		}
	}

	Calc();
	if(mIsGoodEvent && mIsAutoFill) Fill();
}

void MyZdcSmd::InitEvent(int cent, float *zdcsmdadc, TVector3 pv, double wght=1.0) {
	double tmp[nAL];
	for(int iAL=0; iAL<nAL; iAL++) {
		tmp[iAL] = zdcsmdadc[iAL];
	}
	InitEvent(cent, tmp, pv, wght);
}


#ifndef MyZdcSmdPlot
void MyZdcSmd::InitEvent(int cent, StPicoEvent *event, TVector3 pv, double wght=1.0) {
	double zdcsmdadc[nAL];
	for(int iST=0; iST<nST; iST++){
		zdcsmdadc[iST]   = 1.* event->ZdcSmdEastHorizontal(iST);   //  EastHori  // picoDst function i 0-7
		zdcsmdadc[iST+8] = 1.* event->ZdcSmdEastVertical(iST);     //  EastVert
		zdcsmdadc[iST+16]= 1.* event->ZdcSmdWestHorizontal(iST);   //  WestHori
		zdcsmdadc[iST+24]= 1.* event->ZdcSmdWestVertical(iST);     //  WestVert
	}
	InitEvent(cent, zdcsmdadc, pv, wght);
}
#endif


void MyZdcSmd::InitEvent(int cent, double *zdcsmdadc, double wght=1.0) {
	TVector3 pv(0.0, 0.0, 0.0);
	InitEvent(cent, zdcsmdadc, pv, wght);
}


void MyZdcSmd::InitEvent(int cent, float  *zdcsmdadc, double wght=1.0) {
	TVector3 pv(0.0, 0.0, 0.0);
	InitEvent(cent, zdcsmdadc, pv, wght);
}


#ifndef MyZdcSmdPlot
void MyZdcSmd::InitEvent(int cent, StPicoEvent *event, double wght=1.0) {
	TVector3 pv(0.0, 0.0, 0.0);
	InitEvent(cent, event, pv, wght);
}
#endif


void MyZdcSmd::Finish(TFile *fOut) {
	fOut->cd();
	if(mDir!="") {
		TDirectory *dOut = fOut->mkdir(mDir);
		dOut->cd();
	}

	if(mStep > 1) {
		hZdcSmdPede->Write();
		hZdcSmdGain->Write();
		hZdcSmdPeak->Write();
		hZdcRawCenter->Write();
		hZdcCorCenter->Write();
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			for(int iST=0; iST<nST; iST++) {
				hZdcSmdRawAdc[iEW][iVH][iST]->Write();
				//hZdcSmdCorAdc[iEW][iVH][iST]->Write();
			}
		}
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iVH=0; iVH<nVH; iVH++) {
			hRawCenter[iEW][iVH]->Write();
			hCorCenter[iEW][iVH]->Write();
		}
	}

	hAvePvx->Write();
	hAvePvy->Write();

	if(mStep > 1) {
		for(int iCent=0; iCent<nCent; iCent++) {
			hEastRctPsiCos[iCent]->Write();
			hEastRctPsiSin[iCent]->Write();
			hWestRctPsiCos[iCent]->Write();
			hWestRctPsiSin[iCent]->Write();
			hFullRctPsiCos[iCent]->Write();
			hFullRctPsiSin[iCent]->Write();
			hEastPvRctPsiCos[iCent]->Write();
			hEastPvRctPsiSin[iCent]->Write();
			hWestPvRctPsiCos[iCent]->Write();
			hWestPvRctPsiSin[iCent]->Write();
			hFullPvRctPsiCos[iCent]->Write();
			hFullPvRctPsiSin[iCent]->Write();
		}
	}

	for(int iCent=0; iCent<nCent; iCent++) {
		hEastRawPsi[iCent]->Write();
		hEastCorPsi[iCent]->Write();
		hEastRctPsi[iCent]->Write();
		hEastShfPsi[iCent]->Write();
		hWestRawPsi[iCent]->Write();
		hWestCorPsi[iCent]->Write();
		hWestRctPsi[iCent]->Write();
		hWestShfPsi[iCent]->Write();
		hFullRawPsi[iCent]->Write();
		hFullCorPsi[iCent]->Write();
		hFullRctPsi[iCent]->Write();
		hFullShfPsi[iCent]->Write();
		h2SubRawPsi[iCent]->Write();
		h2SubCorPsi[iCent]->Write();
		h2SubRctPsi[iCent]->Write();
		h2SubShfPsi[iCent]->Write();
		hEastPvRctPsi[iCent]->Write();
		hEastPvShfPsi[iCent]->Write();
		hWestPvRctPsi[iCent]->Write();
		hWestPvShfPsi[iCent]->Write();
		hFullPvRctPsi[iCent]->Write();
		hFullPvShfPsi[iCent]->Write();
		h2SubPvRctPsi[iCent]->Write();
		h2SubPvShfPsi[iCent]->Write();
	}

	hResoRawPsi->Write();
	hResoCorPsi->Write();
	hResoRctPsi->Write();
	hResoShfPsi->Write();
	hResoPvRctPsi->Write();
	hResoPvShfPsi->Write();

	h2EastRawQ->Write();
	h2EastCorQ->Write();
	h2WestRawQ->Write();
	h2WestCorQ->Write();
	h2FullRawQ->Write();
	h2FullCorQ->Write();

	fOut->cd();
}
