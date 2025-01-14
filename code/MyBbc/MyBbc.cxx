#include "MyBbc.h"
#include <TDirectoryFile.h>
#include <TRandom.h>

using namespace std;


void MyBbc::ResetEvent(bool isinit = false) {
	mCent = 0;
	mWght = 1.0;
	mIsGoodEvent = true;

	// [East,West][Tiles]
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			mBbcRawAdc[iEW][iTL] = 0;
			mBbcCorAdc[iEW][iTL] = 0;
		}
	}

	for(int iCent=0; iCent<nCent; iCent++) {
		if(mRawBbcEp[iCent] && !isinit) mRawBbcEp[iCent]->ResetEvent();
		if(mCorBbcEp[iCent] && !isinit) mCorBbcEp[iCent]->ResetEvent();
	}
}


void MyBbc::Init() {
	ResetEvent(true);

	for(int iCent=0; iCent<nCent; iCent++) {
		TString StrCent = Form("Cent%d", iCent);
		mRawBbcEp[iCent] = new MyPlanePack(mName+"RawBbcEp"+StrCent);
		mCorBbcEp[iCent] = new MyPlanePack(mName+"CorBbcEp"+StrCent);
		mRawBbcEp[iCent]->SetIsAutoFill(false);
		mCorBbcEp[iCent]->SetIsAutoFill(false);
	}

	mVerbose = 2;
	mGainMod = 0;
	mIsGoodInput = true;
	mIsAutoFill = true; // default true; if false, must call Fill() later

	// record the parameters from the first calculation, keep consistency
	hBbcPede = new TProfile(mName+"BbcPede", mName+"BbcPede", nAL,0,nAL);
	hBbcGain = new TProfile(mName+"BbcGain", mName+"BbcGain", nAL,0,nAL);
	hBbcPeak = new TProfile(mName+"BbcPeak", mName+"BbcPeak", nAL,0,nAL);
	hBbcPede->Sumw2();
	hBbcGain->Sumw2();
	hBbcPeak->Sumw2();

	// [East,West][Vert,Hori][Tiles]
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			mBbcPede[iEW][iTL] = 0;
			mBbcGain[iEW][iTL] = 1.0;
			mBbcPeak[iEW][iTL] = 0;
		}
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			TString StrTmp = Form("EW%dTL%d", iEW, iTL);
			hBbcRawAdc[iEW][iTL] = new TH1D(mName+"BbcRawAdc"+StrTmp, mName+"BbcRawAdc"+StrTmp, 5000,-5,4995);
			//hBbcCorAdc[iEW][iTL] = new TH1D(mName+"BbcCorAdc"+StrTmp, mName+"BbcCorAdc"+StrTmp, 5000,-5,4995);
			hBbcRawAdc[iEW][iTL]->Sumw2();
			//hBbcCorAdc[iEW][iTL]->Sumw2();
		}
	}

	for(int i=1; i<_NMAX+1; i++) {
		hResoRawPsi[i] = new TProfile(mName+Form("BbcResoRawPsi%d",i), mName+Form("BbcResoRawPsi%d",i), nCent,0,nCent);
		hResoCorPsi[i] = new TProfile(mName+Form("BbcResoCorPsi%d",i), mName+Form("BbcResoCorPsi%d",i), nCent,0,nCent);
		hResoShfPsi[i] = new TProfile(mName+Form("BbcResoShfPsi%d",i), mName+Form("BbcResoShfPsi%d",i), nCent,0,nCent);
		hResoRawPsi[i]->Sumw2();
		hResoCorPsi[i]->Sumw2();
		hResoShfPsi[i]->Sumw2();
	}
}


MyBbc::MyBbc() {
	mDir = "";
	mName = "";
	mStep = 1;
	Init();
}


MyBbc::MyBbc(TString name) {
	mDir = "";
	mName = name;
	mStep = 1;
	Init();
}


MyBbc::MyBbc(TString dir, TString name) {
	mDir = dir;
	mName = name;
	mStep = 1;
	Init();
}


MyBbc::~MyBbc() {
}


bool MyBbc::IsGoodTile(int iEW, int iTL) {
	if(iEW<0 || iEW>=nEW) return false;
	if(iTL<0 || iTL>=nTL) return false;
	return true;
}


double MyBbc::GetAngle(int iEW, int iTL) {
	int eastwest = iEW;
	int iTile = iTL;

	const float phi_div= M_PI/6; 
	const float phi_pi = M_PI; 

	float bbc_phi=phi_div;
	switch(iTile) {
		case 0:  bbc_phi=3*phi_div;
			 break;
		case 1:  bbc_phi=phi_div;
			 break;
		case 2:  bbc_phi=-1*phi_div;
			 break;
		case 3:  bbc_phi=-3*phi_div;
			 break;
		case 4:  bbc_phi=-5*phi_div;
			 break;
		case 5:  bbc_phi=5*phi_div;
			 break;
		case 6:  bbc_phi= (gRandom->Rndm()>0.5) ? 2*phi_div:4*phi_div;
			 break;
		case 7:  bbc_phi=3*phi_div;
			 break;
		case 8:  bbc_phi=phi_div;
			 break;
		case 9:  bbc_phi=0.;
			 break;
		case 10: bbc_phi=-phi_div;
			 break;
		case 11: bbc_phi=(gRandom->Rndm()>0.5) ? -2*phi_div:-4*phi_div;
			 break;
		case 12: bbc_phi=-3*phi_div;
			 break;
		case 13: bbc_phi=-5*phi_div;
			 break;
		case 14: bbc_phi=phi_pi;
			 break;
		case 15: bbc_phi=5*phi_div;
			 break;
	}

	if(eastwest==0){if (bbc_phi > -0.001){ bbc_phi = phi_pi-bbc_phi;}
		else {bbc_phi= -phi_pi-bbc_phi;}
	}

	if(bbc_phi<0.      ) bbc_phi +=2*phi_pi;
	if(bbc_phi>2*phi_pi) bbc_phi -=2*phi_pi;

	return bbc_phi;
}


double MyBbc::BbcAdcRawToCor(int iEW, int iTL, double raw) {
	// Pedestal + Gain correction
	double Pede = mBbcPede[iEW][iTL];
	double Gain = mBbcGain[iEW][iTL];
	double cor = (raw - Pede) / Gain;
	if(cor<0) cor = 0;
	return cor;
}


bool MyBbc::InputFile(TString filepath) {
	mInputFilePath = filepath;
	TFile *fInput = new TFile(filepath, "READ");
	cout << "MyBbc: " << mName << endl;
	if(fInput->IsZombie()) {
		cout << "MyBbc: invalid input file! "  << filepath << endl;
		cout << "MyBbc: go back to step 1!" << endl;
		mStep = 1;
		cout << "MyBbc: mStep = " << mStep << "." << endl;
		return false;
	}

	cout << "MyBbc: file loaded from " << filepath << endl;

	TString StrDir = "";
	if(mDir!="") StrDir = mDir+"/";
	TString StrDirName = StrDir + mName;

	TDirectoryFile *dInput = (TDirectoryFile*)fInput->Get(mDir);
	if(!dInput || dInput->IsZombie()) {
		cout << "MyBbc: " << StrDirName << " histograms not found!" << endl;
		if(mRefDir=="") StrDirName = mRefName;
		if(mRefDir!="") StrDirName = mRefDir + "/" + mRefName;
		cout << "MyBbc: search backup histograms: " << StrDirName << endl;
		dInput = (TDirectoryFile*)fInput->Get(mRefDir);
	}
	if(!dInput || dInput->IsZombie()) {
		cout << "MyBbc: invalid input histogram!" << endl;
		mStep = 1;
		cout << "MyBbc: mStep = " << mStep << "." << endl;
		return false;
	}

	TProfile *hInputBbcPede = (TProfile*)fInput->Get(StrDirName+"BbcPede");
	TProfile *hInputBbcGain = (TProfile*)fInput->Get(StrDirName+"BbcGain");
	TProfile *hInputBbcPeak = (TProfile*)fInput->Get(StrDirName+"BbcPeak");
	if(hInputBbcPede) hInputBbcPede->SetName("Input"+mName+"BbcPede");
	if(hInputBbcGain) hInputBbcGain->SetName("Input"+mName+"BbcGain");
	if(hInputBbcPeak) hInputBbcPeak->SetName("Input"+mName+"BbcPeak");
	bool mIsGoodInputGain = false;
	if(hInputBbcGain) if(hInputBbcGain->GetEntries()>0) mIsGoodInputGain = true;
	if(!mIsGoodInputGain) cout << "MyBbc: input rootfile does NOT have gain parameter; mGainMod = " << mGainMod << endl;

	TH1D *hInputBbcRawAdc[nEW][nTL];
	double SumTile = 0;
	double SumGain = 0;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			TString StrTmp = Form("EW%dTL%d", iEW, iTL);
			hInputBbcRawAdc[iEW][iTL] = (TH1D*)fInput->Get(StrDirName+"BbcRawAdc"+StrTmp);
			hInputBbcRawAdc[iEW][iTL]->SetName("Input"+mName+"BbcRawAdc"+StrTmp);
			if(!IsGoodTile(iEW, iTL)) continue;
			if(hInputBbcRawAdc[iEW][iTL]->GetEntries()==0) {
				mIsGoodInput = false;
				continue;
			}
			if((mGainMod==0 && !mIsGoodInputGain) || (mGainMod==3)) {
				int bPeak = hInputBbcRawAdc[iEW][iTL]->GetMaximumBin();
				double Peak = hInputBbcRawAdc[iEW][iTL]->GetBinLowEdge(bPeak);
				double Gain = hInputBbcRawAdc[iEW][iTL]->GetMean();
				mBbcPede[iEW][iTL] = 0; // TODO: we just use default pedestal here. I think it has be corrected in data production, so PicoDst does not need this.
				mBbcGain[iEW][iTL] = Gain;
				mBbcPeak[iEW][iTL] = Peak;
				SumTile += 1.0;
				SumGain += Gain;
			}
		}
	}
	for(int iEW=0; iEW<nEW; iEW++) {
		double x1=mBbcGain[iEW][0];
		double x2=mBbcGain[iEW][1];
		double x3=mBbcGain[iEW][2];
		double x4=mBbcGain[iEW][3];
		double x5=mBbcGain[iEW][4];
		double x6=mBbcGain[iEW][5];
		double mean1=(x1+x2+x3+x4+x5+x6)/6.0;
		double x7 =mBbcGain[iEW][6] /2.0;
		double x8 =mBbcGain[iEW][6] /2.0;
		double x9 =mBbcGain[iEW][9] ;
		double x10=mBbcGain[iEW][11]/2.0;
		double x11=mBbcGain[iEW][11]/2.0;
		double x12=mBbcGain[iEW][14];
		double mean2=(x7+x8+x9+x10+x11+x12)/6.0;
		double x13=mBbcGain[iEW][7] ;
		double x14=mBbcGain[iEW][8] ;
		double x15=mBbcGain[iEW][10];
		double x16=mBbcGain[iEW][12];
		double x17=mBbcGain[iEW][13];
		double x18=mBbcGain[iEW][15];
		double mean3=(x13+x14+x15+x16+x17+x18)/6.0;
		double tmpgain[16] = {x1/mean1,x2/mean1,x3/mean1,x4/mean1,x5/mean1,x6/mean1,x7/mean2,x13/mean3,x14/mean3,x9/mean2,x15/mean3,x10/mean2,x16/mean3,x17/mean3,x12/mean2,x18/mean3};
		for(int iTL=0; iTL<nTL; iTL++) mBbcGain[iEW][iTL] = tmpgain[iTL];
	}
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			if(!IsGoodTile(iEW, iTL)) continue;
			//if((mGainMod==0 && !mIsGoodInputGain) || (mGainMod==3)) mBbcGain[iEW][iTL]; // as it is
			if(mGainMod==1) mBbcGain[iEW][iTL] = 1.0;
			if(mGainMod==0 && mIsGoodInputGain) {
				int iAL = iTL + nTL*iEW;
				mBbcGain[iEW][iTL] = hInputBbcGain->GetBinContent(iAL+1);
			}
		}
	}

	if(mVerbose>=2) {
	cout << "MyBbc: Pedestal & Gain Table" << endl;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			cout << StrEW[iEW] << " " << "Tile" << iTL << " " << flush;
			cout << "Peak=" << mBbcPeak[iEW][iTL] << "   " << flush;
			cout << "Pede=" << mBbcPede[iEW][iTL] << "   " << flush;
			cout << "Gain=" << mBbcGain[iEW][iTL] << endl;
		}
	}
	}

	cout << "MyBbc: Pedestal & Gain done." << endl;

	double WtsRawAdc[nEW][nTL] = {0};
	double WtsCorAdc[nEW][nTL] = {0};
	double SumRawAdc[nEW][nTL] = {0};
	double SumCorAdc[nEW][nTL] = {0};
	double ErrRawAdc[nEW][nTL] = {0};
	double ErrCorAdc[nEW][nTL] = {0};
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			if(!IsGoodTile(iEW, iTL)) continue;
			for(int i=0; i<hInputBbcRawAdc[iEW][iTL]->GetXaxis()->GetNbins(); i++) {
				double Entry  = hInputBbcRawAdc[iEW][iTL]->GetBinContent(i+1);
				double Error  = hInputBbcRawAdc[iEW][iTL]->GetBinError(i+1);
				double RawAdc = hInputBbcRawAdc[iEW][iTL]->GetXaxis()->GetBinLowEdge(i+1);
				double CorAdc = BbcAdcRawToCor(iEW, iTL, RawAdc);
				WtsRawAdc[iEW][iTL] += Entry;
				WtsCorAdc[iEW][iTL] += Entry;
				SumRawAdc[iEW][iTL] += Entry*RawAdc;
				SumCorAdc[iEW][iTL] += Entry*CorAdc;
				ErrRawAdc[iEW][iTL] += pow(Error*RawAdc, 2.0);
				ErrCorAdc[iEW][iTL] += pow(Error*CorAdc, 2.0);
			}
		}
	}

	mStep = 2;
	cout << "MyBbc: mStep = " << mStep << "." << endl;

	bool isshifted = true;
	for(int iCent=0; iCent<nCent; iCent++) isshifted = isshifted && (mCorBbcEp[iCent]->InputHistogram(dInput));
	if(!isshifted) return true;

	mStep = 3;
	//cout << "MyBbc: Shift done." << endl;
	cout << "MyBbc: mStep = " << mStep << "." << endl;
	cout << "MyBbc: Ready!" << endl;

	return true;
}


void MyBbc::Calc() {
	mRawBbcEp[mCent]->SetWeight(mWght);
	mCorBbcEp[mCent]->SetWeight(mWght);

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			//int iAL = iTL + nTL*iEW; 
			double eta = iEW-0.5; // a placeholder of eta: East<0, West>0
			double phi = GetAngle(iEW, iTL);
			mRawBbcEp[mCent]->Add(phi, eta, mBbcRawAdc[iEW][iTL]);
			mCorBbcEp[mCent]->Add(phi, eta, mBbcCorAdc[iEW][iTL]);
		}
	}

	mIsGoodEvent = mIsGoodInput;

	mRawBbcEp[mCent]->Calc();
	mCorBbcEp[mCent]->Calc();

	mIsGoodEvent = mIsGoodEvent && mRawBbcEp[mCent]->GetIsGood() && mCorBbcEp[mCent]->GetIsGood();
}


double MyBbc::GetFullRawPsi(int ord) {
	return mRawBbcEp[mCent]->Full.GetRawPsi(ord);
}
double MyBbc::GetFullCorPsi(int ord) {
	return mCorBbcEp[mCent]->Full.GetRawPsi(ord);
}
double MyBbc::GetFullShfPsi(int ord) {
	return mCorBbcEp[mCent]->Full.GetCorPsi(ord);
}

double MyBbc::GetEastRawPsi(int ord) {
	return mRawBbcEp[mCent]->East.GetRawPsi(ord);
}
double MyBbc::GetEastCorPsi(int ord) {
	return mCorBbcEp[mCent]->East.GetRawPsi(ord);
}
double MyBbc::GetEastShfPsi(int ord) {
	return mCorBbcEp[mCent]->East.GetCorPsi(ord);
}

double MyBbc::GetWestRawPsi(int ord) {
	return mRawBbcEp[mCent]->West.GetRawPsi(ord);
}
double MyBbc::GetWestCorPsi(int ord) {
	return mCorBbcEp[mCent]->West.GetRawPsi(ord);
}
double MyBbc::GetWestShfPsi(int ord) {
	return mCorBbcEp[mCent]->West.GetCorPsi(ord);
}


void MyBbc::Fill() {
	if(!mIsGoodEvent) return;

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			int iAL = iTL + nTL*iEW; 
			hBbcPede->Fill(iAL, mBbcPede[iEW][iTL], mWght);
			hBbcGain->Fill(iAL, mBbcGain[iEW][iTL], mWght);
			hBbcPeak->Fill(iAL, mBbcPeak[iEW][iTL], mWght);
		}
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			hBbcRawAdc[iEW][iTL]->Fill(mBbcRawAdc[iEW][iTL], mWght);
			//hBbcCorAdc[iEW][iTL]->Fill(mBbcCorAdc[iEW][iTL], mWght);
		}
	}

	mRawBbcEp[mCent]->Fill();
	mCorBbcEp[mCent]->Fill();

	for(int i=1; i<_NMAX+1; i++) {
		// looks mismatched due to different meaning of "Cor" here and MyPlane
		hResoRawPsi[i]->Fill(mCent, mRawBbcEp[mCent]->GetRawReso(i), mWght);
		hResoCorPsi[i]->Fill(mCent, mCorBbcEp[mCent]->GetRawReso(i), mWght);
		hResoShfPsi[i]->Fill(mCent, mCorBbcEp[mCent]->GetCorReso(i), mWght);
	}
}


void MyBbc::Finish(TFile *fOut) {
	fOut->cd();
	if(mDir!="") {
		TDirectory *dOut = fOut->mkdir(mDir);
		dOut->cd();
	}

	if(mStep > 1) {
		hBbcPede->Write();
		hBbcGain->Write();
		hBbcPeak->Write();
	}

	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			hBbcRawAdc[iEW][iTL]->Write();
			//hBbcCorAdc[iEW][iTL]->Write();
		}
	}

	for(int i=1; i<_NMAX+1; i++) {
		hResoRawPsi[i]->Write();
		hResoCorPsi[i]->Write();
		hResoShfPsi[i]->Write();
	}

	if(mStep>1) {
		for(int iCent=0; iCent<nCent; iCent++) {
			mRawBbcEp[iCent]->Finish();
			mCorBbcEp[iCent]->Finish();
		}
	}

	fOut->cd();
}


void MyBbc::InitEvent(int cent, double *bbcadc, double wght) {
	ResetEvent();
	mCent = cent;
	mWght = wght;

	// input East 0-15; West 16-31;
	// East 0, West 1; 
	// [East,West][Tiles]
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			int iAL = iTL + nTL*iEW;
			mBbcRawAdc[iEW][iTL] = bbcadc[iAL];
			mBbcCorAdc[iEW][iTL] = BbcAdcRawToCor(iEW, iTL, mBbcRawAdc[iEW][iTL]);
		}
	}

	Calc();
	if(mIsGoodEvent && mIsAutoFill) Fill();
}


void MyBbc::InitEvent(int cent, float *bbcadc, double wght) {
	double tmpbbcadc[nAL];
	for(int iAL=0; iAL<nAL; iAL++) {
		tmpbbcadc[iAL] = bbcadc[iAL];
	}
	InitEvent(cent, tmpbbcadc, wght);
}

#ifndef MyBbcPlot
void MyBbc::InitEvent(int cent, StPicoEvent *event, double wght) {
	double tmpbbcadc[nAL];
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			int iAL = iTL + nTL*iEW;
			if(iEW==0) tmpbbcadc[iAL] = event->bbcAdcEast(iTL);
			if(iEW==1) tmpbbcadc[iAL] = event->bbcAdcWest(iTL);
		}
	}
	InitEvent(cent, tmpbbcadc, wght);
}
#endif


vector<double> MyBbc::GetPede() const {
	vector<double> tmp;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			tmp.push_back(mBbcPede[iEW][iTL]);
		}
	}
	return tmp;
}
void MyBbc::SetPede(vector<double> tmp) {
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			mBbcPede[iEW][iTL] = tmp[iTL+nTL*iEW];
		}
	}
}

vector<double> MyBbc::GetGain() const {
	vector<double> tmp;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			tmp.push_back(mBbcGain[iEW][iTL]);
		}
	}
	return tmp;
}
void MyBbc::SetGain(vector<double> tmp) {
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			mBbcGain[iEW][iTL] = tmp[iTL+nTL*iEW];
		}
	}
}

vector<double> MyBbc::GetPeak() const {
	vector<double> tmp;
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			tmp.push_back(mBbcPeak[iEW][iTL]);
		}
	}
	return tmp;
}
void MyBbc::SetPeak(vector<double> tmp) {
	for(int iEW=0; iEW<nEW; iEW++) {
		for(int iTL=0; iTL<nTL; iTL++) {
			mBbcPeak[iEW][iTL] = tmp[iTL+nTL*iEW];
		}
	}
}


