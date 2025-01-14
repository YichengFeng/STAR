#include "MyPhiCorrection.h"
#include <cmath>
#include <TMath.h>

using namespace std;


//-----------------------------------------------------------------------------
void MyPhiCorrection::Init(TString Name) { // Please only use 0~9, a~Z, _
	// input
	mName = Name;
	mEventWeight = 1.0;
	//vector<double> TmpPtEdge{0.0,0.4,0.6,0.8,1.0,1.4,2.0};
	vector<double> TmpPtEdge{0.0, 0.4, 0.7};
	mPtEdge = TmpPtEdge;
	mEtaL = -1.8; // TODO: iTPC
	mEtaH = +1.8; // TODO: iTPC
	mPhiL = 0.0;
	mPhiH = 2*M_PI;

	// index
	iCent = 0;
	iEventType = 0;
	iCh = 0;
	iPt = 0;
	iEta = 0;
	iPhi = 0;
	iTrackType = 0;
	iEvtrkType = 0;

	// average
	IsAveKept = true;

	// weight cutoff
	IsWeightCutOff = true;
	mWtL = 0.10;
	mWtH = 10.0;

	// step
	CurrentStage = 1;

	// Vz align
	funcVzDEta = TF1(mName+"FuncVzDEta", "0");
	mDEta = 0;

	// weight
	for(int jEventType=0; jEventType<nEventType; jEventType++) {
		for(int jTrackType=0; jTrackType<nTrackType; jTrackType++) {
			for(int jPhi=0; jPhi<nPhi; jPhi++) {
			for(int jEta=0; jEta<nEta; jEta++) {
				PhiEtaWeight[jEventType][jTrackType][jPhi][jEta] = 1.0;
			} // jEta
			} // jPhi
		} // jTrackType
	} // jEventType

	for(int jFlow=0; jFlow<nFlow; jFlow++) {
		for(int jEvtrkType=0; jEvtrkType<nEvtrkType; jEvtrkType++) {
			RawEvtrkSumWts[jFlow][jEvtrkType] = 0.0;
			RawEvtrkAveCos[jFlow][jEvtrkType] = 0.0;
			RawEvtrkAveSin[jFlow][jEvtrkType] = 0.0;
			CorEvtrkSumWts[jFlow][jEvtrkType] = 0.0;
			CorEvtrkAveCos[jFlow][jEvtrkType] = 0.0;
			CorEvtrkAveSin[jFlow][jEvtrkType] = 0.0;
		} // jEvtrkType
		for(int jEventType=0; jEventType<nEventType; jEventType++) {
			for(int jCh=0; jCh<nCh; jCh++) {
				RawEvchgSumWts[jFlow][jEventType][jCh] = 0.0;
				RawEvchgAveCos[jFlow][jEventType][jCh] = 0.0;
				RawEvchgAveSin[jFlow][jEventType][jCh] = 0.0;
				CorEvchgSumWts[jFlow][jEventType][jCh] = 0.0;
				CorEvchgAveCos[jFlow][jEventType][jCh] = 0.0;
				CorEvchgAveSin[jFlow][jEventType][jCh] = 0.0;
			} // jCh
		} // jEventType
		for(int jCent=0; jCent<nCent; jCent++) {
			for(int jCh=0; jCh<nCh; jCh++) {
				RawCechgSumWts[jFlow][jCent][jCh] = 0.0;
				RawCechgAveCos[jFlow][jCent][jCh] = 0.0;
				RawCechgAveSin[jFlow][jCent][jCh] = 0.0;
				CorCechgSumWts[jFlow][jCent][jCh] = 0.0;
				CorCechgAveCos[jFlow][jCent][jCh] = 0.0;
				CorCechgAveSin[jFlow][jCent][jCh] = 0.0;
			} // jCh
		} // jCent
	} // jFlow

	for(int jCent=0; jCent<nCent; jCent++) {
		for(int jCh=0; jCh<nCh; jCh++) {
			for(int jShf=0; jShf<nShf; jShf++) {
				ShfAveCos[jCent][jCh][jShf] = 0;
				ShfAveSin[jCent][jCh][jShf] = 0;
			} // jShf
		} // jCh
	} // jCent

	// histograms/profiles
	for(int jEventType=0; jEventType<nEventType; jEventType++) {
		for(int jTrackType=0; jTrackType<nTrackType; jTrackType++) {
			TString StrType = Form("EventType%dTrackType%d", jEventType, jTrackType);
			hRawPhiEtaDist[jEventType][jTrackType] = new TH2F(mName+"RawPhiEtaDist"+StrType, mName+"RawPhiEtaDist"+StrType, nPhi,mPhiL,mPhiH, nEta,mEtaL,mEtaH);
			hCorPhiEtaDist[jEventType][jTrackType] = new TH2F(mName+"CorPhiEtaDist"+StrType, mName+"CorPhiEtaDist"+StrType, nPhi,mPhiL,mPhiH, nEta,mEtaL,mEtaH);
		} // jTrackType
	} // jEventType
	if(IsAveKept) {
		for(int jFlow=0; jFlow<nFlow; jFlow++) {
			TString StrFlow = Form("Flow%d", jFlow);
			hRawEvtrkAveCos[jFlow] = new TProfile(mName+"RawEvtrkAveCos"+StrFlow, mName+"RawEvtrkAveCos"+StrFlow, nEvtrkType,0,nEvtrkType);
			hRawEvtrkAveSin[jFlow] = new TProfile(mName+"RawEvtrkAveSin"+StrFlow, mName+"RawEvtrkAveSin"+StrFlow, nEvtrkType,0,nEvtrkType);
			hCorEvtrkAveCos[jFlow] = new TProfile(mName+"CorEvtrkAveCos"+StrFlow, mName+"CorEvtrkAveCos"+StrFlow, nEvtrkType,0,nEvtrkType);
			hCorEvtrkAveSin[jFlow] = new TProfile(mName+"CorEvtrkAveSin"+StrFlow, mName+"CorEvtrkAveSin"+StrFlow, nEvtrkType,0,nEvtrkType);
			hRawEvtrkAveCos[jFlow]->Sumw2();
			hRawEvtrkAveSin[jFlow]->Sumw2();
			hCorEvtrkAveCos[jFlow]->Sumw2();
			hCorEvtrkAveSin[jFlow]->Sumw2();
		} // jFlow
		for(int jCent=0; jCent<nCent; jCent++) {
			TString StrCent = Form("Cent%d", jCent);
			for(int jCh=0; jCh<nCh; jCh++) {
				TString StrCh = Form("Ch%d", jCh);
				hShfAveCos[jCent][jCh] = new TProfile(mName+"ShfAveCos"+StrCent+StrCh, mName+"ShfAveCos"+StrCent+StrCh, nShf,0.5,nShf+0.5);
				hShfAveSin[jCent][jCh] = new TProfile(mName+"ShfAveSin"+StrCent+StrCh, mName+"ShfAveSin"+StrCent+StrCh, nShf,0.5,nShf+0.5);
				hShfAveCos[jCent][jCh]->Sumw2();
				hShfAveSin[jCent][jCh]->Sumw2();
			}
		}
	}

	cout << "MyPhiCorrection: " << mName << endl;
}


//-----------------------------------------------------------------------------
MyPhiCorrection::MyPhiCorrection() {
	mDir = "";
	mRefDir = "";
	mRefName = "";
	Init("");
}


//-----------------------------------------------------------------------------
MyPhiCorrection::MyPhiCorrection(TString Name) {
	mDir = "";
	mRefDir = "";
	mRefName = "";
	Init(Name);
}


//-----------------------------------------------------------------------------
MyPhiCorrection::MyPhiCorrection(TString Dir, TString Name) {
	mDir = Dir;
	mRefDir = "";
	mRefName = "";
	Init(Name);
}


//-----------------------------------------------------------------------------
MyPhiCorrection::~MyPhiCorrection() {
}


//-----------------------------------------------------------------------------
bool MyPhiCorrection::LoadPhiEtaWeight(TString InputRootFile) {
	mPathPhiEtaWeight = InputRootFile;

	CurrentStage = 3;

	TFile *fInputPhiEtaWeight = new TFile(InputRootFile, "READ");
	//TFile *fInputPhiEtaWeight = TFile::Open(InputRootFile);
	if(fInputPhiEtaWeight->IsZombie()) {
		cout << "Invalid fInputPhiEtaWeight!!! " << flush;
		cout << "Go back to step 1." << endl;
		CurrentStage = 1;
		return false;
	}

	cout << "MyPhiCorrection: fInputPhiEtaWeight loaded." << endl;

	TString StrDir = "";
	if(mDir!="") StrDir = mDir + "/";
	TString StrDirName = StrDir + mName;

	TH2F *hInputRawPhiEtaDist[nEventType][nTrackType];
	for(int jEventType=0; jEventType<nEventType; jEventType++) {
		for(int jTrackType=0; jTrackType<nTrackType; jTrackType++) {
			TString StrType = Form("EventType%dTrackType%d", jEventType, jTrackType);
			hInputRawPhiEtaDist[jEventType][jTrackType] = (TH2F*)fInputPhiEtaWeight->Get(StrDirName+"RawPhiEtaDist"+StrType);
			if(!hInputRawPhiEtaDist[jEventType][jTrackType]) {
				cout << "Invalid hInputRawPhiEtaDist!!! " << StrDirName << endl;
				if(mRefDir=="") StrDirName = mRefName;
				if(mRefDir!="") StrDirName = mRefDir + "/" + mRefName;
				cout << "Search backup hInputRawPhiEtaDist: " << StrDirName << endl;
				hInputRawPhiEtaDist[jEventType][jTrackType] = (TH2F*)fInputPhiEtaWeight->Get(StrDirName+"RawPhiEtaDist"+StrType);
				if(!hInputRawPhiEtaDist[jEventType][jTrackType]) {
					cout << "No valid input histogram!!! " << flush;
					cout << "Go back to step 1." << endl;
					CurrentStage = 1;
					return false;
				}
			}
			hInputRawPhiEtaDist[jEventType][jTrackType]->SetName(mName+"InputRawPhiEtaDist"+StrType);
		} // jTrackType
	} // jEventType

	cout << "MyPhiCorrection: RawPhiEtaDist loaded." << endl;

	// ReWeight
	for(int jEventType=0; jEventType<nEventType; jEventType++) {
		for(int jTrackType=0; jTrackType<nTrackType; jTrackType++) {
			for(int jEta=0; jEta<nEta; jEta++) {
				double SumEntry = 0;
				double BinEntry = 0;
				for(int jPhi=0; jPhi<nPhi; jPhi++) {
					double Entry = hInputRawPhiEtaDist[jEventType][jTrackType]->GetBinContent(jPhi+1, jEta+1);
					if(Entry>0) {
						SumEntry += Entry;
						BinEntry += 1.0;
					}
				} // jPhi
				double AveEntry = BinEntry>0?(SumEntry/BinEntry):1.0;
				for(int jPhi=0; jPhi<nPhi; jPhi++) {
					double Entry = hInputRawPhiEtaDist[jEventType][jTrackType]->GetBinContent(jPhi+1, jEta+1);
					double Weight = Entry>0?(AveEntry/Entry):0.0;
					if(IsWeightCutOff && Weight>mWtH) Weight = 0; // detector hole set to 0
#ifdef CheckPhiCorPlot
					//if(Weight==0 && SumEntry!=0 && !(jEta==19 && (jPhi>=106 && jPhi<=116)) && !(jEta==0 && (jPhi>=106 && jPhi<=116))) cout << mName << "   " << jEventType << "   " << jTrackType << "   " << jEta << "   " << jPhi << endl;
					if(Weight==0 && SumEntry!=0) cout << mName << "   " << jEventType << "   " << jTrackType << "   " << jEta << "   " << jPhi << endl;
#endif
					PhiEtaWeight[jEventType][jTrackType][jPhi][jEta] = Weight;
				} // jPhi
			} // jEta
		} // jTrackType
	} // jEventType

	cout << "MyPhiCorrection: PhiEtaWeight calculated." << endl;

	if(IsAveKept) {
		TProfile *hInputRawEvtrkAveCos[nFlow];
		TProfile *hInputRawEvtrkAveSin[nFlow];
		TProfile *hInputCorEvtrkAveCos[nFlow];
		TProfile *hInputCorEvtrkAveSin[nFlow];
		for(int jFlow=0; jFlow<nFlow; jFlow++) {
			TString StrFlow = Form("Flow%d", jFlow);
			hInputRawEvtrkAveCos[jFlow] = (TProfile*)fInputPhiEtaWeight->Get(StrDirName+"RawEvtrkAveCos"+StrFlow);
			hInputRawEvtrkAveSin[jFlow] = (TProfile*)fInputPhiEtaWeight->Get(StrDirName+"RawEvtrkAveSin"+StrFlow);
			hInputCorEvtrkAveCos[jFlow] = (TProfile*)fInputPhiEtaWeight->Get(StrDirName+"CorEvtrkAveCos"+StrFlow);
			hInputCorEvtrkAveSin[jFlow] = (TProfile*)fInputPhiEtaWeight->Get(StrDirName+"CorEvtrkAveSin"+StrFlow);
			if(hInputRawEvtrkAveCos[jFlow]) hInputRawEvtrkAveCos[jFlow]->SetName(mName+"InputRawEvtrkAveCos"+StrFlow);
			if(hInputRawEvtrkAveSin[jFlow]) hInputRawEvtrkAveSin[jFlow]->SetName(mName+"InputRawEvtrkAveSin"+StrFlow);
			if(hInputCorEvtrkAveCos[jFlow]) hInputCorEvtrkAveCos[jFlow]->SetName(mName+"InputCorEvtrkAveCos"+StrFlow);
			if(hInputCorEvtrkAveSin[jFlow]) hInputCorEvtrkAveSin[jFlow]->SetName(mName+"InputCorEvtrkAveSin"+StrFlow);
			if(hInputCorEvtrkAveCos[jFlow]) if(hInputCorEvtrkAveCos[jFlow]->GetEntries()==0) CurrentStage = 2;
		} // jFlow
		
		for(int jFlow=0; jFlow<nFlow; jFlow++) {
			for(int jEvtrkType=0; jEvtrkType<nEvtrkType; jEvtrkType++) {
				if(hInputRawEvtrkAveCos[jFlow]) RawEvtrkSumWts[jFlow][jEvtrkType] = hInputRawEvtrkAveCos[jFlow]->GetBinEntries(jEvtrkType+1);
				if(hInputRawEvtrkAveCos[jFlow]) RawEvtrkAveCos[jFlow][jEvtrkType] = hInputRawEvtrkAveCos[jFlow]->GetBinContent(jEvtrkType+1);
				if(hInputRawEvtrkAveSin[jFlow]) RawEvtrkAveSin[jFlow][jEvtrkType] = hInputRawEvtrkAveSin[jFlow]->GetBinContent(jEvtrkType+1);
				if(hInputCorEvtrkAveCos[jFlow]) CorEvtrkSumWts[jFlow][jEvtrkType] = hInputCorEvtrkAveCos[jFlow]->GetBinEntries(jEvtrkType+1);
				if(hInputCorEvtrkAveCos[jFlow]) CorEvtrkAveCos[jFlow][jEvtrkType] = hInputCorEvtrkAveCos[jFlow]->GetBinContent(jEvtrkType+1);
				if(hInputCorEvtrkAveSin[jFlow]) CorEvtrkAveSin[jFlow][jEvtrkType] = hInputCorEvtrkAveSin[jFlow]->GetBinContent(jEvtrkType+1);
				// XXX: Index: EvtrkType = EventType:TrackType = Cent:Pt:Ch
				int jEventType = jEvtrkType / nTrackType;
				int jCent = jEventType % nCent;
				int jTrackType = jEvtrkType % nTrackType;
				//int jPt = jTrackType / nCh;
				int jCh = jTrackType % nCh;
				RawEvchgSumWts[jFlow][jEventType][jCh] += RawEvtrkSumWts[jFlow][jEvtrkType];
				RawEvchgAveCos[jFlow][jEventType][jCh] += RawEvtrkSumWts[jFlow][jEvtrkType] * RawEvtrkAveCos[jFlow][jEvtrkType];
				RawEvchgAveSin[jFlow][jEventType][jCh] += RawEvtrkSumWts[jFlow][jEvtrkType] * RawEvtrkAveSin[jFlow][jEvtrkType];
				CorEvchgSumWts[jFlow][jEventType][jCh] += CorEvtrkSumWts[jFlow][jEvtrkType];
				CorEvchgAveCos[jFlow][jEventType][jCh] += CorEvtrkSumWts[jFlow][jEvtrkType] * CorEvtrkAveCos[jFlow][jEvtrkType];
				CorEvchgAveSin[jFlow][jEventType][jCh] += CorEvtrkSumWts[jFlow][jEvtrkType] * CorEvtrkAveSin[jFlow][jEvtrkType];
				RawCechgSumWts[jFlow][jCent][jCh] += RawEvtrkSumWts[jFlow][jEvtrkType];
				RawCechgAveCos[jFlow][jCent][jCh] += RawEvtrkSumWts[jFlow][jEvtrkType] * RawEvtrkAveCos[jFlow][jEvtrkType];
				RawCechgAveSin[jFlow][jCent][jCh] += RawEvtrkSumWts[jFlow][jEvtrkType] * RawEvtrkAveSin[jFlow][jEvtrkType];
				CorCechgSumWts[jFlow][jCent][jCh] += CorEvtrkSumWts[jFlow][jEvtrkType];
				CorCechgAveCos[jFlow][jCent][jCh] += CorEvtrkSumWts[jFlow][jEvtrkType] * CorEvtrkAveCos[jFlow][jEvtrkType];
				CorCechgAveSin[jFlow][jCent][jCh] += CorEvtrkSumWts[jFlow][jEvtrkType] * CorEvtrkAveSin[jFlow][jEvtrkType];
			} // jEvtrkType
			for(int jEventType=0; jEventType<nEventType; jEventType++) {
				for(int jCh=0; jCh<nCh; jCh++) {
					if(RawEvchgSumWts[jFlow][jEventType][jCh]>0) {
						RawEvchgAveCos[jFlow][jEventType][jCh] /= RawEvchgSumWts[jFlow][jEventType][jCh];
						RawEvchgAveSin[jFlow][jEventType][jCh] /= RawEvchgSumWts[jFlow][jEventType][jCh];
					} else {
						RawEvchgAveCos[jFlow][jEventType][jCh] = 0;
						RawEvchgAveSin[jFlow][jEventType][jCh] = 0;
					}
					if(CorEvchgSumWts[jFlow][jEventType][jCh]>0) {
						CorEvchgAveCos[jFlow][jEventType][jCh] /= CorEvchgSumWts[jFlow][jEventType][jCh];
						CorEvchgAveSin[jFlow][jEventType][jCh] /= CorEvchgSumWts[jFlow][jEventType][jCh];
					} else {
						CorEvchgAveCos[jFlow][jEventType][jCh] = 0;
						CorEvchgAveSin[jFlow][jEventType][jCh] = 0;
					}
				} // jCh
			} // jEventType
			for(int jCent=0; jCent<nCent; jCent++) {
				for(int jCh=0; jCh<nCh; jCh++) {
					if(RawCechgSumWts[jFlow][jCent][jCh]>0) {
						RawCechgAveCos[jFlow][jCent][jCh] /= RawCechgSumWts[jFlow][jCent][jCh];
						RawCechgAveSin[jFlow][jCent][jCh] /= RawCechgSumWts[jFlow][jCent][jCh];
					} else {
						RawCechgAveCos[jFlow][jCent][jCh] = 0;
						RawCechgAveSin[jFlow][jCent][jCh] = 0;
					}
					if(CorCechgSumWts[jFlow][jCent][jCh]>0) {
						CorCechgAveCos[jFlow][jCent][jCh] /= CorCechgSumWts[jFlow][jCent][jCh];
						CorCechgAveSin[jFlow][jCent][jCh] /= CorCechgSumWts[jFlow][jCent][jCh];
					} else {
						CorCechgAveCos[jFlow][jCent][jCh] = 0;
						CorCechgAveSin[jFlow][jCent][jCh] = 0;
					}
				} // jCh
			} // jCent
		} // jFlow

		TProfile *hInputShfAveCos[nCent][nCh];
		TProfile *hInputShfAveSin[nCent][nCh];
		for(int jCent=0; jCent<nCent; jCent++) {
			TString StrCent = Form("Cent%d", jCent);
			for(int jCh=0; jCh<nCh; jCh++) {
				TString StrCh = Form("Ch%d", jCh);
				hInputShfAveCos[jCent][jCh] = (TProfile*)fInputPhiEtaWeight->Get(StrDirName+"ShfAveCos"+StrCent+StrCh);
				hInputShfAveSin[jCent][jCh] = (TProfile*)fInputPhiEtaWeight->Get(StrDirName+"ShfAveSin"+StrCent+StrCh);
				if(hInputShfAveCos[jCent][jCh]) hInputShfAveCos[jCent][jCh]->SetName(mName+"InputShfAveCos"+StrCent+StrCh);
				if(hInputShfAveSin[jCent][jCh]) hInputShfAveSin[jCent][jCh]->SetName(mName+"InputShfAveSin"+StrCent+StrCh);
				if(hInputShfAveCos[jCent][jCh]) if(hInputShfAveCos[jCent][jCh]->GetEntries()==0 && jCent>0) CurrentStage = 2;
			} // jCh
		} // jCent

		for(int jCent=0; jCent<nCent; jCent++) {
			for(int jCh=0; jCh<nCh; jCh++) {
				for(int jShf=1; jShf<nShf; jShf++) {
					if(hInputShfAveCos[jCent][jCh]) ShfAveCos[jCent][jCh][jShf] = hInputShfAveCos[jCent][jCh]->GetBinContent(jShf);
					if(hInputShfAveSin[jCent][jCh]) ShfAveSin[jCent][jCh][jShf] = hInputShfAveSin[jCent][jCh]->GetBinContent(jShf);
				} // jShf
			} // jCh
		} // jCent

		cout << "MyPhiCorrection: AveCos/Sin loaded." << endl;
		if(CurrentStage == 2) cout << "MyPhiCorrection: corrected AveCos/Sin are not ready. CurrentStage = " << CurrentStage << endl;
		if(CurrentStage == 3) cout << "MyPhiCorrection: AveCos/Sin read out." << endl;
	}

	fInputPhiEtaWeight->Close("R");

	return true;
}


//-----------------------------------------------------------------------------
bool MyPhiCorrection::LoadVzDEta(TString InputRootFile) {
	mPathVzDEta = InputRootFile;

	TFile *fInputVzDEta = new TFile(InputRootFile, "READ");
	if(fInputVzDEta->IsZombie()) {
		cout << "Invalid fInputVzDEta!!! " << flush;
		cout << "No Vz alignment this time." << endl;
		return false;
	}

	// VzDEta only depends on the collision type. (NO Reg, EvenType, TrackType, ...)
	TGraph *gParaVzDEta = (TGraph*)fInputVzDEta->Get(mName+"ParaVzDEta");
	if(!gParaVzDEta) gParaVzDEta = (TGraph*)fInputVzDEta->Get(mRefName+"ParaVzDEta");
	if(gParaVzDEta) {
		funcVzDEta = TF1(mName+"FuncVzDEta", gParaVzDEta->GetTitle());
		funcVzDEta.SetParameters(gParaVzDEta->GetY()); 
	} else {
		cout << "Invalid gParaVzDEta!!! " << flush;
		cout << "No Vz alignment this time." << endl;
		return false;
	}

	fInputVzDEta->Close("R");

	return true;
}

//-----------------------------------------------------------------------------
bool MyPhiCorrection::InputFile(TString fphietaweight, TString fvzdeta) {
	return (LoadPhiEtaWeight(fphietaweight) && LoadVzDEta(fvzdeta));
}


//-----------------------------------------------------------------------------
void MyPhiCorrection::InitEvent(int Cent, double Vz, double EventWeight) {
	mCent = Cent;
	iCent = Cent;

	mVz = Vz;

	mEventWeight = EventWeight;

	// XXX: Index: EvtrkType = EventType:TrackType = Cent:Pt:Ch
	iEventType = iCent;

	// XXX: Vz align Eta distribution
	mDEta = funcVzDEta.Eval(mVz);
}


//-----------------------------------------------------------------------------
void MyPhiCorrection::InitTrack(int Ch, double Pt, double Phi, double Eta, double Wts=1.0) {
	mCh = Ch;
	iCh = Ch<0?0:1;

	mPt = Pt;
	for(int jPt=0; jPt<nPt; jPt++) {
		if(Pt>=mPtEdge[jPt]) iPt = jPt;
		else break;
	}

	mPhi = Phi;
	while(mPhi<0) mPhi += 2*M_PI;
	while(mPhi>=2*M_PI) mPhi -= 2*M_PI;
	iPhi = (int)((mPhi-mPhiL)/(mPhiH-mPhiL)*nPhi);
	if(iPhi<0) iPhi = 0;
	if(iPhi>=nPhi) iPhi = nPhi - 1;

	mEta = Eta + mDEta;
	iEta = (int)((mEta-mEtaL)/(mEtaH-mEtaL)*nEta);
	if(iEta<0) iEta = 0;
	if(iEta>=nEta) iEta = nEta - 1;

	// XXX: Index: EvtrkType = EventType:TrackType = Cent:Pt:Ch
	iTrackType = iPt*nCh + iCh;

	iEvtrkType = iEventType*nTrackType + iTrackType;

	// Shift
	mShfPhi = mPhi;
	if(IsAveKept) {
		for(int jShf=1; jShf<nShf; jShf++) mShfPhi += 2.0/jShf*(-ShfAveSin[iCent][iCh][jShf]*cos(jShf*mPhi) + ShfAveCos[iCent][iCh][jShf]*sin(jShf*mPhi));
		while(mShfPhi<0) mShfPhi += 2*M_PI;
		while(mShfPhi>=2*M_PI) mShfPhi -= 2*M_PI;
	}

	// Input track weight; e.g., 1/efficiency
	mTrackWeight = Wts;
}


//-----------------------------------------------------------------------------
void MyPhiCorrection::Fill() {
	double TmpWeight = PhiEtaWeight[iEventType][iTrackType][iPhi][iEta];
	hRawPhiEtaDist[iEventType][iTrackType]->Fill(mPhi, mEta, mEventWeight*mTrackWeight);
	if(IsAveKept) {
		for(int jFlow=0; jFlow<nFlow; jFlow++) {
			double TmpCos = cos((jFlow+1.0)*mPhi);
			double TmpSin = sin((jFlow+1.0)*mPhi);
			hRawEvtrkAveCos[jFlow]->Fill(iEvtrkType, TmpCos, mEventWeight*mTrackWeight);
			hRawEvtrkAveSin[jFlow]->Fill(iEvtrkType, TmpSin, mEventWeight*mTrackWeight);
		} // jFlow
		if(CurrentStage>=2) {
			//hCorPhiEtaDist[iEventType][iTrackType]->Fill(mPhi, mEta, TmpWeight*mEventWeight*mTrackWeight);
			for(int jFlow=0; jFlow<nFlow; jFlow++) {
				double TmpCos = cos((jFlow+1.0)*mPhi);
				double TmpSin = sin((jFlow+1.0)*mPhi);
				hCorEvtrkAveCos[jFlow]->Fill(iEvtrkType, TmpCos, TmpWeight*mEventWeight*mTrackWeight);
				hCorEvtrkAveSin[jFlow]->Fill(iEvtrkType, TmpSin, TmpWeight*mEventWeight*mTrackWeight);
			} // jFlow
			for(int jShf=1; jShf<nShf; jShf++) {
				double TmpCos = cos(jShf*mPhi);
				double TmpSin = sin(jShf*mPhi);
				hShfAveCos[iCent][iCh]->Fill(jShf, TmpCos, TmpWeight*mEventWeight*mTrackWeight);
				hShfAveSin[iCent][iCh]->Fill(jShf, TmpSin, TmpWeight*mEventWeight*mTrackWeight);
			}
		}
	}
}


//-----------------------------------------------------------------------------
double MyPhiCorrection::GetPhiEtaWeight() const {
	return (double)PhiEtaWeight[iEventType][iTrackType][iPhi][iEta];
}


//-----------------------------------------------------------------------------
void MyPhiCorrection::Finish(TFile *fOut) {
	fOut->cd();
	if(mDir!="") {
		TDirectory *dOut = fOut->mkdir(mDir);
		dOut->cd();
	}

	for(int jEventType=0; jEventType<nEventType; jEventType++) {
		for(int jTrackType=0; jTrackType<nTrackType; jTrackType++) {
			hRawPhiEtaDist[jEventType][jTrackType]->Write();
			//hCorPhiEtaDist[jEventType][jTrackType]->Write();
		} // jTrackType
	} // jEventType
	if(IsAveKept) {
		for(int jFlow=0; jFlow<nFlow; jFlow++) {
			hRawEvtrkAveCos[jFlow]->Write();
			hRawEvtrkAveSin[jFlow]->Write();
			hCorEvtrkAveCos[jFlow]->Write();
			hCorEvtrkAveSin[jFlow]->Write();
		} // jFlow
		for(int jCent=0; jCent<nCent; jCent++) {
			for(int jCh=0; jCh<nCh; jCh++) {
				hShfAveCos[jCent][jCh]->Write();
				hShfAveSin[jCent][jCh]->Write();
			}
		} // jCent
	}

	fOut->cd();
}


//-----------------------------------------------------------------------------
// Get average cos/sin
//-----------------------------------------------------------------------------
bool MyPhiCorrection::IsGoodOrder(int order) const {
	if(order<1 || order>nFlow) {
		cout << "Invalid order in Get...AveCos/Sin!!!" << endl;
		return false;
	}
	return true;
}	

// ------- Evtrk ------- //
double MyPhiCorrection::GetRawEvtrkSumWts(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawEvtrkSumWts[iFlow][iEvtrkType];
}

double MyPhiCorrection::GetRawEvtrkAveCos(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawEvtrkAveCos[iFlow][iEvtrkType];
}

double MyPhiCorrection::GetRawEvtrkAveSin(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawEvtrkAveSin[iFlow][iEvtrkType];
}

double MyPhiCorrection::GetCorEvtrkSumWts(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorEvtrkSumWts[iFlow][iEvtrkType];
}

double MyPhiCorrection::GetCorEvtrkAveCos(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorEvtrkAveCos[iFlow][iEvtrkType];
}

double MyPhiCorrection::GetCorEvtrkAveSin(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorEvtrkAveSin[iFlow][iEvtrkType];
}

// ------- Evchg ------- //
double MyPhiCorrection::GetRawEvchgSumWts(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawEvchgSumWts[iFlow][iEventType][iCh];
}

double MyPhiCorrection::GetRawEvchgAveCos(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawEvchgAveCos[iFlow][iEventType][iCh];
}

double MyPhiCorrection::GetRawEvchgAveSin(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawEvchgAveSin[iFlow][iEventType][iCh];
}

double MyPhiCorrection::GetCorEvchgSumWts(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorEvchgSumWts[iFlow][iEventType][iCh];
}

double MyPhiCorrection::GetCorEvchgAveCos(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorEvchgAveCos[iFlow][iEventType][iCh];
}

double MyPhiCorrection::GetCorEvchgAveSin(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorEvchgAveSin[iFlow][iEventType][iCh];
}

// ------- Cechg ------- //
double MyPhiCorrection::GetRawCechgSumWts(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawCechgSumWts[iFlow][iCent][iCh];
}

double MyPhiCorrection::GetRawCechgAveCos(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawCechgAveCos[iFlow][iCent][iCh];
}

double MyPhiCorrection::GetRawCechgAveSin(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)RawCechgAveSin[iFlow][iCent][iCh];
}

double MyPhiCorrection::GetCorCechgSumWts(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorCechgSumWts[iFlow][iCent][iCh];
}

double MyPhiCorrection::GetCorCechgAveCos(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorCechgAveCos[iFlow][iCent][iCh];
}

double MyPhiCorrection::GetCorCechgAveSin(int order) const {
	if(!IsGoodOrder(order)) return 0;
	int iFlow = order - 1;
	return (double)CorCechgAveSin[iFlow][iCent][iCh];
}
