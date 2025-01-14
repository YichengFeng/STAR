/***************************************************************************
 * Author: Yicheng Feng
 * Email: feng216@purdue.edu, fengyich@outlook.com
 * Date: 2023/09/17
 * Note: To reweight phi so that the phi distribution is uniform.
 *       The reweighting is done for various EventType and TrackType.
 *       depending on their VertexZ, Centrality; Charge, Pt, Eta.
 *       Eta distribution is aligned to VertexZ=0. 
 ***************************************************************************/


#ifndef MyPhiCorrection_H
#define MyPhiCorrection_H

#include <iostream>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TF1.h>


class MyPhiCorrection
{
private:
	// constants
	static const int nFlow = 3;
	static const int nCent = 10;
	static const int nEventType = 10; // (nCent)
	static const int nCh = 2;
	static const int nPt = 3;
	static const int nTrackType = 6; // (nPt:nCh)
	static const int nPhi = 120;
	static const int nEta = 36; // TODO: iTPC
	static const int nEvtrkType = 60; // (nEventType:nTrackType)

	static const int nShf = 13;

	// input
	TString mDir;
	TString mName;
	double mEventWeight;
	std::vector<double> mPtEdge; // size = nPt
	double mPhiL, mPhiH;
	double mEtaL, mEtaH;

	TString mRefDir;
	TString mRefName;

	TString mPathPhiEtaWeight;
	TString mPathVzDEta;

	// index
	// XXX: Index: EvtrkType = EventType:TrackType = Cent:Pt:Ch
	int iCent;
	int iEventType;
	int iCh;
	int iPt;
	int iTrackType;
	int iPhi;
	int iEta;
	int sEta;
	int iEvtrkType;

	// values
	int    mCent;
	double mVz;
	int    mCh;
	double mPt;
	double mPhi;
	double mEta;
	double mShfPhi;
	double mTrackWeight;

	// average
	bool IsAveKept;

	// weight cutoff
	bool IsWeightCutOff;
	double mWtL;
	double mWtH;

	// step
	int CurrentStage; // 1 for step 1; 2 for step 2; 3 for step 3;

	// Vz align
	TF1 funcVzDEta;
	double mDEta;

	// weight
	float PhiEtaWeight[nEventType][nTrackType][nPhi][nEta];
	// average cos/sin
	float RawEvtrkSumWts[nFlow][nEvtrkType];
	float RawEvtrkAveCos[nFlow][nEvtrkType];
	float RawEvtrkAveSin[nFlow][nEvtrkType];
	float CorEvtrkSumWts[nFlow][nEvtrkType];
	float CorEvtrkAveCos[nFlow][nEvtrkType];
	float CorEvtrkAveSin[nFlow][nEvtrkType];
	float RawEvchgSumWts[nFlow][nEventType][nCh];
	float RawEvchgAveCos[nFlow][nEventType][nCh];
	float RawEvchgAveSin[nFlow][nEventType][nCh];
	float CorEvchgSumWts[nFlow][nEventType][nCh];
	float CorEvchgAveCos[nFlow][nEventType][nCh];
	float CorEvchgAveSin[nFlow][nEventType][nCh];
	float RawCechgSumWts[nFlow][nCent][nCh];
	float RawCechgAveCos[nFlow][nCent][nCh];
	float RawCechgAveSin[nFlow][nCent][nCh];
	float CorCechgSumWts[nFlow][nCent][nCh];
	float CorCechgAveCos[nFlow][nCent][nCh];
	float CorCechgAveSin[nFlow][nCent][nCh];
	// shift
	float ShfAveCos[nCent][nCh][nShf];
	float ShfAveSin[nCent][nCh][nShf];

	// histograms
	TH2F *hRawPhiEtaDist[nEventType][nTrackType];
	TH2F *hCorPhiEtaDist[nEventType][nTrackType];
	// profiles
	TProfile *hRawEvtrkAveCos[nFlow];
	TProfile *hRawEvtrkAveSin[nFlow];
	TProfile *hCorEvtrkAveCos[nFlow];
	TProfile *hCorEvtrkAveSin[nFlow];
	TProfile *hShfAveCos[nCent][nCh];
	TProfile *hShfAveSin[nCent][nCh];

public:
	void Init(TString Name);
	MyPhiCorrection();
	MyPhiCorrection(TString Name); // Please only use 0~9, a~Z, _
	MyPhiCorrection(TString Dir, TString Name); // Please only use 0~9, a~Z, _
	~MyPhiCorrection();

	// input
	//  The constructors should be enough!
	//  If no necessary, don't use the following methods.
	void SetName(TString Name) { mName = Name; }
	void SetDir(TString Dir) { mDir = Dir; }
	void SetRefName(TString RefName) { mRefName = RefName; }
	void SetRefDir(TString RefDir) { mRefDir = RefDir; }
	void SetEventWeight(double EventWeight) { mEventWeight = EventWeight; }
	void SetPtEdge(std::vector<double> PtEdge) {
		if(PtEdge.size() == nPt) {
			mPtEdge = PtEdge;
		} else {
			std::cout << "PtEdge size does NOT match!!! " << std::flush;
			std::cout << "The default PtEdge is used." << std::endl;
		}
	}
	void SetPhiRange(double PhiL, double PhiH) { mPhiL = PhiL; mPhiH = PhiH; }
	void SetPhiRange(std::vector<double> tmp) { mPhiL = tmp[0]; mPhiH = tmp[1]; }
	void SetEtaRange(double EtaL, double EtaH) { mEtaL = EtaL; mEtaH = EtaH; }
	void SetEtaRange(std::vector<double> tmp) { mEtaL = tmp[0]; mEtaH = tmp[1]; }
	void SetAveKept(bool isavekept) { IsAveKept = isavekept; }
	void SetWeightCutOff(bool isweightcutoff, double mwtl=0.10, double mwth=10.0) { IsWeightCutOff = isweightcutoff; mWtL = mwtl; mWtH = mwth; }
	void SetWeightCutOff(std::vector<double> tmp) { IsWeightCutOff = (tmp[0]>0); mWtL = tmp[1]; mWtH = tmp[2]; }

	// output
	TString GetName() const { return mName; }
	TString GetDir() const { return mDir; }
	TString GetRefName() const { return mRefName; }
	TString GetRefDir() const { return mRefDir; }
	double GetEventWeight() const { return mEventWeight; }
	std::vector<double> GetPtEdge() const { return mPtEdge; }
	std::vector<double> GetPhiRange() const { std::vector<double> PhiRange{mPhiL, mPhiH}; return PhiRange; }
	std::vector<double> GetEtaRange() const { std::vector<double> EtaRange{mEtaL, mEtaH}; return EtaRange; }
	bool GetAveKept() const { return IsAveKept; }
	std::vector<double> GetWeightCutOff() const { std::vector<double> tmp{IsWeightCutOff?1.0:-1.0, mWtL, mWtH}; return tmp; }
	double GetDEta() const { return mDEta; }

	int GetIdCent     () const { return iCent;      }
	int GetIdEventType() const { return iEventType; }
	int GetIdCh       () const { return iCh;        }
	int GetIdPt       () const { return iPt;        }
	int GetIdPhi      () const { return iPhi;       }
	int GetIdEta      () const { return iEta;       }
	int GetIdTrackType() const { return iTrackType; }

	// work
	bool LoadPhiEtaWeight(TString InputRootFile);
	TString GetPathPhiEtaWeight() const { return mPathPhiEtaWeight; }
	bool LoadVzDEta(TString InputRootFile);
	TString GetPathVzDEta() const { return mPathVzDEta; }
	bool InputFile(TString fphietaweight, TString fvzdeta);
	void InitEvent(int Cent, double Vz, double EventWeight);
	void InitTrack(int Ch, double Pt, double Phi, double Eta, double Wts);
	void Fill();
	double GetPhiEtaWeight() const;
	void Finish(TFile *fOut);

	bool IsGoodOrder(int order) const;
	double GetRawEvtrkSumWts(int order) const;
	double GetRawEvtrkAveCos(int order) const;
	double GetRawEvtrkAveSin(int order) const;
	double GetCorEvtrkSumWts(int order) const;
	double GetCorEvtrkAveCos(int order) const;
	double GetCorEvtrkAveSin(int order) const;
	double GetRawEvchgSumWts(int order) const;
	double GetRawEvchgAveCos(int order) const;
	double GetRawEvchgAveSin(int order) const;
	double GetCorEvchgSumWts(int order) const;
	double GetCorEvchgAveCos(int order) const;
	double GetCorEvchgAveSin(int order) const;
	double GetRawCechgSumWts(int order) const;
	double GetRawCechgAveCos(int order) const;
	double GetRawCechgAveSin(int order) const;
	double GetCorCechgSumWts(int order) const;
	double GetCorCechgAveCos(int order) const;
	double GetCorCechgAveSin(int order) const;

	double GetShfPhi() const { return mShfPhi; }
	double GetEta() const { return mEta; }

	static const int GetnFlow() { return nFlow; }
	static const int GetnCent() { return nCent; }
	static const int GetnEventType() { return nEventType; }
	static const int GetnCh() { return nCh; }
	static const int GetnPt() { return nPt; }
	static const int GetnTrackType() { return nTrackType; }
	static const int GetnPhi() { return nPhi; }
	static const int GetnEta() { return nEta; }
};

#endif
