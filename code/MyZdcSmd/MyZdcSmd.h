/***************************************************************************
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com, feng216@purdue.edu
 * Date: 2023/03/20
 * Note: this class reads in the ZDCSMD information
 *       calibrates the ZDCSMD (gain rate, recenter, shift correction)
 *       then gets the first-order Q-vector and plane Psi
 *       it also makes QA plots
 ***************************************************************************/

#ifndef MyZdcSmd_H
#define MyZdcSmd_H

#include <iostream>
#include <TString.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
//#include "StPicoDstMaker/StPicoDstMaker.h"
//#include "StPicoEvent/StPicoDst.h"
//#include "StPicoEvent/StPicoTrack.h"
#ifndef MyZdcSmdPlot
	#include "StPicoEvent/StPicoEvent.h"
#endif
//#include "StPicoEvent/StPicoBTofPidTraits.h"


class MyZdcSmd
{
private:
	// constant
	static const int nCent = 10;
	static const int nEW = 2; // [East,West]
	static const int nVH = 2; // [Vert,Hori]=[X,Y]
	static const int nST = 8; // [Slats]
	static const int nAL = 32; // nEW*nVH*nST;
	const TString StrEW[nEW] = {"East", "West"};
	const TString StrVH[nVH] = {"Vert", "Hori"};
	const TString StrXY[nVH] = {"X", "Y"};
	static const int nShf = 13;

	double Position[nEW][nVH][nST];
	double Poserror[nEW][nVH][nST];

	// current step
	int mStep; // 1: PedeGain+Recenter; 2: Shfit; 3: Ready;

	// file I/O
	TH1D *hZdcSmdRawAdc[nEW][nVH][nST];
	TH1D *hZdcSmdCorAdc[nEW][nVH][nST];
	TProfile *hZdcSmdPede;
	TProfile *hZdcSmdGain;
	TProfile *hZdcSmdPeak;
	TProfile *hZdcRawCenter;
	TProfile *hZdcCorCenter;
	TProfile *hRawCenter[nEW][nVH]; // [Vert,Hori]=[X,Y]
	TProfile *hCorCenter[nEW][nVH]; // [Vert,Hori]=[X,Y]
	double RawCenter[nEW][nVH]; // [Vert,Hori]=[X,Y]
	double CorCenter[nEW][nVH]; // [Vert,Hori]=[X,Y]
	double RefCenter[nEW][nVH]; // [Vert,Hori]=[X,Y]
	TProfile *hEastRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hEastRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hWestRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hWestRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hFullRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hFullRctPsiSin[nCent]; // Psi=Psi1
	double EastRctPsiCos[nCent][nShf];
	double EastRctPsiSin[nCent][nShf];
	double WestRctPsiCos[nCent][nShf];
	double WestRctPsiSin[nCent][nShf];
	double FullRctPsiCos[nCent][nShf];
	double FullRctPsiSin[nCent][nShf];

	TProfile *hEastPvRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hEastPvRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hWestPvRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hWestPvRctPsiSin[nCent]; // Psi=Psi1
	TProfile *hFullPvRctPsiCos[nCent]; // Psi=Psi1
	TProfile *hFullPvRctPsiSin[nCent]; // Psi=Psi1
	double EastPvRctPsiCos[nCent][nShf];
	double EastPvRctPsiSin[nCent][nShf];
	double WestPvRctPsiCos[nCent][nShf];
	double WestPvRctPsiSin[nCent][nShf];
	double FullPvRctPsiCos[nCent][nShf];
	double FullPvRctPsiSin[nCent][nShf];

	TProfile *hAvePvx;
	TProfile *hAvePvy;
	double AvePvx;
	double AvePvy;

	// input
	TString mDir;
	TString mName;
	TString mInputFilePath;
	int mCent;
	double mWght;
	TVector3 mPv; // TPC primary vertex
	double mPvx;
	double mPvy;

	TString mRefDir;
	TString mRefName;

	int mVerbose;
	int mMerge;
	static const int nGainMod = 4;
	int mGainMod; 
	// 0: parameter in input rootfile, or from fit if input does not have parameter; (default, invalid modes go back to this)
	// 1: fixed to 1; 
	// 2: correction from SetGain(); 
	// 3: must come from fit
	static const int nCenterMod = 3;
	int mCenterMod; 
	// 0: read > calc > ref; (default, invalid modes go back to this)
	// 1: fixed to reference
	// 2: gaussian center

	bool mIsAutoFill; // default: true; if false, must call Fill() later to fill this event

	bool mIsGoodInput;
	bool mIsGoodEvent;

	// data
	double mZdcSmdRawAdc[nEW][nVH][nST]; // [East,West][Vert,Hori][Slats]
	double mZdcSmdCorAdc[nEW][nVH][nST]; // [East,West][Vert,Hori][Slats]
	double mZdcSmdPede[nEW][nVH][nST]; // [East,West][Vert,Hori][Slats]
	double mZdcSmdGain[nEW][nVH][nST]; // [East,West][Vert,Hori][Slats]
	double mZdcSmdPeak[nEW][nVH][nST]; // [East,West][Vert,Hori][Slats]

	// output
	// Q=Q1, Psi=Psi1
	TVector2 mEastRawQ; // Raw
	TVector2 mEastCorQ; // Raw + PedeGain
	TVector2 mEastRctQ; // Raw + PedeGain + Recenter
	TVector2 mEastShfQ; // Raw + PedeGain + Recenter + Shift
	TVector2 mWestRawQ;
	TVector2 mWestCorQ;
	TVector2 mWestRctQ;
	TVector2 mWestShfQ;
	TVector2 mFullRawQ;
	TVector2 mFullCorQ;
	TVector2 mFullRctQ;
	TVector2 mFullShfQ;
	double mEastRawPsi; // Raw
	double mEastCorPsi; // Raw + PedeGain
	double mEastRctPsi; // Raw + PedeGain + Recenter
	double mEastShfPsi; // Raw + PedeGain + Recenter + Shift
	double mWestRawPsi;
	double mWestCorPsi;
	double mWestRctPsi;
	double mWestShfPsi;
	double mFullRawPsi;
	double mFullCorPsi;
	double mFullRctPsi;
	double mFullShfPsi;
	TH1D *hEastRawPsi[nCent];
	TH1D *hEastCorPsi[nCent];
	TH1D *hEastRctPsi[nCent];
	TH1D *hEastShfPsi[nCent];
	TH1D *hWestRawPsi[nCent];
	TH1D *hWestCorPsi[nCent];
	TH1D *hWestRctPsi[nCent];
	TH1D *hWestShfPsi[nCent];
	TH1D *hFullRawPsi[nCent];
	TH1D *hFullCorPsi[nCent];
	TH1D *hFullRctPsi[nCent];
	TH1D *hFullShfPsi[nCent];
	TH2D *h2SubRawPsi[nCent];
	TH2D *h2SubCorPsi[nCent];
	TH2D *h2SubRctPsi[nCent];
	TH2D *h2SubShfPsi[nCent];
	TProfile *hResoRawPsi;
	TProfile *hResoCorPsi;
	TProfile *hResoRctPsi;
	TProfile *hResoShfPsi;
	TH2D *h2EastRawQ;
	TH2D *h2EastCorQ;
	TH2D *h2WestRawQ;
	TH2D *h2WestCorQ;
	TH2D *h2FullRawQ;
	TH2D *h2FullCorQ;

	TVector2 mEastPvRctQ; // Raw + PedeGain + PvRecenter
	TVector2 mEastPvShfQ; // Raw + PedeGain + PvRecenter + Shift
	TVector2 mWestPvRctQ;
	TVector2 mWestPvShfQ;
	TVector2 mFullPvRctQ;
	TVector2 mFullPvShfQ;
	double mEastPvRctPsi; // Raw + PedeGain + PvRecenter
	double mEastPvShfPsi; // Raw + PedeGain + PvRecenter + Shift
	double mWestPvRctPsi;
	double mWestPvShfPsi;
	double mFullPvRctPsi;
	double mFullPvShfPsi;
	TH1D *hEastPvRctPsi[nCent];
	TH1D *hEastPvShfPsi[nCent];
	TH1D *hWestPvRctPsi[nCent];
	TH1D *hWestPvShfPsi[nCent];
	TH1D *hFullPvRctPsi[nCent];
	TH1D *hFullPvShfPsi[nCent];
	TH2D *h2SubPvRctPsi[nCent];
	TH2D *h2SubPvShfPsi[nCent];
	TProfile *hResoPvRctPsi;
	TProfile *hResoPvShfPsi;

public:
	void ResetEvent();
	void Init();
	MyZdcSmd();
	MyZdcSmd(TString name);
	MyZdcSmd(TString dir, TString name);
	~MyZdcSmd();

	bool GetIsAutoFill() const { return mIsAutoFill; }
	void SetIsAutoFill(bool isautofill) { mIsAutoFill = isautofill; }

	static bool IsGoodSlat(int iEW, int iHV, int iST);
	static double GetPosition(int iEW, int iVH, int iST);
	double GetRctPosition(int iEW, int iVH, int iST);
	double GetPvRctPosition(int iEW, int iVH, int iST);
	double ZdcSmdAdcRawToCor(int iEW, int iVH, int iST, double raw);
	static int GetLastBin(TH1D *h);
	bool InputFile(TString filepath);
	TString GetInputFilePath() const { return mInputFilePath; }
	static double Phi02Pi(double phi);
	void CalcAve();
	void CalcGau();
	void CalcPsi();
	void Calc();
	void Fill();
	void Finish(TFile *fOut);

	void InitEvent(int cent, double *zdcsmdadc,  double wght);
	void InitEvent(int cent, float  *zdcsmdadc,  double wght);
#ifndef MyZdcSmdPlot
	void InitEvent(int cent, StPicoEvent *event, double wght);
#endif
	void InitEvent(int cent, double *zdcsmdadc,  TVector3 pv, double wght);
	void InitEvent(int cent, float  *zdcsmdadc,  TVector3 pv, double wght);
#ifndef MyZdcSmdPlot
	void InitEvent(int cent, StPicoEvent *event, TVector3 pv, double wght);
#endif

	// I/O
	void SetName(TString name) { mName = name; }
	void SetDir(TString dir) { mDir = dir; }
	TString GetName() const { return mName; }
	TString GetDir() const { return mDir; }
	void SetRefName(TString name) { mRefName = name; }
	void SetRefDir(TString dir) { mRefDir = dir; }
	TString GetRefName() const { return mRefName; }
	TString GetRefDir() const { return mRefDir; }
	void SetWght(double mwght) { mWght = mwght; } // event weight, default=1.0;
	double GetWght() const { return mWght; }
	void SetCent(int cent) { mCent = cent; }
	int GetCent() const { return mCent; }

	void SetVerbose(int verbose) { mVerbose = verbose; } // set the print level: 0,1,2 minimum->maximum
	void SetMerge(int merge) { mMerge = merge; } // how to merge sub to full: 0: weighted sum Q, 1: directly sum aveQ.
	void SetGainMod(int gainmod) {
		if(gainmod>=0 && gainmod<nGainMod) {
			mGainMod = gainmod;
		} else {
			mGainMod = 0;
			std::cout << "MyZdcSmd: invalid gain mode, reset to default mGainMod = " << mGainMod << std::endl;
		}
	}
	void SetCenterMod(int centermod) {
		if(centermod>=0 && centermod<nCenterMod) {
			mCenterMod = centermod;
		} else {
			mCenterMod = 0;
			std::cout << "MyZdcSmd: invalid center mode, reset to default mCenterMod = " << mCenterMod << std::endl;
		}
	}
	int GetVerbose() const { return mVerbose; }
	int GetMerge() const { return mMerge; }
	int GetGainMod() const { return mGainMod; }
	int GetCenterMod() const { return mCenterMod; }

	std::vector<double> GetPede() const;
	std::vector<double> GetGain() const;
	std::vector<double> GetPeak() const;
	void SetPede(std::vector<double> tmp);
	void SetGain(std::vector<double> tmp);
	void SetPeak(std::vector<double> tmp);

	bool IsGoodEvent() const { return mIsGoodEvent; }

	TVector2 GetEastRawCenter() const { return TVector2(RawCenter[0][0], RawCenter[0][1]); }
	TVector2 GetWestRawCenter() const { return TVector2(RawCenter[1][0], RawCenter[1][1]); }
	TVector2 GetEastCorCenter() const { return TVector2(CorCenter[0][0], CorCenter[0][1]); }
	TVector2 GetWestCorCenter() const { return TVector2(CorCenter[1][0], CorCenter[1][1]); }
	TVector2 GetEastRefCenter() const { return TVector2(RefCenter[0][0], RefCenter[0][1]); }
	TVector2 GetWestRefCenter() const { return TVector2(RefCenter[1][0], RefCenter[1][1]); }
	void SetEastRefCenter(double x, double y) { RefCenter[0][0] = x; RefCenter[0][1] = y; }
	void SetWestRefCenter(double x, double y) { RefCenter[1][0] = x; RefCenter[1][1] = y; }
	void SetEastRefCenter(const TVector2 v) { RefCenter[0][0] = v.X(); RefCenter[0][1] = v.Y(); }
	void SetWestRefCenter(const TVector2 v) { RefCenter[1][0] = v.X(); RefCenter[1][1] = v.Y(); }

	double GetAvePvx() const { return AvePvx; }
	double GetAvePvy() const { return AvePvy; }
	double GetPvx() const { return mPvx; }
	double GetPvy() const { return mPvy; }

	double GetZdcSmdRawAdc(int iEW, int iVH, int iST) const { if(IsGoodSlat(iEW, iVH, iST)) { return mZdcSmdRawAdc[iEW][iVH][iST]; } else { return -999; } }
	double GetZdcSmdCorAdc(int iEW, int iVH, int iST) const { if(IsGoodSlat(iEW, iVH, iST)) { return mZdcSmdCorAdc[iEW][iVH][iST]; } else { return -999; } }
	double GetZdcSmdPede(int iEW, int iVH, int iST) const { if(IsGoodSlat(iEW, iVH, iST)) { return mZdcSmdPede[iEW][iVH][iST]; } else { return -999; } }
	double GetZdcSmdGain(int iEW, int iVH, int iST) const { if(IsGoodSlat(iEW, iVH, iST)) { return mZdcSmdGain[iEW][iVH][iST]; } else { return -999; } }

	TVector2 GetEastRawQ() const { return mEastRawQ; }
	TVector2 GetEastCorQ() const { return mEastCorQ; }
	TVector2 GetEastRctQ() const { return mEastRctQ; }
	TVector2 GetEastShfQ() const { return mEastShfQ; }
	TVector2 GetWestRawQ() const { return mWestRawQ; }
	TVector2 GetWestCorQ() const { return mWestCorQ; }
	TVector2 GetWestRctQ() const { return mWestRctQ; }
	TVector2 GetWestShfQ() const { return mWestShfQ; }
	TVector2 GetFullRawQ() const { return mFullRawQ; }
	TVector2 GetFullCorQ() const { return mFullCorQ; }
	TVector2 GetFullRctQ() const { return mFullRctQ; }
	TVector2 GetFullShfQ() const { return mFullShfQ; }
	double GetEastRawPsi() const { return mEastRawPsi; }
	double GetEastCorPsi() const { return mEastCorPsi; }
	double GetEastRctPsi() const { return mEastRctPsi; }
	double GetEastShfPsi() const { return mEastShfPsi; }
	double GetWestRawPsi() const { return mWestRawPsi; }
	double GetWestCorPsi() const { return mWestCorPsi; }
	double GetWestRctPsi() const { return mWestRctPsi; }
	double GetWestShfPsi() const { return mWestShfPsi; }
	double GetFullRawPsi() const { return mFullRawPsi; }
	double GetFullCorPsi() const { return mFullCorPsi; }
	double GetFullRctPsi() const { return mFullRctPsi; }
	double GetFullShfPsi() const { return mFullShfPsi; }

	TVector2 GetEastPvRctQ() const { return mEastPvRctQ; }
	TVector2 GetEastPvShfQ() const { return mEastPvShfQ; }
	TVector2 GetWestPvRctQ() const { return mWestPvRctQ; }
	TVector2 GetWestPvShfQ() const { return mWestPvShfQ; }
	TVector2 GetFullPvRctQ() const { return mFullPvRctQ; }
	TVector2 GetFullPvShfQ() const { return mFullPvShfQ; }
	double GetEastPvRctPsi() const { return mEastPvRctPsi; }
	double GetEastPvShfPsi() const { return mEastPvShfPsi; }
	double GetWestPvRctPsi() const { return mWestPvRctPsi; }
	double GetWestPvShfPsi() const { return mWestPvShfPsi; }
	double GetFullPvRctPsi() const { return mFullPvRctPsi; }
	double GetFullPvShfPsi() const { return mFullPvShfPsi; }
};

#endif
