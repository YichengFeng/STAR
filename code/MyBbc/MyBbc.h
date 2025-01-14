/***************************************************************************
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com, feng216@purdue.edu
 * Date: 2023/10/06
 * Note: this class reads in the BBC information
 *       calibrates the BBC (gain rate, recenter, shift correction)
 *       then gets the first-order Q-vector and plane Psi
 *       it also makes QA plots
 ***************************************************************************/

#ifndef MyBbc_H
#define MyBbc_H

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
#include "../MyPlaneMaker/MyPlaneInfo.h"
//#include "../MyPlaneMaker/MyPlaneInfo.cxx"
#include "../MyPlaneMaker/MyPlanePack.h"
//#include "../MyPlaneMaker/MyPlanePack.cxx"
//#include "StPicoDstMaker/StPicoDstMaker.h"
//#include "StPicoEvent/StPicoDst.h"
//#include "StPicoEvent/StPicoTrack.h"
#ifndef MyBbcPlot
	#include "StPicoEvent/StPicoEvent.h"
#endif
//#include "StPicoEvent/StPicoBTofPidTraits.h"


class MyBbc
{
private:
	// constant
	static const int nCent = 10;
	static const int nEW = 2; // [East,West]
	static const int nTL = 16; // [Tiles]
	static const int nAL = 32; // nEW*nTL;
	const TString StrEW[nEW] = {"East", "West"};
	static const int nShf = 13;

	// current step
	int mStep; // 1: PedeGain; 2: Shfit; 3: Ready;

	// file I/O
	TH1D *hBbcRawAdc[nEW][nTL];
	TH1D *hBbcCorAdc[nEW][nTL];
	TProfile *hBbcPede;
	TProfile *hBbcGain;
	TProfile *hBbcPeak;

	// input
	TString mDir;
	TString mName;
	TString mInputFilePath;
	int mCent;
	double mWght;

	TString mRefDir;
	TString mRefName;

	int mVerbose;
	static const int nGainMod = 4;
	int mGainMod; 
	// 0: parameter in input rootfile, or from fit if input does not have parameter; (default, invalid modes go back to this)
	// 1: fixed to 1; 
	// 2: correction from SetGain(); 
	// 3: must come from fit

	bool mIsAutoFill; // default: true; if false, must call Fill() later to fill this event

	bool mIsGoodInput;
	bool mIsGoodEvent;

	// data
	double mBbcRawAdc[nEW][nTL]; // [East,West][Tiles]
	double mBbcCorAdc[nEW][nTL]; // [East,West][Tiles]
	double mBbcPede[nEW][nTL]; // [East,West][Tiles]
	double mBbcGain[nEW][nTL]; // [East,West][Tiles]
	double mBbcPeak[nEW][nTL]; // [East,West][Tiles]

	// plane
	MyPlanePack *mRawBbcEp[nCent]; // only gives Raw, no input
	MyPlanePack *mCorBbcEp[nCent]; // gives both Cor and Shf, with input

	TProfile *hResoRawPsi[_NMAX+1];
	TProfile *hResoCorPsi[_NMAX+1];
	TProfile *hResoShfPsi[_NMAX+1];

public:
	void ResetEvent(bool isinit);
	void Init();
	MyBbc();
	MyBbc(TString name);
	MyBbc(TString dir, TString name);
	~MyBbc();

	bool GetIsAutoFill() const { return mIsAutoFill; }
	void SetIsAutoFill(bool isautofill) { mIsAutoFill = isautofill; }

	static bool IsGoodTile(int iEW, int iTL);
	static double GetAngle(int iEW, int iTL);
	double BbcAdcRawToCor(int iEW, int iTL, double raw);
	bool InputFile(TString filepath);
	TString GetInputFilePath() const { return mInputFilePath; }
	void Calc();
	void Fill();
	void Finish(TFile *fOut);

	void InitEvent(int cent, double *bbcadc,  double wght);
	void InitEvent(int cent, float  *bbcadc,  double wght);
#ifndef MyBbcPlot
	void InitEvent(int cent, StPicoEvent *event, double wght);
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
	void SetGainMod(int gainmod) {
		if(gainmod>=0 && gainmod<nGainMod) {
			mGainMod = gainmod;
		} else {
			mGainMod = 0;
			std::cout << "MyBbc: invalid gain mode, reset to default mGainMod = " << mGainMod << std::endl;
		}
	}
	int GetVerbose() const { return mVerbose; }
	int GetGainMod() const { return mGainMod; }

	std::vector<double> GetPede() const;
	std::vector<double> GetGain() const;
	std::vector<double> GetPeak() const;
	void SetPede(std::vector<double> tmp);
	void SetGain(std::vector<double> tmp);
	void SetPeak(std::vector<double> tmp);

	bool IsGoodEvent() const { return mIsGoodEvent; }

	double GetBbcRawAdc(int iEW, int iTL) const { if(IsGoodTile(iEW, iTL)) { return mBbcRawAdc[iEW][iTL]; } else { return -999; } }
	double GetBbcCorAdc(int iEW, int iTL) const { if(IsGoodTile(iEW, iTL)) { return mBbcCorAdc[iEW][iTL]; } else { return -999; } }
	double GetBbcPede(int iEW, int iTL) const { if(IsGoodTile(iEW, iTL)) { return mBbcPede[iEW][iTL]; } else { return -999; } }
	double GetBbcGain(int iEW, int iTL) const { if(IsGoodTile(iEW, iTL)) { return mBbcGain[iEW][iTL]; } else { return -999; } }

	double GetFullRawPsi(int ord);
	double GetFullCorPsi(int ord);
	double GetFullShfPsi(int ord);
	double GetEastRawPsi(int ord);
	double GetEastCorPsi(int ord);
	double GetEastShfPsi(int ord);
	double GetWestRawPsi(int ord);
	double GetWestCorPsi(int ord);
	double GetWestShfPsi(int ord);
};

#endif
