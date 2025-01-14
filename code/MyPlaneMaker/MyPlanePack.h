#ifndef MYPLANEPACK_H
#define MYPLANEPACK_H

#include <TVector2.h>
#include <TString.h>
#include <TProfile2D.h>
#include "MyPlaneInfo.h"
#include <iostream>
#include <TDirectoryFile.h>


class MyPlanePack
{
private:
	TString Name;

	bool IsAutoFill;

	double Weight;

	bool IsGood;
	double RawReso[_NMAX+1];
	double CorReso[_NMAX+1];

	TProfile2D *hCosFull;
	TProfile2D *hSinFull;
	TProfile2D *hCosEast;
	TProfile2D *hSinEast;
	TProfile2D *hCosWest;
	TProfile2D *hSinWest;

	TH1D *hRawFullDist[_NMAX+1];
	TH1D *hRawEastDist[_NMAX+1];
	TH1D *hRawWestDist[_NMAX+1];
	TH2D *hRawSubsDist[_NMAX+1];
	TH1D *hCorFullDist[_NMAX+1];
	TH1D *hCorEastDist[_NMAX+1];
	TH1D *hCorWestDist[_NMAX+1];
	TH2D *hCorSubsDist[_NMAX+1];

	bool IsGoodInput;

	void ClaimHistogram();

public:
	MyPlanePack();
	MyPlanePack(TString name);
	~MyPlanePack();

	MyPlaneInfo Full;
	MyPlaneInfo East; // sub-event eta<0
	MyPlaneInfo West; // sub-event eta>0

	void ResetEvent();
	bool InputHistogram(TString inputpath);
	bool InputHistogram(TDirectoryFile *fInput);

	void Add(double phi, double eta, double wts);
	void Calc();
	void Fill();
	void Finish();

	void SetName(TString name);
	TString GetName() const;
	void SetIsAutoFill(bool isautofill);
	bool GetIsAutoFill() const;

	void SetWeight(double weight);
	double GetWeight() const;

	bool GetIsGood() const;
	double GetRawReso(int i) const;
	double GetCorReso(int i) const;
};


inline void MyPlanePack::SetName(TString name) { Name = name; }
inline TString MyPlanePack::GetName() const { return Name; }
inline void MyPlanePack::SetIsAutoFill(bool isautofill) { IsAutoFill = isautofill; }
inline bool MyPlanePack::GetIsAutoFill() const { return IsAutoFill; }

inline void MyPlanePack::SetWeight(double weight) { Weight = weight; }
inline double MyPlanePack::GetWeight() const { return Weight; }

inline bool MyPlanePack::GetIsGood() const { return IsGood; }
inline double MyPlanePack::GetRawReso(int i) const {
	if(i>=1 && i<=_NMAX) {
		return RawReso[i];
	} else {
		std::cout << "ERROR: invalid Q-order!" << std::endl;
		return -999;
	}
}
inline double MyPlanePack::GetCorReso(int i) const {
	if(i>=1 && i<=_NMAX) {
		return CorReso[i];
	} else {
		std::cout << "ERROR: invalid Q-order!" << std::endl;
		return -999;
	}
}


#endif
