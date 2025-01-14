#ifndef MyEventMixer_H
#define MyEventMixer_H

#include <iostream>
#include <cmath>
#include <vector>
#include "MyEventToMix.h"

using namespace std;


class MyEventMixer
{
private:
	static const int EM_MAXEVENT = 1;
	static const int EM_NCENT = 10;
	static const int EM_NVZ = 60;
	static const int EM_NEVENTTYPE = 600; // (EM_NCENT*EM_NVZ)

	int mCent;
	double mVz;
	double mVzL, mVzH, mVzM;
	double mEventWeight;

	int iCent;
	int iVz;
	int iEventType;

	int CurrentMixSize[EM_NEVENTTYPE]; // no more than EM_MAXEVENT
	MyEventToMix mEventToMix[EM_NEVENTTYPE][EM_MAXEVENT];
	int CurrentMixEvent;

public:
	MyEventMixer();
	~MyEventMixer();

	void SetVzRange(double VzL, double VzH) { mVzL = VzL; mVzH = VzH; mVzM = 0.5*(mVzL+mVzH); }
	vector<double> GetVzRange() const { vector<double> VzRange{mVzL, mVzH}; return VzRange; }

	// For the input/initialized event
	void InitEvent(int Cent, double Vz, double EventWeight);

	int GetCent() const { return mCent; }
	double GetVz() const { return mVz; }
	double GetEventWeight() const { return mEventWeight; }

	int GetIdCent() const { return iCent; }
	int GetIdVz() const { return iVz; }
	int GetIdEventType() const { return iEventType; }

	// MyEventToMix
	int GetCurrentMixSize() const;
	MyEventToMix GetEventToMix(int iEvent) const;
	void SetEventToMix(MyEventToMix EventToMix);
	void SetEventToMix(int n, int cent, double vz, double weight, int *chg, double *pt, double *phi, double *eta, double *wts);
	void SetEventToMix(int n, int *chg, double *pt, double *phi, double *eta, double *wts);

	void SetEventToMixAddBackup(double x);
};

#endif
