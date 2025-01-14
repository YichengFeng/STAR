#include "MyEventMixer.h"
#include <TRandom.h>

using namespace std;


//-----------------------------------------------------------------------------
MyEventMixer::MyEventMixer() {
	mCent = 0;
	mVz = 0;
	mVzL = -70;
	mVzH = +70;
	mVzM = (mVzL + mVzH)/2.;
	mEventWeight = 1.0;

	//iCent = 0;
	//iVz = (EM_NVZ-2)/2+1;
	//iEventType = iVz*EM_NCENT + iCent;
	InitEvent(mCent, mVz, mEventWeight);

	for(int jEventType=0; jEventType<EM_NEVENTTYPE; jEventType++) {
		CurrentMixSize[jEventType] = 0;
	}
	//MyEventToMix mEventToMix[EM_NEVENTTYPE][EM_MAXEVENT];
	CurrentMixEvent = 0;
}

//-----------------------------------------------------------------------------
MyEventMixer::~MyEventMixer() {
}

//-----------------------------------------------------------------------------
void MyEventMixer::InitEvent(int Cent, double Vz, double EventWeight) {
	mCent = Cent;
	iCent = Cent;

	mVz = Vz;
	//if(Vz<mVzL) iVz = 0;
	//else if(Vz<mVzM) iVz = 1;
	//else if(Vz<mVzH) iVz = 2;
	//else iVz = 3;
	iVz = (int)((Vz-mVzL)/(mVzH-mVzL)*EM_NVZ);
	if(iVz<0) iVz = 0;
	if(iVz>=EM_NVZ) iVz = EM_NVZ - 1;

	mEventWeight = EventWeight;

	iEventType = iVz*EM_NCENT + iCent;
}

//-----------------------------------------------------------------------------
int MyEventMixer::GetCurrentMixSize() const {
	return CurrentMixSize[iEventType];
}

//-----------------------------------------------------------------------------
MyEventToMix MyEventMixer::GetEventToMix(int iEvent) const {
	if(iEvent<0) return MyEventToMix();
	if(iEvent>=EM_MAXEVENT) return MyEventToMix();
	if(iEvent>=CurrentMixSize[iEventType]) return MyEventToMix();
	return mEventToMix[iEventType][iEvent];
}

//-----------------------------------------------------------------------------
void MyEventMixer::SetEventToMix(MyEventToMix EventToMix) {
	InitEvent(EventToMix.Cent, EventToMix.Vz, EventToMix.Weight);
	if(CurrentMixSize[iEventType]<0) CurrentMixSize[iEventType] = 0;
	if(CurrentMixSize[iEventType]<EM_MAXEVENT) {
		int iEvent = CurrentMixSize[iEventType];
		mEventToMix[iEventType][iEvent] = EventToMix;
		CurrentMixSize[iEventType] ++;
		CurrentMixEvent = iEvent;
	} else {
		int iEvent = gRandom->Integer(EM_MAXEVENT);
		mEventToMix[iEventType][iEvent] = EventToMix;
		CurrentMixEvent = iEvent;
	}
}

//-----------------------------------------------------------------------------
void MyEventMixer::SetEventToMix(int n, int cent, double vz, double weight, int *chg, double *pt, double *phi, double *eta, double *wts) {
	InitEvent(cent, vz, weight);
	if(CurrentMixSize[iEventType]<0) CurrentMixSize[iEventType] = 0;
	if(CurrentMixSize[iEventType]<EM_MAXEVENT) {
		int iEvent = CurrentMixSize[iEventType];
		mEventToMix[iEventType][iEvent].SetEventTrack(n, cent, vz, weight, chg, pt, phi, eta, wts);
		CurrentMixSize[iEventType] ++;
		CurrentMixEvent = iEvent;
	} else {
		int iEvent = gRandom->Integer(EM_MAXEVENT);
		mEventToMix[iEventType][iEvent].SetEventTrack(n, cent, vz, weight, chg, pt, phi, eta, wts);
		CurrentMixEvent = iEvent;
	}
}

//-----------------------------------------------------------------------------
void MyEventMixer::SetEventToMix(int n, int *chg, double *pt, double *phi, double *eta, double *wts) {
	if(CurrentMixSize[iEventType]<0) CurrentMixSize[iEventType] = 0;
	if(CurrentMixSize[iEventType]<EM_MAXEVENT) {
		int iEvent = CurrentMixSize[iEventType];
		mEventToMix[iEventType][iEvent].SetEventTrack(n, mCent, mVz, mEventWeight, chg, pt, phi, eta, wts);
		CurrentMixSize[iEventType] ++;
		CurrentMixEvent = iEvent;
	} else {
		int iEvent = gRandom->Integer(EM_MAXEVENT);
		mEventToMix[iEventType][iEvent].SetEventTrack(n, mCent, mVz, mEventWeight, chg, pt, phi, eta, wts);
		CurrentMixEvent = iEvent;
	}
}

//-----------------------------------------------------------------------------
void MyEventMixer::SetEventToMixAddBackup(double x) {
	mEventToMix[iEventType][CurrentMixEvent].Backup.push_back(x);
}
