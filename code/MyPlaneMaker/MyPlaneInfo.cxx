#include "MyPlaneInfo.h"
#include <TRandom.h>
#include <iostream>

using namespace std;


void MyPlaneInfo::Reset() {
	WtsQ = 0;
	IsGood = false;
	for(int i=0; i<_NMAX+1; i++) {
		VctQ[i].Set(0.0, 0.0);
		AveVctQ[i].Set(0.0, 0.0);
		RawPsi[i] = 0;
		for(int j=0; j<_KMAX+1; j++) {
			AveCosKPsi[i][j] = 0;
			AveSinKPsi[i][j] = 0;
		}
	}
}


void MyPlaneInfo::ResetEvent() {
	WtsQ = 0;
	IsGood = false;
	for(int i=0; i<_NMAX+1; i++) {
		VctQ[i].Set(0.0, 0.0);
		AveVctQ[i].Set(0.0, 0.0);
		RawPsi[i] = 0;
	}
}


MyPlaneInfo::MyPlaneInfo() {
	Name = "";
	Reset();
}


MyPlaneInfo::MyPlaneInfo(TString name) {
	Name = name;
	Reset();
}


MyPlaneInfo::~MyPlaneInfo() {
}


void MyPlaneInfo::InputShift(const TProfile2D *hcos, const TProfile2D *hsin) {
	for(int i=1; i<_NMAX+1; i++) {
		for(int j=1; j<_KMAX+1; j++) {
			AveCosKPsi[i][j] = hcos->GetBinContent(i, j);
			AveSinKPsi[i][j] = hsin->GetBinContent(i, j);
		}
	}
}


void MyPlaneInfo::OutputShift(TProfile2D *hcos, TProfile2D *hsin) {
	for(int i=1; i<_NMAX+1; i++) {
		for(int j=1; j<_KMAX+1; j++) {
			hcos->Fill(i, j, cos(i*j*RawPsi[i]));
			hsin->Fill(i, j, sin(i*j*RawPsi[i]));
		}
	}
}


void MyPlaneInfo::OutputShift(TProfile2D *hcos, TProfile2D *hsin, double weight) {
	for(int i=1; i<_NMAX+1; i++) {
		for(int j=1; j<_KMAX+1; j++) {
			hcos->Fill(i, j, cos(i*j*RawPsi[i]), weight);
			hsin->Fill(i, j, sin(i*j*RawPsi[i]), weight);
		}
	}
}


void MyPlaneInfo::Add(double phi, double wts) {
	WtsQ += wts;
	for(int i=0; i<_NMAX+1; i++) {
		VctQ[i] += wts*TVector2(cos(i*phi), sin(i*phi));
	}
}


void MyPlaneInfo::Calc() {
	// 0th order
	AveVctQ[0] = VctQ[0]/WtsQ;

	// 1th or higher order
	for(int i=1; i<_NMAX+1; i++) {
		RawPsi[i] = VctQ[i].Phi()/i;
		AveVctQ[i] = VctQ[i]/WtsQ;
		double dpsi = 0;
		for(int j=1; j<_KMAX+1; j++) {
			dpsi += 2.0/i/j*(-AveSinKPsi[i][j]*cos(i*j*RawPsi[i])+AveCosKPsi[i][j]*sin(i*j*RawPsi[i]));
		}
		CorPsi[i] = RawPsi[i]+dpsi;
		double period = 2.0*M_PI/i;
		if(CorPsi[i]<0)       CorPsi[i] += period;
		if(CorPsi[i]>=period) CorPsi[i] -= period;
	}
	IsGood = (WtsQ!=0);
	if(!IsGood) {
		for(int i=1; i<_NMAX+1; i++) {
			double period = 2.0*M_PI/i;
			RawPsi[i] = gRandom->Uniform(0, period);
			CorPsi[i] = RawPsi[i];
		}
	}
}

