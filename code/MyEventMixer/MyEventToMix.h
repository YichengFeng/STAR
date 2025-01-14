#ifndef MyEventToMix_H
#define MyEventToMix_H

#include <iostream>
#include <vector>


class MyEventToMix
{
private:
	static const int EM_MAXTRACK = 1000;

public:
	int N;
	int Cent;
	double Vz;
	double Weight;
	std::vector<double> Backup;
	
	int    Chg[EM_MAXTRACK];
	double Pt [EM_MAXTRACK];
	double Phi[EM_MAXTRACK];
	double Eta[EM_MAXTRACK];
	double Wts[EM_MAXTRACK];
	double Diy[EM_MAXTRACK];
	double Diy2[EM_MAXTRACK];

	MyEventToMix() {
		N = 0;
		Cent = 0;
		Vz = 0;
		Weight = 1.0;
	}

	MyEventToMix(int n, int cent, double vz, double weight, int *chg, double *pt, double *phi, double *eta, double *wts) {
		SetEventTrack(n, cent, vz, weight, chg, pt, phi, eta, wts);
	}

	MyEventToMix(int n, int cent, double vz, double weight, int *chg, double *pt, double *phi, double *eta, double *wts, double *diy) {
		SetEventTrack(n, cent, vz, weight, chg, pt, phi, eta, wts, diy);
	}

	MyEventToMix(int n, int cent, double vz, double weight, int *chg, double *pt, double *phi, double *eta, double *wts, double *diy, double *diy2) {
		SetEventTrack(n, cent, vz, weight, chg, pt, phi, eta, wts, diy, diy2);
	}

	~MyEventToMix() { }

	void SetEventTrack(int n, int cent, double vz, double weight, int *chg, double *pt, double *phi, double *eta, double *wts) {
		N = n<EM_MAXTRACK?n:EM_MAXTRACK;
		Cent = cent;
		Vz = vz;
		Weight = weight;
		for(int i=0; i<N; i++) {
			Chg[i] = chg[i];
			Pt [i] = pt [i];
			Phi[i] = phi[i];
			Eta[i] = eta[i];
			Wts[i] = wts[i];
			Diy[i] = 0;
			Diy2[i] = 0;
		}
		Backup.clear();
	}

	void SetEventTrack(int n, int cent, double vz, double weight, int *chg, double *pt, double *phi, double *eta, double *wts, double *diy) {
		N = n<EM_MAXTRACK?n:EM_MAXTRACK;
		Cent = cent;
		Vz = vz;
		Weight = weight;
		for(int i=0; i<N; i++) {
			Chg[i] = chg[i];
			Pt [i] = pt [i];
			Phi[i] = phi[i];
			Eta[i] = eta[i];
			Wts[i] = wts[i];
			Diy[i] = diy[i];
			Diy2[i] = 0;
		}
		Backup.clear();
	}

	void SetEventTrack(int n, int cent, double vz, double weight, int *chg, double *pt, double *phi, double *eta, double *wts, double *diy, double *diy2) {
		N = n<EM_MAXTRACK?n:EM_MAXTRACK;
		Cent = cent;
		Vz = vz;
		Weight = weight;
		for(int i=0; i<N; i++) {
			Chg[i] = chg[i];
			Pt [i] = pt [i];
			Phi[i] = phi[i];
			Eta[i] = eta[i];
			Wts[i] = wts[i];
			Diy[i] = diy[i];
			Diy2[i] = diy2[i];
		}
		Backup.clear();
	}
};

#endif
