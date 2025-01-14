#ifndef MYPLANEINFO_H
#define MYPLANEINFO_H

#define _NMAX 3
#define _KMAX 20

#include <TVector2.h>
#include <TString.h>
#include <TProfile2D.h>
#include <iostream>


class MyPlaneInfo
{
private:
	TString Name;

	double WtsQ;
	TVector2 VctQ[_NMAX+1];
	TVector2 AveVctQ[_NMAX+1];
	bool IsGood;
	double RawPsi[_NMAX+1];
	double CorPsi[_NMAX+1];
	double AveCosKPsi[_NMAX+1][_KMAX+1];
	double AveSinKPsi[_NMAX+1][_KMAX+1];

public:
	MyPlaneInfo();
	MyPlaneInfo(TString name);
	~MyPlaneInfo();

	void Reset();
	void ResetEvent();

	void Add(double phi, double wts);
	void InputShift(const TProfile2D *hcos, const TProfile2D *hsin);
	void OutputShift(TProfile2D *hcos, TProfile2D *hsin);
	void OutputShift(TProfile2D *hcos, TProfile2D *hsin, double weight);
	void Calc();

	void SetName(TString name);
	TString GetName() const;

	double GetWtsQ() const;
	TVector2 GetVctQ(int i) const;
	double GetQx(int i) const;
	double GetQy(int i) const;

	TVector2 GetAveVctQ(int i) const;
	double GetAveQx(int i) const;
	double GetAveQy(int i) const;
	bool GetIsGood() const;

	double GetRawPsi(int i) const;
	double GetCorPsi(int i) const;
};


inline void MyPlaneInfo::SetName(TString name) { Name = name; }
inline TString MyPlaneInfo::GetName() const { return Name; }

inline double MyPlaneInfo::GetWtsQ() const { return WtsQ; }

inline TVector2 MyPlaneInfo::GetVctQ(int i) const {
	if(i>=1 && i<=_NMAX) {
		return VctQ[i];
	} else {
		std::cout << "ERROR: invalid Q-order!" << std::endl;
		return TVector2(0.0, 0.0);
	}
}
inline double MyPlaneInfo::GetQx(int i) const { return GetVctQ(i).X(); }
inline double MyPlaneInfo::GetQy(int i) const { return GetVctQ(i).Y(); }

inline TVector2 MyPlaneInfo::GetAveVctQ(int i) const {
	if(i>=1 && i<=_NMAX) {
		return AveVctQ[i];
	} else {
		std::cout << "ERROR: invalid Q-order!" << std::endl;
		return TVector2(0.0, 0.0);
	}
}
inline double MyPlaneInfo::GetAveQx(int i) const { return GetAveVctQ(i).X(); }
inline double MyPlaneInfo::GetAveQy(int i) const { return GetAveVctQ(i).Y(); }

inline bool MyPlaneInfo::GetIsGood() const { return IsGood; }

inline double MyPlaneInfo::GetRawPsi(int i) const {
	if(i>=1 && i<=_NMAX) {
		return RawPsi[i];
	} else {
		std::cout << "ERROR: invalid Q-order!" << std::endl;
		return 0;
	}
}

inline double MyPlaneInfo::GetCorPsi(int i) const {
	if(i>=1 && i<=_NMAX) {
		return CorPsi[i];
	} else {
		std::cout << "ERROR: invalid Q-order!" << std::endl;
		return 0;
	}
}


#endif
