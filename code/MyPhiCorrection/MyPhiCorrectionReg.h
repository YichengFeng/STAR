/***************************************************************************
 * Author: Yicheng Feng
 * Email: feng216@purdue.edu, fengyich@outlook.com
 * Date: 2023/08/07
 * Note: Reg stands for regions, can be used for run-by-run, day-by-day, etc.
 ***************************************************************************/

#ifndef MyPhiCorrectionReg_H
#define MyPhiCorrectionReg_H

#include "MyPhiCorrection.h"


class MyPhiCorrectionReg
{
private:
	bool IsAveKept;
	bool IsWeightCutOff;
	double mWtL, mWtH;

	int nReg;
	MyPhiCorrection *mAll;
	std::vector<MyPhiCorrection*> mReg;

public:
	MyPhiCorrectionReg() {
		nReg = 0;
		mAll = new MyPhiCorrection();
	}
	MyPhiCorrectionReg(TString name, int nreg) {
		nReg = nreg;
		mReg.clear();
		mAll = new MyPhiCorrection(name);
		for(int iReg=0; iReg<nReg; iReg++) {
			MyPhiCorrection *tmp = nullptr; // new MyPhiCorrection(mName+Form("Reg%d", iReg));
			mReg.push_back(tmp);
		}
	}
	MyPhiCorrectionReg(TString dir, TString name, int nreg) {
		nReg = nreg;
		mReg.clear();
		mAll = new MyPhiCorrection(dir, name);
		for(int iReg=0; iReg<nReg; iReg++) {
			MyPhiCorrection *tmp = nullptr; // new MyPhiCorrection(mName+Form("Reg%d", iReg));
			mReg.push_back(tmp);
		}
	}

	void Finish(TFile *fOut) {
		mAll->Finish(fOut);
		for(int iReg=0; iReg<(int)mReg.size(); iReg++) {
			if(mReg[iReg]!=nullptr) mReg[iReg]->Finish(fOut);
		}
	}

	int CheckReg(int ireg) const {
		int iReg = ireg;
		if(iReg<0) {
			iReg = 0;
			std::cout << "MyPhiCorrectionReg: bad iReg = " << ireg << ", set to " << iReg << std::endl;
		} else if(iReg>=(int)mReg.size()) {
			iReg = mReg.size() - 1;
			std::cout << "MyPhiCorrectionReg: bad iReg = " << ireg << ", set to " << iReg << std::endl;
		}
		return iReg;
	}

	void InitReg(int ireg) {
		int iReg = CheckReg(ireg);
		TString StrReg = Form("Reg%d", iReg);
		mReg[iReg] = new MyPhiCorrection(mAll->GetDir()+StrReg, mAll->GetName()+StrReg);
		mReg[iReg]->SetRefDir(mAll->GetDir());
		mReg[iReg]->SetRefName(mAll->GetName());
		mReg[iReg]->SetPtEdge(mAll->GetPtEdge());
		mReg[iReg]->SetPhiRange(mAll->GetPhiRange());
		mReg[iReg]->SetEtaRange(mAll->GetEtaRange());
		mReg[iReg]->SetAveKept(mAll->GetAveKept());
		mReg[iReg]->SetWeightCutOff(mAll->GetWeightCutOff());
		mReg[iReg]->LoadPhiEtaWeight(mAll->GetPathPhiEtaWeight());
		mReg[iReg]->LoadVzDEta(mAll->GetPathVzDEta());
	}

	MyPhiCorrection* GetAll() { return mAll; }
	MyPhiCorrection* GetReg(int ireg) {
		int iReg = CheckReg(ireg);
		if(mReg[iReg]==nullptr) InitReg(iReg);
		return mReg[iReg];
	}
};


#endif
