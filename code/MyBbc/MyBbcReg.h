/***************************************************************************
 * Author: Yicheng Feng
 * Email: feng216@purdue.edu, fengyich@outlook.com
 * Date: 2023/10/06
 * Note: Reg stands for regions, can be used for run-by-run, day-by-day, etc.
 ***************************************************************************/

#ifndef MyBbcReg_H
#define MyBbcReg_H

#include <vector>
#include "MyBbc.h"


class MyBbcReg
{
private:
	int nReg;

	MyBbc *mAll;
	std::vector<MyBbc*> mReg;

	void InitReg(int ireg) {
		int iReg = CheckReg(ireg);
		TString StrReg = Form("Reg%d", iReg);
		mReg[iReg] = new MyBbc(mAll->GetDir()+StrReg, mAll->GetName()+StrReg);
		mReg[iReg]->SetRefName(mAll->GetName());
		mReg[iReg]->SetRefDir(mAll->GetDir());
		mReg[iReg]->SetVerbose(mAll->GetVerbose());
		mReg[iReg]->SetGainMod(2);
		mReg[iReg]->SetGain(mAll->GetGain());
		mReg[iReg]->SetIsAutoFill(mAll->GetIsAutoFill());
		mReg[iReg]->InputFile(mAll->GetInputFilePath());
	}

public:
	MyBbcReg() {
		nReg = 0;
		mAll = new MyBbc();
	}
	MyBbcReg(TString name, int nreg) {
		nReg = nreg;
		mReg.clear();
		mAll = new MyBbc(name);
		for(int iReg=0; iReg<nReg; iReg++) {
			MyBbc *tmp = nullptr; // = new MyBbc(mName+Form("Reg%d", iReg));
			mReg.push_back(tmp);
		}
	}
	MyBbcReg(TString dir, TString name, int nreg) {
		nReg = nreg;
		mReg.clear();
		mAll = new MyBbc(dir, name);
		for(int iReg=0; iReg<nReg; iReg++) {
			MyBbc *tmp = nullptr; // = new MyBbc(mName+Form("Reg%d", iReg));
			mReg.push_back(tmp);
		}
	}

	void Finish(TFile *fOut) {
		mAll->Finish(fOut);
		for(int iReg=0; iReg<(int)mReg.size(); iReg++) {
			if(mReg[iReg]!=nullptr) mReg[iReg]->Finish(fOut);
		}
	}

	int CheckReg(int ireg) {
		int iReg = ireg;
		if(iReg<0) {
			iReg = 0;
			std::cout << "MyBbcReg: bad iReg = " << ireg << ", set to " << iReg << std::endl;
		} else if(iReg>=(int)mReg.size()) {
			iReg = mReg.size() - 1;
			std::cout << "MyBbcReg: bad iReg = " << ireg << ", set to " << iReg << std::endl;
		}
		return iReg;
	}

	MyBbc* GetAll() { return mAll; }
	MyBbc* GetReg(int ireg) {
		int iReg = CheckReg(ireg);
		if(mReg[iReg]==nullptr) InitReg(iReg);
		return mReg[iReg];
	}
};


#endif
