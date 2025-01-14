/***************************************************************************
 * Author: Yicheng Feng
 * Email: feng216@purdue.edu, fengyich@outlook.com
 * Date: 2023/08/07
 * Note: Reg stands for regions, can be used for run-by-run, day-by-day, etc.
 ***************************************************************************/

#ifndef MyZdcSmdReg_H
#define MyZdcSmdReg_H

#include <vector>
#include "MyZdcSmd.h"


class MyZdcSmdReg
{
private:
	int nReg;

	MyZdcSmd *mAll;
	MyZdcSmd *mGau; // Gaussian center
	std::vector<MyZdcSmd*> mReg;

	void InitGau() {
		TString StrGau = Form("Gau");
		mGau = new MyZdcSmd(mAll->GetDir()+StrGau, mAll->GetName()+StrGau);
		mGau->SetRefName(mAll->GetName());
		mGau->SetRefDir(mAll->GetDir());
		mGau->SetCenterMod(2);
		mGau->SetEastRefCenter(mAll->GetEastRefCenter());
		mGau->SetWestRefCenter(mAll->GetWestRefCenter());
		mGau->SetVerbose(mAll->GetVerbose());
		mGau->SetMerge(mAll->GetMerge());
		mGau->SetGainMod(2);
		mGau->SetGain(mAll->GetGain());
		mGau->SetIsAutoFill(mAll->GetIsAutoFill());
		mGau->InputFile(mAll->GetInputFilePath());
	}

	void InitReg(int ireg) {
		int iReg = CheckReg(ireg);
		TString StrReg = Form("Reg%d", iReg);
		mReg[iReg] = new MyZdcSmd(mAll->GetDir()+StrReg, mAll->GetName()+StrReg);
		mReg[iReg]->SetRefName(mAll->GetName());
		mReg[iReg]->SetRefDir(mAll->GetDir());
		mReg[iReg]->SetEastRefCenter(mAll->GetEastRefCenter());
		mReg[iReg]->SetWestRefCenter(mAll->GetWestRefCenter());
		mReg[iReg]->SetVerbose(mAll->GetVerbose());
		mReg[iReg]->SetMerge(mAll->GetMerge());
		mReg[iReg]->SetGainMod(2);
		mReg[iReg]->SetGain(mAll->GetGain());
		mReg[iReg]->SetIsAutoFill(mAll->GetIsAutoFill());
		mReg[iReg]->InputFile(mAll->GetInputFilePath());
	}

public:
	MyZdcSmdReg() {
		nReg = 0;
		mAll = new MyZdcSmd();
		mGau = nullptr;
	}
	MyZdcSmdReg(TString name, int nreg) {
		nReg = nreg;
		mReg.clear();
		mAll = new MyZdcSmd(name);
		mGau = nullptr;
		for(int iReg=0; iReg<nReg; iReg++) {
			MyZdcSmd *tmp = nullptr; // = new MyZdcSmd(mName+Form("Reg%d", iReg));
			mReg.push_back(tmp);
		}
	}
	MyZdcSmdReg(TString dir, TString name, int nreg) {
		nReg = nreg;
		mReg.clear();
		mAll = new MyZdcSmd(dir, name);
		mGau = nullptr;
		for(int iReg=0; iReg<nReg; iReg++) {
			MyZdcSmd *tmp = nullptr; // = new MyZdcSmd(mName+Form("Reg%d", iReg));
			mReg.push_back(tmp);
		}
	}

	void Finish(TFile *fOut) {
		mAll->Finish(fOut);
		if(mGau!=nullptr) mGau->Finish(fOut);
		for(int iReg=0; iReg<(int)mReg.size(); iReg++) {
			if(mReg[iReg]!=nullptr) mReg[iReg]->Finish(fOut);
		}
	}

	int CheckReg(int ireg) {
		int iReg = ireg;
		if(iReg<0) {
			iReg = 0;
			std::cout << "MyZdcSmdReg: bad iReg = " << ireg << ", set to " << iReg << std::endl;
		} else if(iReg>=(int)mReg.size()) {
			iReg = mReg.size() - 1;
			std::cout << "MyZdcSmdReg: bad iReg = " << ireg << ", set to " << iReg << std::endl;
		}
		return iReg;
	}

	MyZdcSmd* GetAll() { return mAll; }
	MyZdcSmd* GetGau() {
		if(mGau==nullptr) InitGau();
		return mGau;
	}
	MyZdcSmd* GetReg(int ireg) {
		int iReg = CheckReg(ireg);
		if(mReg[iReg]==nullptr) InitReg(iReg);
		return mReg[iReg];
	}
};


#endif
