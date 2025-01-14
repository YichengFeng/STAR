/***************************************************************************
 * StEpdEpFinder Author: Mike Lisa
 * This file from: Yicheng Feng
 * Date: 2023/08/07
 * Note: Reg stands for regions, can be used for run-by-run, day-by-day, etc.
 *       StEpdEpFinder has I/O rewritten. 
 ***************************************************************************/

#ifndef StEpdEpFinderReg_H
#define StEpdEpFinderReg_H

#include "StEpdEpFinder.h"


class StEpdEpFinderReg
{
private:
	int nReg;
	int nEventTypeBin;
	StEpdEpFinder *mAll;
	std::vector<StEpdEpFinder*> mReg;
	std::vector<int> mRegEventNumber;

public:
	StEpdEpFinderReg() {
		nReg = 0;
		nEventTypeBin = 10;
		mAll = new StEpdEpFinder();
	}
	StEpdEpFinderReg(TString dir, TString name, int nreg, int neventtypebin) {
		nReg = nreg;
		nEventTypeBin = neventtypebin;
		mReg.clear();
		mRegEventNumber.clear();
		mAll = new StEpdEpFinder(dir, name, nEventTypeBin);
		for(int iReg=0; iReg<nReg; iReg++) {
			StEpdEpFinder *tmp = nullptr; // new StEpdEpFinder(mName+Form("Reg%d", iReg));
			mReg.push_back(tmp);
			mRegEventNumber.push_back(0);
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
			std::cout << "StEpdEpFinderReg: bad iReg = " << ireg << ", set to " << iReg << std::endl;
		} else if(iReg>=(int)mReg.size()) {
			iReg = mReg.size() - 1;
			std::cout << "StEpdEpFinderReg: bad iReg = " << ireg << ", set to " << iReg << std::endl;
		}
		return iReg;
	}

	void InitReg(int ireg) {
		int iReg = CheckReg(ireg);
		TString StrReg = Form("Reg%d", iReg);
		mReg[iReg] = new StEpdEpFinder(mAll->GetDir()+StrReg, mAll->GetName()+StrReg, nEventTypeBin);
		mReg[iReg]->SetRefDir(mAll->GetDir());
		mReg[iReg]->SetRefName(mAll->GetName());
		mReg[iReg]->LoadCorrectionFile(mAll->GetCorrectionFilePath());
		mReg[iReg]->SetnMipThreshold(mAll->GetnMipThreshold());
		mReg[iReg]->SetMaxTileWeight(mAll->GetMaxTileWeight());
		mReg[iReg]->SetEpdHitFormat(mAll->GetEpdHitFormat());
		for(int order=1; order<=_EpOrderMax; order++) {
			if(mAll->GetWeightingScheme() == 1) { // weighted by eta
				if(mAll->GetIsEtaWeighted(order)) mReg[iReg]->SetEtaWeights(order, mAll->GetEtaWeights(order));
			}
			if(mAll->GetWeightingScheme() == 2) { // weighted by ring
				if(mAll->GetIsRingWeighted(order)) mReg[iReg]->SetRingWeights(order, mAll->GetRingWeights(order));
			}
		}
	}

	StEpdEpFinder* GetAll() { return mAll; }
	StEpdEpFinder* GetReg(int ireg) {
		int iReg = CheckReg(ireg);
		if(mReg[iReg]==nullptr) InitReg(iReg);
		mRegEventNumber[iReg] ++;
		return mReg[iReg];
	}

	int GetRegEventNumber(int ireg) const {
		int iReg = CheckReg(ireg);
		return mRegEventNumber[iReg];
	}
};


#endif
