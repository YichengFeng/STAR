//#define CheckPhiCorPlot

#include "MyPhiCorrection.h"
#include "MyPhiCorrection.cxx"
#include "MyPhiCorrectionReg.h"

using namespace std;


void CheckPhiCor() {

	const int nSyst = 9;

	const int nReg = 1;

	TString StrType = "";

	double mVzL = -30;
	double mVzH = +30;

	MyPhiCorrectionReg *mPhiCorReg[nSyst];
	for(int iSyst=0; iSyst<nSyst; iSyst++) {
		TString StrSyst = Form("Syst%d", iSyst);
		mPhiCorReg[iSyst] = new MyPhiCorrectionReg(StrSyst, StrSyst, nReg); // (Dir, Name, nReg)
		mPhiCorReg[iSyst]->GetAll()->SetAveKept(true); // if not save <cos>, <sin> values to speed up and save disk&memory
		mPhiCorReg[iSyst]->GetAll()->LoadPhiEtaWeight("/home/fengyich/Research/BesPPRP/rootfile/Run11AuAu200GeV/SystPbyp/rootfile/prepare.root"); // TODO: TPC phi-eta weights and average values
		mPhiCorReg[iSyst]->GetAll()->LoadVzDEta("/home/fengyich/Research/BesPPRP/rootfile/Run11AuAu200GeV/Pbyp/vzdeta.root"); // TODO: TPC phi-eta weights and average values
		mPhiCorReg[iSyst]->GetAll()->SetWeightCutOff(true, 0.10, 10.0); // false to make it very flat
		for(int iReg=0; iReg<nReg; iReg++) {
			mPhiCorReg[iSyst]->InitReg(iReg);
		}
	}

	vector<double> vz;
	vector<double> deta;
	for(int i=0; i<59; i++) {
		double tmpvz = -30 + (i+0.5);
		mPhiCorReg[0]->GetAll()->InitEvent(5, tmpvz, 1.0);
		vz.push_back(tmpvz);
		deta.push_back(mPhiCorReg[0]->GetAll()->GetDEta());
	}

	TGraph *og = new TGraph(59, &vz[0], &deta[0]);
	TCanvas *cTmp = new TCanvas("Tmp", "Tmp");
	og->Draw("PLA");
	
}
