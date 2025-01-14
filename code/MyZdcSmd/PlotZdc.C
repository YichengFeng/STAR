#define MyZdcSmdPlot
#include "MyZdcSmd.h"
#include "MyZdcSmd.cxx"

using namespace std;


void PlotZdc() {
	MyZdcSmd *mZdcSmd = new MyZdcSmd("Syst0", "Syst0");
	mZdcSmd->SetGainMod(3);
	mZdcSmd->InputFile("../../haddsum/output/zdcsmd.root");
}
