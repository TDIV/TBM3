/*-----------------------------------------------------------|
| Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,            |
|                    Jian-Xin Zhu                            |
|                                                            |
| This file is distributed under the terms of the BSD        |
| Berkeley Software Distribution. See the file `LICENSE' in  |
| the root directory of the present distribution, or         |
| https://en.wikipedia.org/wiki/BSD_licenses.                |
|                                                            |
|-----------------------------------------------------------*/
//
//  main-tbm-run.cpp
//  TBM^3
//

#include "header.hpp"

#include <iostream>
#include <random>

class TBModel: public tbm::TBModelBase, public tbm::TBClassicalSpinBase{
public:
	TBModel(string filename):
		tbm::TBModelBase(filename),
		tbm::TBClassicalSpinBase(tbd)
		{
			iteration_max = tbd.Lat.parameter.VAR("max_iter", 5000).real();
			iteration_steps = 0;
		}

	void	initOrder()		override	{
		loadTBM();
		Lat.createAtomList("on","normal", false);
		
		tbd.order.clear();
		tbd.order.load(tbm.initOrder.orderOperationList);
		tbd.order.save();
	}
	void	Hamiltonian()	override	{
		constructHamList();
	}
	void	run()			override	{
		// Set the cuda device for calculation.
		cudaSetDevice((int)Lat.parameter.VAR("cudaSetDevice", 0).real());
		
		cout<<endl<<"Starting..."<<endl<<endl;
		
		
		if( Lat.parameter.VAR("isCalculateMu", 0).real() == 1 ){
			
			if( Lat.parameter.VAR("disable_quantum", 0).real() != 0 ){
				cout<< "Warning, due to flag 'disable_quantum' turned on."
					<< "The chemical potential calculation will be ignored."<<endl<<endl;
			} else {
				cout<<">> Calculating the chemical potential, Mu."<<endl;
				TBModelBase::calculateChemicalPotential(true);
				
				TBModelBase::KHamEvd(tbd);
				cout<<endl;
				cout<<"With spin:"<< Lat.parameter.STR("spin")<<endl;
				cout<<"And space:"<< Lat.parameter.STR("space")<<endl;
				cout<<"Total electron count: "<<TBModelBase::calculateElectronFilling(tbd);
				cout<<endl<<endl;
			}
		}
		
		if( Lat.parameter.VAR("isCalculateVar", 0).real() == 1 ){
			calculateVar();
		}
		
		if( Lat.parameter.VAR("isCalculateBand", 0).real() == 1 ){
			cout<<endl<<">> Calculating the Band structure."<<endl;
			calculateBandStructure(tbd, abs( Lat.parameter.VAR("bandPoints", 30).real()) );
			cout<<endl;
		}
	
		if( Lat.parameter.VAR("isCalculateLDOS", 0).real() == 1 ){
			cout<<endl<<">> Calculating LDOS."<<endl;
			calculateLDOS(tbd);
			cout<<endl;
		}
		
		if( Lat.parameter.VAR("isCalculateSpinX", 0).real() == 1 ){
			cout<<endl<<">> Calculating Spin susceptibility."<<endl;
			//calculateSpinSusceptibility(tbd);
			calculateSpinSusceptibility(stbd);
			cout<<endl;
		}
	
		
	}
	
	bool	iterationStepIncr()			{
		iteration_steps += 1;
		return iteration_steps <= iteration_max;
	}
	
	double	calculateDenMeanField()		{
		
		if( Lat.parameter.VAR("disable_quantum", 1).real() == 1 ){ return 0; }
			
		double den_diff = 1;
		double den_diff_bound = abs(Lat.parameter.VAR("den_diff", 0.001).real());
		den_diff_bound = 0.1;
		
		while( den_diff > den_diff_bound and iterationStepIncr() ){
			
			tbd.order.load();
			tbd.order.save("previous");
			KHamEvd(tbd);
			den_diff = iterateDenOrder(tbd.order, Lat.parameter.VAR("den_mix",0.1).real());
			cout<< gmt::fformat(iteration_steps, 5) <<"  Den-diff>> "<< gmt::fformat(den_diff,16)<<" ";
			double TotalE = 0;
			for( auto & iter: tbd.energyMap ){
				if( abs(iter.second) > 0.0000001 ){
					cout<< gmt::fformat(iter.first+":", 7)<<" "<< gmt::fformat(iter.second, 10)<<" ";
					TotalE += iter.second;
				}
			}
			cout<< gmt::fformat("Total:", 7) << gmt::fformat(TotalE,10)<<" ";
			cout<< gmt::fformat("Mu:", 3)<<" "<< gmt::fformat(tbd.Lat.parameter.VAR("Mu").real());
			cout<<endl;
		}
		return den_diff;
	}
	double	calculateSpinVar()			{
		
		tbd.order.load();
		tbd.order.save("previous");
		
		KHamEvd(tbd);
		auto diff = iterateSpinOrder(tbd.order);
		
		cout<< gmt::fformat(iteration_steps, 5) <<" Spin-diff>> "<< gmt::fformat(diff,16)<<" ";
		
		double TotalE = 0;
		for( auto & iter: tbd.energyMap ){
			if( abs(iter.second) > 0.0000001 ){
				cout<< gmt::fformat(iter.first+":", 7)<<" "<< gmt::fformat(iter.second, 10)<<" ";
				TotalE += iter.second;
			}
		}
		cout<< gmt::fformat("Total:", 6)<<" "<< gmt::fformat(TotalE,10)<<" ";
		cout<< gmt::fformat("Mu:", 3)<<" "<< gmt::fformat(tbd.Lat.parameter.VAR("Mu").real());
		cout<<endl;
		return diff;
	}
	void	calculateVar()				{
		
		double spin_diff = 1;
		double den_diff = 1;
		
		// Self-consistent loop.
		while(	(
				spin_diff > abs(Lat.parameter.VAR("spin_diff", 0.001).real())	or
				den_diff > abs(Lat.parameter.VAR("den_diff", 0.001).real())
				)		and
				iterationStepIncr()){
		//while( iterationStepIncr() ){
		
			den_diff	= calculateDenMeanField();
			spin_diff	= calculateSpinVar();
		}
	}
	
private:
	unsigned iteration_steps;
	unsigned iteration_max;
};

int main(int argc, char *argv[]) {
	
	/* -------------------------------------------
	 Organizing/Collect arguments into a container.
	 ---------------------------------------------*/
	vector<string> args;
	for (int i=0; i<argc; i++) args.push_back(string(argv[i]));
	
	vector<string>	operationList;
	for(unsigned i=1 ; i<argc ; i++){ operationList.push_back(args[i]); }

	// If there is no input arguments, give a default input filename.
	if( operationList.size() == 0 )	{
		operationList.push_back("default.lat");
		operationList.push_back("-run");
	}
	
	// If there is no second argument, make it run.
	if( operationList.size() == 1 )	{
		operationList.push_back("-run");
	}
	
	string filename = operationList[0];
	
	// -----------------------------------------------------------------
	/* --Performing several operations of the program.--*/
	// -----------------------------------------------------------------
	// Run the job!
	// -----------------------------------------------------------------
	if		( operationList[1] == "-run")	{
		TBModel model(filename);
		model.render();
	}
	
	// -----------------------------------------------------------------
	// Expand the lattice.
	// -----------------------------------------------------------------
	else if	( operationList[1] == "-expand"){
		cout<<endl<<"Expanding the lattice of '"<<filename<<"' to a larger one :";
		vector<unsigned> N;
		for(unsigned i=2 ; i<operationList.size() ; i++){ N.push_back(tbm::StrToInt(operationList[i])); }
		
		TBModel model(filename);
		if( N.size() == 1){ model.saveExpandedLattice(N[0], 1, 1); }
		if( N.size() == 2){ model.saveExpandedLattice(N[0], N[1], 1); }
		if( N.size() == 3){ model.saveExpandedLattice(N[0], N[1], N[2]); }
	}
	
	// -----------------------------------------------------------------
	// Initialize the order parameter.
	// -----------------------------------------------------------------
	else if	( operationList[1] == "-init")	{
		cout<<endl<<"Initialize the order parameter into:"<<filename<<".ord"<<endl<<endl;
		TBModel model(filename);
		model.initOrder();
	}
	
	// -----------------------------------------------------------------
	// Generate vesta file formate.
	// -----------------------------------------------------------------
	else if	( operationList[1] == "-ovesta"){
		cout<<endl<<"Convert to the VESTA file formate:"<<filename<<".vesta"<<endl<<endl;
		TBModel model(filename);
		model.convertTo_VESTA(operationList);
	}
	
	// -----------------------------------------------------------------
	// Shift all atoms with a given vector.
	// -----------------------------------------------------------------
	else if	( operationList[1] == "-shift"){
		cout<<endl<<"Shift the atom coordinate."<<endl<<endl;
		
		vector<double> shiftXYZ;
		for(unsigned i=2 ; i<operationList.size() ; i++){ shiftXYZ.push_back(tbm::StrToDouble(operationList[i])); }
		
		while( shiftXYZ.size() < 3 ){
			shiftXYZ.push_back(0);
		}
		
		TBModel model(filename);
		model.shiftXYZ(shiftXYZ[0], shiftXYZ[1], shiftXYZ[2]);
	}
	
	// -----------------------------------------------------------------
	// Change the atom name by some selection rule.
	// -----------------------------------------------------------------
	else if	( operationList[1] == "-changeAtom"){
		cout<<endl<<"Change the atom name."<<endl<<endl;
		
		vector<string> changeAtomOptList;
		for(unsigned i=2 ; i<operationList.size() ; i++){ changeAtomOptList.push_back(operationList[i]); }
		
		TBModel model(filename);
		model.changeAtom(changeAtomOptList);
	}
	
	else	{
		tbm::ErrorMessage("Error, input operation not found:" +operationList[1]);
	}
	cout<<"Finished."<<endl<<endl;

	
	return 0;
}







