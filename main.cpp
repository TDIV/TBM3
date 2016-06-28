/*-----------------------------------------------------------|
| Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,            |
|                    Jian-Xin Zhu                            |
|                                                            |
| This file is distributed under the terms of the GNU        |
| General Public License. See the file `LICENSE' in          |
| the root directory of the present distribution, or         |
| http://www.gnu.org/copyleft/gpl.txt .                      |
|                                                            |
|-----------------------------------------------------------*/

//
//  main.cpp
//  TBM^3
//

#include "tbm3.hpp"

#include <iostream>
#include <random>

class TBModel: public tbm::TBModelBase, public tbm::TBOrderIteration{
public:
	TBModel(string filename):
		tbm::TBModelBase(filename),
		tbm::TBOrderIteration(tbd)
		{  }

	void initOrder()	override {
		tbd.order.clear();
		tbd.order.load(Lat.initOrder.orderOperationList);
		tbd.order.save();
	}

	void Hamiltonian()	override {
		
		constructHamList();
	}

	void run()			override {
		
		cout<<endl<<"Starting..."<<endl<<endl;
		
		if( Lat.parameter.VAR("isCalculateMu", 0).real() == 1 ){
			cout<<">> Calculating the chemical potential, Mu."<<endl;
			
			if( Lat.parameter.VAR("disable_quantum", 0).real() != 0 ){
				cout<< "Warning, due to flag 'disable_quantum' turned on."
					<< "The chemical potential calculation will be ignored."<<endl<<endl;
			} else {
				calculateChemicalPotential(true);
				
				TBModelBase::KHamEvd(tbd);
				cout<<endl;
				cout<<"With spin:"<< Lat.parameter.STR("spin")<<endl;
				cout<<"And space:"<< Lat.parameter.STR("space")<<endl;
				cout<<"Total electron count: "<<calculateElectronFilling(tbd);
				cout<<endl<<endl;
			}
		}
		
		if( Lat.parameter.VAR("isCalculateVar", 0).real() == 1 ){
			cout<<">> Calculating LLG spin dynamic."<<endl;
			
			double max_spin_diff = 1;
			
			// Initial iteration, regardless the convergence criteria.
			for( unsigned i=0 ; i<4 ; i++){
				
				TBModelBase::KHamEvd(tbd);
				TBOrderIteration::iterateOrder(tbd.order);
			}
			
			// Full iteration, depends on the convergence criteria.
			while( max_spin_diff > abs(Lat.parameter.VAR("spin_diff", 0).real()) ){
				
				TBModelBase::KHamEvd(tbd);
				max_spin_diff = TBOrderIteration::iterateOrder(tbd.order);
				cout<<"diff:"<< gmt::fformat(max_spin_diff,16)<<" ";
				double TotalE = 0;
				for( auto & iter: tbd.energyMap ){
					cout<< gmt::fformat(iter.first+":", 10)<<" "<< gmt::fformat(iter.second, 10)<<" ";
					TotalE += iter.second;
				}
				cout<< gmt::fformat("Total:", 6) << gmt::fformat(TotalE,6);
				cout<<endl;
			}
			cout<<endl;
		}
		
		if( Lat.parameter.VAR("isCalculateBand", 0).real() == 1 ){
			cout<<">> Calculating the Band structure."<<endl;
			calculateBandStructure(tbd);
			cout<<endl;
		}
		
		cout<<"Finished."<<endl<<endl;
		
	}
};

int main(int argc, char *argv[]) {
	
	vector<string> args;
	for (int i=0; i<argc; i++) args.push_back(string(argv[i]));
	
	// Construct lattice structure(class) from input: xxx.lat.
	vector<string>	operationList;
	set<string>		defaultOperations;
	
	defaultOperations.insert("-ivasp");
	defaultOperations.insert("-iwein2k");
	
	for(unsigned i=1 ; i<argc ; i++){ operationList.push_back(args[i]); }
	
	if( operationList.size() == 0 )	{ operationList.push_back("mos.lat");}
	
	if( defaultOperations.find(operationList[0]) == defaultOperations.end() ){
		string filename = operationList[0];
		
		TBModel model(filename);
		if( operationList.size() == 1){
			model.render();
		}
		
		if( operationList.size() >= 2 ){
			if( operationList[1] == "-expand"){
				cout<<endl<<"Expanding the lattice of '"<<filename<<"' to a larget one :";
				vector<unsigned> N;
				for(unsigned i=2 ; i<operationList.size() ; i++){ N.push_back(tbm::StrToInt(operationList[i])); }
				
				if( N.size() == 1){ model.saveExpandedLattice(N[0], 1, 1); }
				if( N.size() == 2){ model.saveExpandedLattice(N[0], N[1], 1); }
				if( N.size() == 3){ model.saveExpandedLattice(N[0], N[1], N[2]); }
			}
			else if( operationList[1] == "-init"){
				cout<<endl<<"Initialize the order parameter into:"<<filename<<".ord"<<endl<<endl;
				model.initOrder();
			}
			else {
				tbm::ErrorMessage("Error, input operation not found:" +operationList[1]);
			}
			cout<<"Finished."<<endl<<endl;
		}
	}
	else{
		cout<<"I will later implement the following operations:"<<endl;
		cout<<operationList[0]<<endl;
	}
	
	
	return 0;
}



















