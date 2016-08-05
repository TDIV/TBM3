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
//  TBModelBase.cpp
//  TBM^3
//


void	TBModelBase::calculateSpinSusceptibility	(TBDataSource & rtbd){
	
	auto isAtom = rtbd.Lat.getAtom(vec(0,0,0));
	
	//if( isAtom.first ){
	//	auto allIndex = isAtom.second.allIndexList();
	//	for( auto & orbitalIndex: allIndex){
	//		//cout<<orbitalIndex.first<<" ";
	//		auto parse = split(orbitalIndex.first, ".");
	//		//for( auto & elem: parse){
	//		//	cout<<elem<<" ";
	//		//}
	//		if( parse[1] == "u")
	//		cout<<endl;
	//		cout<<orbitalIndex.second<<endl;
	//	}
	//}
	
	while(rtbd.Lat.iterate()){
		auto atom = rtbd.Lat.getAtom();
		cout<<atom.atomName<<" "<<atom.pos<<endl;
		for(auto & kvv: rtbd.KEigenValVec){
			cout<<"EVal"<<endl;
			cout<<kvv.eigenValue<<endl;
			cout<<"EVec"<<endl;
			cout<<kvv.eigenVector<<endl;
			
			cout<<endl;
			for(unsigned ii=0 ; ii<kvv.eigenVector.cols() ; ii++){
				for(unsigned jj=0 ; jj<kvv.eigenVector.rows() ; jj++){
					cout<<kvv.eigenVector(ii,jj)<<" ";
				}
				cout<<endl;
			}
		}
		
	}
	
}