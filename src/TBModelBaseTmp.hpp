/*-----------------------------------------------------------|
| Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,            |
|                    Jian-Xin Zhu                            |
|                                                            |
| This file is distributed under the terms of the BSD        |
| Berkeley Software Distribution. See the file `LICENSE' in  |
| the root directory of the present distribution.            |
|                                                            |
|-----------------------------------------------------------*/
//
//  TBModelBaseTmp.hpp
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
	
	double Mu = 0;
	double eng_diff = 0.01;
	vector<EigVec>	selectedKVV;
	//
	//for(auto & kvv: rtbd.KEigenValVec){
	//	//cout<<"EVal"<<endl;
	//	//cout<<kvv.eigenValue<<endl;
	//	//cout<<"EVec"<<endl;
	//	//cout<<kvv.eigenVector<<endl;
	//	
	//	//cout<<endl;
	//	for(unsigned ii=0 ; ii<kvv.eigenValue.size() ; ii++){
	//		auto Eng = kvv.eigenValue[ii];
	//		if( abs(Eng - Mu) < eng_diff ){
	//			selectedKVV.push_back(kvv);
	//			break;
	//		}
	//	}
	//}
	//
	
	auto Nb = Lat.parameter.VEC("Nb");
	double max_Nb = 0;
	for(unsigned iNb = 0 ; iNb < Nb.size() ; iNb++){
		if( Nb[iNb].real() > max_Nb ) max_Nb = Nb[iNb].real();
	}
	
	map<double, pair<r_mat, double> > indexedQMap;
	
	for(auto & kvvI: selectedKVV)
	for(auto & kvvJ: selectedKVV){
		
		auto q = kvvI.kPoint - kvvJ.kPoint;
		
		double dictIndex = 0;
		double IncrementNb = 1;
		for( unsigned i=0 ; i< Nb.size() ; i++){
			IncrementNb *= 2*max_Nb;
			dictIndex += IncrementNb * (q[i]+max_Nb);
		}
		
		pair<r_mat, double> previouslyCalculated;
		if( indexedQMap.find(dictIndex) != indexedQMap.end() ){
			previouslyCalculated = indexedQMap[dictIndex];
		}
		
		double valueCalculated = 0 + previouslyCalculated.second;
		auto calculated = make_pair(q, valueCalculated);
		indexedQMap[dictIndex] = calculated;
		
	}
	
}