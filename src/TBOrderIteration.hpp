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
//  TBOrderIteration.hpp
//  TBM^3
//

class TBOrderIteration{
private:
	TBDataSource & TBD;
	
	vector<boost::tuple<Atom, x_mat, r_var, string> >			HundCouplingList;
	vector<boost::tuple<Atom, x_mat, x_mat, r_var, string> >	SuperExchangeList;
	vector<boost::tuple<Atom, x_mat, x_mat, x_mat, string> >	DMExchangeList;
	
public:
	
	TBOrderIteration(TBDataSource & _tbd): TBD(_tbd)		{ }
	
	void addSuperExchange	(string opt, string svar)		{ // svar === " @:cspin, Jse "
		replaceAll(opt, "\t", " ");
		auto  dummyParser = split(opt, " ");
		if( dummyParser.size() != 1){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		removeSpace(opt);
		auto parser = split(opt, ":");
		if( parser.size() != 3){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto pair = TBD.Lat.getPair(parser[2]);
		
		if( pair.atomI.atomName != parser[0])	return;
		if( pair.atomJ.atomName != parser[1])	return;
		if( !pair.withinRange() )				return;
		
		removeSpace(svar);
		auto svarParser = split(svar, "*");
		auto orderKey = svarParser[0];
		r_var Jsx = 1.0;
		if( svarParser.size() == 2) Jsx = TBD.parseSiteString(pair.atomI, svarParser[1])[0].real();
		
		auto ordS_I = TBD.order.findOrder(pair.atomI, orderKey);
		auto ordS_J = TBD.order.findOrder(pair.atomJ, orderKey);
		if( !ordS_I.first )	return;
		if( !ordS_J.first )	return;
		
		SuperExchangeList.push_back(boost::make_tuple(pair.atomI, ordS_I.second, ordS_J.second, Jsx, orderKey));
	}
	void addDMExchange		(string opt, string svar)		{ // svar === " @:cspin, Jse "
		//replaceAll(opt, "\t", " ");
		auto  dummyParser = split(opt, " ");
		if( dummyParser.size() != 1){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		removeSpace(opt);
		auto parser = split(opt, ":");
		if( parser.size() != 3){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto pair = TBD.Lat.getPair(parser[2]);
		
		if( pair.atomI.atomName != parser[0])	return;
		if( pair.atomJ.atomName != parser[1])	return;
		if( !pair.withinRange() )				return;
		
		removeSpace(svar);
		auto svarParser = split(svar, "*");
		auto orderKey = svarParser[0];
		
		auto ordS_I = TBD.order.findOrder(pair.atomI, orderKey);
		auto ordS_J = TBD.order.findOrder(pair.atomJ, orderKey);
		if( !ordS_I.first )	return;
		if( !ordS_J.first )	return;
		
		x_mat Dij(1,3);
		r_mat bondIJ = pair.bondIJ();
		Dij[0] = bondIJ[0];
		Dij[1] = bondIJ[1];
		Dij[2] = bondIJ[2];
		Dij = Dij * (1/abs(cdot(Dij, Dij)));
		
		if	( svarParser.size() == 2){
			auto tmpXVec = TBD.parseSiteString(pair.atomI, svarParser[1]);
			if		( tmpXVec.size() == 1 ){ Dij = Dij * tmpXVec[0].real(); }
			else if	( tmpXVec.size() == 3 ){ Dij = tmpXVec; }
		}
		
		cout<<Dij<<endl;
		
		DMExchangeList.push_back(boost::make_tuple(pair.atomI, ordS_I.second, ordS_J.second, Dij, orderKey));
	}
	
	void constructHamList()									{
		TBD.initOrder();
		HundCouplingList.clear();
		SuperExchangeList.clear();
		DMExchangeList.clear();
		
		while( TBD.Lat.iterate() ){
			for( auto & iter : TBD.Lat.hamParser.getOperationList("superEx"))	addSuperExchange(iter.first, iter.second);
			for( auto & iter : TBD.Lat.hamParser.getOperationList("dmEx")	)	addDMExchange(iter.first, iter.second);
		}
		
		
	}
	
	void iterateOrder(OrderParameter & newOrder)			{
		
		map<unsigned, pair<x_mat, string> > FieldJH, FieldSE, FieldDM;
		
		// Construct the Field-Term for variational method of Super Exchange.
		for ( auto & elem :SuperExchangeList)	{
			auto atomI = elem.get<0>();
			auto si = elem.get<1>();
			auto sj = elem.get<2>();
			auto Js = elem.get<3>();
			auto optKey = elem.get<4>();
			
			x_mat sumSj(1,3);
			if (si.size() == 3 and sj.size() == 3) sumSj = Js * sj;
			
			if( FieldSE.find(atomI.atomIndex) == FieldSE.end() ) {
				FieldSE[atomI.atomIndex] = make_pair(x_mat(1,3), optKey);
				FieldSE[atomI.atomIndex].first = FieldSE[atomI.atomIndex].first + sumSj;
			}
			else{
				FieldSE[atomI.atomIndex].first = FieldSE[atomI.atomIndex].first + sumSj;
			}
		}
		
		// Construct the Field-Term for variational method of DM Exchange.
		for ( auto & elem :DMExchangeList)	{
			auto atomI = elem.get<0>();
			auto si = elem.get<1>();
			auto sj = elem.get<2>();
			auto d = elem.get<3>();
			auto optKey = elem.get<4>();
			
			x_mat sumDj(1,3);
			if (si.size() == 3 and sj.size() == 3 and d.size() == 3){
				sumDj[0] = - d[1]* sj[2] + d[2]* sj[1] ;
				sumDj[1] = - d[2]* sj[0] + d[0]* sj[2] ;
				sumDj[2] = - d[0]* sj[1] + d[1]* sj[0] ;
			}                           
			
			if( FieldDM.find(atomI.atomIndex) == FieldDM.end() ) {
				FieldDM[atomI.atomIndex] = make_pair(x_mat(1,3), optKey);
				FieldDM[atomI.atomIndex].first = FieldDM[atomI.atomIndex].first + sumDj;
			}
			else{
				FieldSE[atomI.atomIndex].first = FieldSE[atomI.atomIndex].first + sumDj;
			}
		}
	}
	
private:
};























