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
//  HamiltonianPreprocessor.hpp
//  TBM^3
//
//  Created by Yuan-Yen Tai on 9/21/16.
//
struct PreprocessorInfo{
	
	string			filename;
	unsigned		lineNumber;
	string			line;
	
	deque<string>	optList;
	deque<string>	varList;
	
	PreprocessorInfo(){}

	PreprocessorInfo(const LineStorage & lineStorage, deque<string> & _optList, deque<string> & _varList){
		filename	= lineStorage.filename;
		lineNumber	= lineStorage.lineNumber;
		line		= lineStorage.line;
	
		optList		= _optList;
		varList		= _varList;
	}
	
	~PreprocessorInfo(){
		filename.clear();
		line.clear();
		optList.clear();
		varList.clear();
	}
};
// Preprocess the #Hamiltonian section
class HamiltonianPreprocessor:	public ParserBase{
public:
	vector<PreprocessorInfo>	list_OrbitalEng;
	vector<PreprocessorInfo>	list_HundSpin;
	vector<PreprocessorInfo>	list_SiteCouple;
	vector<PreprocessorInfo>	list_SiteCoupleHc;
	vector<PreprocessorInfo>	list_HoppingInt;
	vector<PreprocessorInfo>	list_HoppingIntHc;
	vector<PreprocessorInfo>	list_BondCouple;
	vector<PreprocessorInfo>	list_BondCoupleHc;
	vector<PreprocessorInfo>	list_ScreenCoulomb;
	vector<PreprocessorInfo>	list_IntraHubbard;
	vector<PreprocessorInfo>	list_SuperEx;
	vector<PreprocessorInfo>	list_DMEx;
	vector<PreprocessorInfo>	list_ClassicalFieldB;
	vector<PreprocessorInfo>	list_QuantumFieldB;
	
	vector<PreprocessorInfo>	list_PairingS;
	vector<PreprocessorInfo>	list_PairingU;
	vector<PreprocessorInfo>	list_PairingD;
	vector<PreprocessorInfo>	list_PairingT;
	
	HamiltonianPreprocessor(): ParserBase("#Hamiltonian"){}
	
	void	append(const LineStorage & lineStorage){
		
		auto lineSept= split(lineStorage.line, ">");
		
		if( lineSept.size() != 3){
			ErrorMessage(lineStorage.filename,
						 lineStorage.lineNumber,
						 " \""+lineStorage.line+"\"\n"+
						 " is not a valid input formate.");
		}
		
		removeSpace( lineSept[0] );
		removeSpaceTopToe(lineSept[1]);
		replaceAll( lineSept[1], "\t", " ");
		removeSpace( lineSept[2] );
		
		auto optSept = split(lineSept[1], " ");
		auto varSept = split(lineSept[2], "*");
		
		replaceAll( lineSept[1], ":", " ");
		auto optList = split(lineSept[1], " ");
		
		if		(lineSept[0] == "orbital")	{
			/*****************************************
			 * orbital > Fe 1 > 1
			 *****************************************/
			if( optSept.size() != 2){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_OrbitalEng.push_back(PreprocessorInfo(lineStorage, optList, varSept));
		}
		else if	(lineSept[0] == "hundSpin")	{
			/*****************************************
			 * hund > Fe 1 > @:cspin * Jh
			 *****************************************/
			if( optSept.size() != 2){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_HundSpin.push_back(PreprocessorInfo(lineStorage, optList, varSept));
		}
		else if	(lineSept[0] == "intraHubbard")	{
			/*****************************************
			 * intraHubbard > Fe 1 > U
			 *****************************************/
			if( optSept.size() != 2){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_IntraHubbard.push_back(PreprocessorInfo(lineStorage, optList, varSept));
		}
		else if	(lineSept[0] == "site")		{
			/*****************************************
			 * site > Fe 1u:2d		> var
			 * site > Fe px.u:py.d	> var
			 *****************************************/
			if( optSept.size() != 2 or optList.size() != 3){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_SiteCouple.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "siteHc")	{
			/*****************************************
			 * siteHc > Fe 1u:2d		> @:cspin * Jh
			 * siteHc > Fe px.u:py.d	> @:cspin * Jh
			 *****************************************/
			if( optSept.size() != 2 or optList.size() != 3){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_SiteCoupleHc.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "hopping")	{
			/*****************************************
			 * hopping > Fe:O:+1+0+0# 1:2 >  0.5
			 * hopping > Fe:O:+1+0+0# 1:px >  0.5
			 *****************************************/
			if( optSept.size() != 2 or optList.size() != 5){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_HoppingInt.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "hoppingHc"){
			/*****************************************
			 * hoppingHc > Fe:O:+1+0+0# 1:2 >  0.5
			 * hoppingHc > Fe:O:+1+0+0# 1:px >  0.5
			 *****************************************/
			if( optSept.size() != 2 or optList.size() != 5){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_HoppingIntHc.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "bond")		{
			/*****************************************
			 * bond > Fe:O:+1+0+0# 1:2 >  0.5
			 * bond > Fe:O:+1+0+0# 1:px >  0.5
			 *****************************************/
			if( optSept.size() != 2 or optList.size() != 5){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_BondCouple.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "bondHc")	{
			/*****************************************
			 * bondHc > Fe:O:+1+0+0# 1:2 >  0.5
			 * bondHc > Fe:O:+1+0+0# 1:px >  0.5
			 *****************************************/
			if( optSept.size() != 2 or optList.size() != 5){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_BondCoupleHc.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "screenCoulomb"){
			/*****************************************
			 * screenCoulomb > Fe ~2 >  @:den * alpha
			 *****************************************/
			if( optSept.size() != 2 or optList.size() != 2 or optList[1][0] != '~'){
				ErrorMessage(lineStorage.filename, lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			if( varSept.size()	< 2){
				ErrorMessage(lineStorage.filename, lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " should have two input variables.");
			}
			
			list_ScreenCoulomb.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "superEx")	{
			/*****************************************
			 * superEx	> Fe:Fe:+1+0+0 > @:cspin * Jse
			 *****************************************/
			if( optSept.size() != 1 or optList.size() != 3){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_SuperEx.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "dmEx")		{
			/*****************************************
			 * dmEx	> Fe:Fe:+1+0+0 > @:cspin * Jse
			 *****************************************/
			if( optSept.size() != 1 or optList.size() != 3){
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
			
			list_DMEx.push_back(PreprocessorInfo(lineStorage, optList , varSept));
		}
		else if(lineSept[0] == "fieldB")	{
			/*****************************************
			 * fieldB > Fe   > [x,y,z] * Jse
			 * fieldB > Fe 1 > [x,y,z] * Jse
			 *****************************************/
			
			if( optSept.size() == 1 and optList.size() == 1){
				list_ClassicalFieldB.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else if( optSept.size() == 2 and optList.size() == 2 ){
				list_QuantumFieldB.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else{
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
		}
		else if(lineSept[0] == "pairingS")	{
			/*****************************************
			 * pairingS	> Fe 1:1			> 1		% On-site Singlet pairing
			 * pairingS	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Singlet pairing
			 *****************************************/
			
			if( optSept.size() == 2 and optList.size() == 3){
				list_PairingS.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else if( optSept.size() == 2 and optList.size() == 5 ){
				list_PairingS.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else{
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
		}
		else if(lineSept[0] == "pairingU")	{
			/*****************************************
			 * pairingU	> Fe 1:1			> 1		% On-site Singlet pairing
			 * pairingU	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Singlet pairing
			 *****************************************/
			
			if( optSept.size() == 2 and optList.size() == 3){
				list_PairingU.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else if( optSept.size() == 2 and optList.size() == 5 ){
				list_PairingU.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else{
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
		}
		else if(lineSept[0] == "pairingD")	{
			/*****************************************
			 * pairingD	> Fe 1:1			> 1		% On-site Singlet pairing
			 * pairingD	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Singlet pairing
			 *****************************************/
			
			if( optSept.size() == 2 and optList.size() == 3){
				list_PairingD.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else if( optSept.size() == 2 and optList.size() == 5 ){
				list_PairingD.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else{
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
		}
		else if(lineSept[0] == "pairingT")	{
			/*****************************************
			 * pairingT	> Fe 1:1			> 1		% On-site Singlet pairing
			 * pairingT	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Singlet pairing
			 *****************************************/
			
			if( optSept.size() == 2 and optList.size() == 3){
				list_PairingT.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else if( optSept.size() == 2 and optList.size() == 5 ){
				list_PairingT.push_back(PreprocessorInfo(lineStorage, optList , varSept));
			}
			else{
				ErrorMessage(lineStorage.filename,
							 lineStorage.lineNumber,
							 " \""+lineStorage.line+"\"\n"+
							 " is not a valid input formate.");
			}
		}
	}

};














