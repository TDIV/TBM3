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
//  import.hpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 9/19/16.
//

/*######################################################
 The following class will combine all the files through
 the #Import blocks.
 ######################################################*/
// Use the TBMImportParser to setup multiple file with the #Import blocks.

class TBMImporter{
public:
	set<string>			tbmFilenameStorage;
	vector<LineStorage>	tbmLineStorage;
	
	TBMImporter(){ }
	~TBMImporter(){
		tbmFilenameStorage.clear();
		tbmLineStorage.clear();
	}
	
	void	append(string filename, string upperLevelFilename)	{
		
		tbmFilenameStorage.insert(filename);
		
        ifstream infile;
		infile.open(filename);
		
		if ( infile.good() ) {
			string line;
			string header = "";
			string flag	= "";
			unsigned lineNumber = 0;
			
            while ( getline(infile, line) ) {
				lineNumber++;
				
				deleteComment(line);
				istringstream iss(line);
				iss >> header;
				
				removeSpace(header);
				
				if	( header == "#Import")	{ flag = header; continue; }
				else if(header[0] == '#')	{ flag = header; }
				
				if	( flag	== "#Import")		{
					removeSpaceTopToe(line);
					if( line[0] != '\"' and line[line.size()-1] != '\"'){
						continue;
					}
					replaceAll(line, "\"", "");

					// Make sure that this filename is not stored inside the pool.
					if( tbmFilenameStorage.find(line) == tbmFilenameStorage.end() ){
						tbmFilenameStorage.insert(line);

						TBMImporter * importer = new TBMImporter();
						importer->append(line, filename);

						for( auto & lineStorage: importer->tbmLineStorage){
							tbmLineStorage.push_back(lineStorage);
						}
						
						for( auto & imp_filename: importer->tbmFilenameStorage){
							tbmFilenameStorage.insert(imp_filename);
						}
					}
					else{
						ErrorMessage("Error: from \""+filename+"\"\n\""
										+ line +"\" is multiple imported.\n"
										+"Please check your file import structure.");
					}
				}
				else{
					tbmLineStorage.push_back(LineStorage(filename, lineNumber, line));
				}
			}
		}
		else{
			ErrorMessage("Error: from "+upperLevelFilename+" \n Cannot find import file: "+ filename);
		}
		
		infile.close();
	}
};


class TBMParser{
public:
	
	set<string>			tbmFilenameStorage;
	
	// Read from xxx.lat.tbm
	Parameter			parameter;
	KSymmetryPoint		kSymmPointParser;
	BondVector			bondVector;
	InitOrder			initOrder;
	CoreCharge			coreCharge;
	LDOSList			ldosList;
	
	HamiltonianParser	hamParser;
	HamiltonianPreprocessor	hamPreprocessor;
	
	TBMParser(){
		line = "";
		header = "";
		flag = "";
		sub_flag = "0";
	}
	
	void	open(string filename, string upperLevelFilename, bool isWithHamiltonianBlock = true)	{
		
		tbmFilenameStorage.insert(filename);
		
        ifstream infile;
		infile.open(filename);
		
		if ( infile.good() ) {
			string line;
			string header = "";
			string flag	= "";
			unsigned mainLineNumber = 0;
			
            while ( getline(infile, line) ) {
				mainLineNumber++;
				
				deleteComment(line);
				istringstream iss(line);
				iss >> header;
				
				removeSpace(header);
				
				auto mainLineStorage = LineStorage(filename, mainLineNumber, line);
				
				if	( header == "#Import")	{ flag = header; continue; }
				else if(header[0] == '#')	{ flag = header; }
				
				if	( flag	== "#Import")		{
					
					removeSpaceTopToe(line);
					if( line[0] != '\"' and line[line.size()-1] != '\"'){
						continue;
					}
					replaceAll(line, "\"", "");

					// Make sure that this filename is not stored inside the pool.
					if( tbmFilenameStorage.find(line) == tbmFilenameStorage.end() ){
						tbmFilenameStorage.insert(line);

						TBMImporter * importer = new TBMImporter();
						importer->append(line, filename);

						for( auto & importLineStorage: importer->tbmLineStorage){
							parseLine(importLineStorage, isWithHamiltonianBlock);
						}
						
						for( auto & imp_filename: importer->tbmFilenameStorage){
							tbmFilenameStorage.insert(imp_filename);
						}
					}
					else{ // If this filename already exist, print an error message.
						ErrorMessage("Error: from \""+filename+"\"\n\""
										+ line +"\" is multiple imported.\n"
										+"Please check your file import structure.");
					}
				}
				else{
					parseLine(mainLineStorage, isWithHamiltonianBlock);
				}
			}
		}
		else{
			ErrorMessage("Error: from "+upperLevelFilename+" \n Cannot find import file: "+ filename);
		}
		
		infile.close();
	}
	
private:

	string line;
	string header = "";
	string flag = "";
	string sub_flag = "";

	void	parseLine(LineStorage & lineStorage, bool withHamiltonianBlock = true){
		
		auto & line = lineStorage.line;
		
		deleteComment(line); // Clean the commented words
		removeSpaceTopToe(line);
		if(line.empty()) return;
		
		istringstream iss(line);
		iss >> header;
		if	( header == parameter())		{flag = header; return;}
		if	( header == ldosList())			{flag = header; return;}
		if	( header == coreCharge())		{flag = header;	return;}
		if	( header == initOrder())		{flag = header;	return;}
		if	( header == kSymmPointParser())	{flag = header; return;}
		if	( header == bondVector())			{
			flag = header;
			iss>>sub_flag;
			return;
		}
		if	( header == hamParser())		{flag = header;	return;}
		
		if	( flag == parameter())			{ parameter.append(line);	return; }
		if	( flag == ldosList())			{ ldosList.append(line);	return; }
		if	( flag == coreCharge())			{ coreCharge.append(line);	return; }
		if	( flag == initOrder())			{ initOrder.append(line);	return; }
		if	( flag == kSymmPointParser())	{ kSymmPointParser.append(line);				return;	}
		if	( flag == bondVector())			{ bondVector.append(StrToInt(sub_flag) ,line);	return;	}
		
		if	( flag == hamParser() and withHamiltonianBlock)		{
			hamParser.append(line);
			hamPreprocessor.append(lineStorage);	return;
		}
		
	}
};


















