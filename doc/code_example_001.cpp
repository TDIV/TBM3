void Hamiltonian(){
	add_Chemical_Potential();
	while( pair_iteration() ){
		add_bond_hc( "Fe:1u  O:1u  +x  t1x");
		add_bond_hc( "Fe:1u  O:1u  +y  t1y");
		add_bond_hc( "Fe:1u  O:1u  +z  t1z");
		add_bond_hc( "Fe:1d  O:1d  +x  t1x");
		add_bond_hc( "Fe:1d  O:1d  +y  t1y");
		add_bond_hc( "Fe:1d  O:1d  +z  t1z");
		add_bond_hc( "Fe:2u  O:1u  +x  t2x");
		add_bond_hc( "Fe:2u  O:1u  +y  t2y");
		add_bond_hc( "Fe:2u  O:1u  +z  t2z");
		add_bond_hc( "Fe:2d  O:1d  +x  t2x");
		add_bond_hc( "Fe:2d  O:1d  +y  t2y");
		add_bond_hc( "Fe:2d  O:1d  +z  t2z");
		add_bond_hc( "Fe:1u  O:1u  -x  t1x");
		add_bond_hc( "Fe:1u  O:1u  -y  t1y");
		add_bond_hc( "Fe:1u  O:1u  -z  t1z");
		add_bond_hc( "Fe:1d  O:1d  -x  t1x");
		add_bond_hc( "Fe:1d  O:1d  -y  t1y");
		add_bond_hc( "Fe:1d  O:1d  -z  t1z");
		add_bond_hc( "Fe:2u  O:1u  -x  t2x");
		add_bond_hc( "Fe:2u  O:1u  -y  t2y");
		add_bond_hc( "Fe:2u  O:1u  -z  t2z");
		add_bond_hc( "Fe:2d  O:1d  -x  t2x");
		add_bond_hc( "Fe:2d  O:1d  -y  t2y");
		add_bond_hc( "Fe:2d  O:1d  -z  t2z");
	}
}

