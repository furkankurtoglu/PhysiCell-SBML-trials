void simulate_SBML_for_cell(Cell* pCell, Phenotype& phenotype , double dt)
{   
	rrc::RRVectorPtr vptr;
	rrc::RRCDataPtr result;  // start time, end time, and number of points

	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
    {
        // SBML indices
        static int SBML_idx_glucose = 0;
        static int SBML_idx_oxygen = 1;
        static int SBML_idx_energy = 2;
        static int SBML_idx_lactate = 3;
        static int SBML_idx_glutamine = 4;

        // rrc::RRVectorPtr vptr;
        // rrc::RRCDataPtr result;

        // BioFVM indices
        static int i_Oxy = microenvironment.find_density_index( "oxygen" ); 
        static int i_Glu = microenvironment.find_density_index( "glucose" );
        static int i_Glt = microenvironment.find_density_index( "glutamine" );
        static int i_Lac = microenvironment.find_density_index( "lactate" );
        
        // Internal Amounts
        double internal_oxygen = phenotype.molecular.internalized_total_substrates[i_Oxy];
        double internal_glucose = phenotype.molecular.internalized_total_substrates[i_Glu];
        double internal_glutamine = phenotype.molecular.internalized_total_substrates[i_Glt];
        double internal_lactate = phenotype.molecular.internalized_total_substrates[i_Lac];
        //std::cout << internal_oxygen << "," << phenotype.volume.total << std::endl;
        
        // Custom Data indices
        static int i_Oxy_i = pCell->custom_data.find_variable_index( "oxygen_i_conc" );
        static int i_Glu_i = pCell->custom_data.find_variable_index( "glucose_i_conc" );
        static int i_Glt_i = pCell->custom_data.find_variable_index( "glutamine_i_conc" );
        static int i_Lac_i = pCell->custom_data.find_variable_index( "lactate_i_conc" );
        static int energy_vi = pCell->custom_data.find_variable_index( "energy" );
        double cell_volume = phenotype.volume.total;

        // Calculating internal concentrations & Updating cell data
        pCell->custom_data[i_Oxy_i] = internal_oxygen / cell_volume;
        pCell->custom_data[i_Glu_i] = internal_glucose / cell_volume;
        pCell->custom_data[i_Lac_i] = internal_lactate / cell_volume;
        pCell->custom_data[i_Glt_i] = internal_glutamine / cell_volume;
/*         if (pCell->custom_data[i_Oxy_i] < 0)
        { std::cout <<  pCell->custom_data[i_Oxy_i]  << std::endl; }
        if (pCell->custom_data[i_Glu_i] < 0)
        { std::cout <<  pCell->custom_data[i_Glu_i]  << std::endl; }
        if (pCell->custom_data[i_Lac_i] < 0)
        { std::cout <<  pCell->custom_data[i_Lac_i]  << std::endl; }
        if (pCell->custom_data[i_Glt_i] < 0)
        { std::cout <<  pCell->custom_data[i_Glt_i]  << std::endl; }
         */
        // ! NO Energy Update is required !
        //std::cout <<  "Internal Oxygen Amount: " << internal_oxygen  << std::endl;
        //std::cout <<  "Internal Oxygen Concentration: " << pCell->custom_data[i_Oxy_i]  << std::endl;

        //std::cout <<  "Internal Glucose Amount: " << internal_glucose  << std::endl;
        //std::cout <<  "Internal Glucose Concentration: " << pCell->custom_data[i_Glu_i]  << std::endl;

        // Geting molecular model structure
        vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
        
        // Setting New Values to SBML
        vptr->Data[SBML_idx_oxygen] = pCell->custom_data[i_Oxy_i];
        vptr->Data[SBML_idx_glucose] = pCell->custom_data[i_Glu_i];
        vptr->Data[SBML_idx_lactate] = pCell->custom_data[i_Lac_i];
        vptr->Data[SBML_idx_glutamine] = pCell->custom_data[i_Glt_i];
        vptr->Data[SBML_idx_energy] = pCell->custom_data[energy_vi];
        
        rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);
        
        //std::cout << "Before Simulation Glucose: " << vptr->Data[SBML_idx_glucose] << std::endl;
        // SBML Simulation
        result = rrc::simulateEx (pCell->phenotype.molecular.model_rr, 0, 0.01, 2);  // start time, end time, and number of points
        

        
        //std::cout << result->ColumnHeaders[0] << result->Data[6] << std::endl;
        
        //std::cout << "After Simulation Energy: " << result->Data[8] << std::endl;
        
        
    /*     std::cout << "--- after simulation:" << std::endl;
        for (int idx=0; idx<vptr->Count; idx++)
        {
            std::cout << idx << ", " << vptr->Data[idx] << std::endl;
        } */
        /////////---------------Simulation Demo-----------------------///////


        //int idx = (result->RSize - 1) * result->CSize + 1;
        //std::cout << "Cell ID 0) Saving last energy value (cell custom var) = " << result->Data[idx] << std::endl;
    /*     for (int idx=0; idx<20; idx++)
        {
            std::cout << idx << ", " << result->Data[8] << std::endl;
        } */

        // Result Indicing!!!!!



        pCell->custom_data[i_Glu_i]  = result->Data[7];
        pCell->custom_data[i_Oxy_i]  = result->Data[8];
        pCell->custom_data[energy_vi]  = result->Data[9];
        //std::cout << "Energy: " << pCell->custom_data[energy_vi] << std::endl;
        pCell->custom_data[i_Lac_i] = result->Data[10];
        pCell->custom_data[i_Glt_i] = result->Data[11];
        
        
        phenotype.molecular.internalized_total_substrates[i_Glu] = pCell->custom_data[i_Glu_i]*cell_volume;
        phenotype.molecular.internalized_total_substrates[i_Oxy] = pCell->custom_data[i_Oxy_i]*cell_volume;
        phenotype.molecular.internalized_total_substrates[i_Lac] = pCell->custom_data[i_Lac_i]*cell_volume;
        phenotype.molecular.internalized_total_substrates[i_Glt] = pCell->custom_data[i_Glt_i]*cell_volume;
        freeRRCData (result);
    }