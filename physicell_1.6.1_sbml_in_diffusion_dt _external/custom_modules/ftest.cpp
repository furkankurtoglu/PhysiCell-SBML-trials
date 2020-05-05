/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/


#include "./ftest.h"

// int oxygen_i, glucose_i, energy_i; 
int oxygen_i, glucose_i; 
int energy_vi; 

// These are for C
// #define STATIC_RRC
#include "rrc_api.h"
#include "rrc_types.h"
// #include "rrc_utilities.h"
extern "C" rrc::RRHandle createRRInstance();

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
// #include <vector>
#include <string>


void create_cell_types( void )
{
	// housekeeping 
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	cell_defaults.type = 0; 
	cell_defaults.name = "cell"; 

	cell_defaults.custom_data.add_variable( "energy" , "dimensionless", 0.0 ); 
	energy_vi = cell_defaults.custom_data.find_variable_index( "energy" ); 
	cell_defaults.functions.update_phenotype = NULL; 
	
    
    cell_defaults.custom_data.add_variable( "oxygen_i_conc" , "mmHg", 0.0 ); 
    cell_defaults.custom_data.add_variable( "glucose_i_conc" , "mol", 0.0 ); 
    cell_defaults.custom_data.add_variable( "lactate_i_conc" , "mol", 0.0 ); 


	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int glucose_substrate_index = microenvironment.find_density_index( "glucose" ); 
	int lactate_substrate_index = microenvironment.find_density_index( "lactate" ); 


	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0.1; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 0.0;
	
	cell_defaults.phenotype.secretion.uptake_rates[glucose_substrate_index] = 0.1;
	cell_defaults.phenotype.secretion.secretion_rates[glucose_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[glucose_substrate_index] = 0.0; 
	
	cell_defaults.phenotype.secretion.uptake_rates[lactate_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[lactate_substrate_index] = 0.3/8/3; 
	cell_defaults.phenotype.secretion.saturation_densities[lactate_substrate_index] = 10.0; 




	// // Set phenotype as you want

	// // turn off birth
	int start_index = live.find_phase_index( PhysiCell_constants::live );
	int end_index = live.find_phase_index( PhysiCell_constants::live );
	cell_defaults.phenotype.cycle.data.transition_rate(start_index,end_index) = 0.0; 

	// // turn off death 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	cell_defaults.phenotype.death.rates[ apoptosis_model_index ] = 0.0; 
	cell_defaults.phenotype.death.rates[ necrosis_model_index ] = 0.0; 

	// // turn off adhesion and repulsion 
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0; 
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = 0.0; 
	
	// // set motilty parameters 
	cell_defaults.phenotype.motility.is_motile = false;

	return; 
}


void setup_microenvironment( void )
{
	// domain parameters read from XML config file

	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// initialize BioFVM 
	initialize_microenvironment(); 	

	oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	glucose_i = microenvironment.find_density_index( "glucose" ); 
    int lactate_i = microenvironment.find_density_index( "lactate" ); 
 
	std::cout << "---------- setup_microenv\n";
	std::cout << "    oxygen_i = " << oxygen_i << std::endl;
	std::cout << "    glucose_i = " << glucose_i << std::endl;
	// std::cout << "    energy_i = " << energy_i << std::endl;

	std::vector<double> dirichlet_conditions(3);
	dirichlet_conditions[oxygen_i]=8.0;
	dirichlet_conditions[glucose_i]=20.0;
	dirichlet_conditions[lactate_i]=0.0;

	std::cout << "setup_microenvironment: num voxels= " << microenvironment.number_of_voxels() << std::endl;
    
	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
        microenvironment.add_dirichlet_node( n,dirichlet_conditions );
	}
    
    microenvironment.set_substrate_dirichlet_activation(1,true);
    microenvironment.set_substrate_dirichlet_activation(2,true);
    microenvironment.set_substrate_dirichlet_activation(3,false);    
	
	return; 
}



void setup_tissue( void )
{
	Cell* pC;
    pC = create_cell(); 
    pC->assign_position( 0.0, 0.0, 0.0 ); 
    pC->set_total_volume( pC->get_total_volume() * 1.0); 

    // Model_t *mm = SBMLDocument_getModel(sbml_doc);
    // std::cout << "mm =" << mm << std::endl;
    // pC->phenotype.molecular.molecular_model = mm;  // assign the intracellular model to each cell

    std::cerr << "------------->>>>>  Creating rrHandle, loadSBML file\n\n";
    rrc::RRHandle rrHandle = createRRInstance();
    if (!rrc::loadSBML (rrHandle, "Toy_Intracellular_Model_for_ext.xml")) {
        std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
    // 	printf ("Error message: %s\n", getLastError());
    // 	getchar ();
    // 	exit (0);
    }
    pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
    int r = rrc::getNumberOfReactions(rrHandle);
    int m = rrc::getNumberOfFloatingSpecies(rrHandle);
    int b = rrc::getNumberOfBoundarySpecies(rrHandle);
    int p = rrc::getNumberOfGlobalParameters(rrHandle);
    int c = rrc::getNumberOfCompartments(rrHandle);

    std::cerr << "Number of reactions = " << r << std::endl;
    std::cerr << "Number of floating species = " << m << std::endl;  // 4
    std::cerr << "Number of boundary species = " << b << std::endl;  // 0
    std::cerr << "Number of compartments = " << c << std::endl;  // 1

    std::cerr << "Floating species names:\n";
    std::cerr << "-----------------------\n";
    std::cerr << stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle)) <<"\n"<< std::endl;


    std::cerr << "Global Parameters names:\n";
    std::cerr << "-----------------------\n";
    std::cerr << stringArrayToString(rrc::getGlobalParameterIds(rrHandle)) <<"\n"<< std::endl;
    // std::cout << "I ARRIVEEEDD HERE" << std::endl;

	return; 
}

void simulate_SBML_for_cell(Cell* pCell, Phenotype& phenotype , double dt)
{   
    // SBML indices
	static int SBML_idx_glucose = 0;
    static int SBML_idx_oxygen = 1;
	static int SBML_idx_energy = 2;
    static int SBML_idx_lactate = 3;

	rrc::RRVectorPtr vptr;
	rrc::RRCDataPtr result;

    // BioFVM indices
    static int i_Oxy = microenvironment.find_density_index( "oxygen" ); 
	static int i_Glu = microenvironment.find_density_index( "glucose" );
	static int i_Lac = microenvironment.find_density_index( "lactate" );
	
    // Internal Amounts
    double internal_oxygen = phenotype.molecular.internalized_total_substrates[i_Oxy];
	double internal_glucose = phenotype.molecular.internalized_total_substrates[i_Glu];
	double internal_lactate = phenotype.molecular.internalized_total_substrates[i_Lac];
    //std::cout << internal_oxygen << "," << phenotype.volume.total << std::endl;
    
    // Custom Data indices
    double i_Oxy_i = pCell->custom_data.find_variable_index( "oxygen_i_conc" );
    double i_Glu_i = pCell->custom_data.find_variable_index( "glucose_i_conc" );
    double i_Lac_i = pCell->custom_data.find_variable_index( "lactate_i_conc" );
    double energy_vi = pCell->custom_data.find_variable_index( "energy" );
    double cell_volume = phenotype.volume.total;
    
    // Calculating internal concentrations & Updating cell data
    pCell->custom_data[i_Oxy_i] = internal_oxygen / cell_volume;
    pCell->custom_data[i_Glu_i] = internal_glucose / cell_volume;
    pCell->custom_data[i_Lac_i] = internal_lactate / cell_volume;
    // ! NO Energy Update is required !
    
    //std::cout <<  "Internal Oxygen Amount: " << internal_oxygen  << std::endl;
    //std::cout <<  "Internal Oxygen Concentration: " << pCell->custom_data[i_Oxy_i]  << std::endl;
    
    std::cout <<  "Internal Glucose Amount: " << internal_glucose  << std::endl;
    std::cout <<  "Internal Glucose Concentration: " << pCell->custom_data[i_Glu_i]  << std::endl;

    // Geting molecular model structure
    vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	
    // Setting New Values to SBML
    vptr->Data[SBML_idx_oxygen] = pCell->custom_data[i_Oxy_i];
    vptr->Data[SBML_idx_glucose] = pCell->custom_data[i_Glu_i];
    vptr->Data[SBML_idx_lactate] = pCell->custom_data[i_Lac_i];
    vptr->Data[SBML_idx_energy] = pCell->custom_data[energy_vi];
    
    rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);
    
    //std::cout << "Before Simulation: " << vptr->Data[SBML_idx_oxygen] << std::endl;
    // SBML Simulation
	result = rrc::simulateEx (pCell->phenotype.molecular.model_rr, 0, 0.01, 2);  // start time, end time, and number of points
    
    //std::cout << "After Simulation Oxygen: " << result->Data[7] << std::endl;
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



    pCell->custom_data[i_Glu_i]  = result->Data[6];
    pCell->custom_data[i_Oxy_i]  = result->Data[7];
	pCell->custom_data[energy_vi]  = result->Data[8];
    // std::cout << "Energy: " << pCell->custom_data[energy_vi] << std::endl;
    pCell->custom_data[i_Lac_i] = result->Data[7];
    
    
    phenotype.molecular.internalized_total_substrates[i_Glu] = pCell->custom_data[i_Glu_i]*cell_volume;
    phenotype.molecular.internalized_total_substrates[i_Oxy] = pCell->custom_data[i_Oxy_i]*cell_volume;
    phenotype.molecular.internalized_total_substrates[i_Lac] = pCell->custom_data[i_Lac_i]*cell_volume;    
    
    
}


void simulate_SBML_for_all_cells(void) 
{
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        //std::cout << "simulate_SBML_for_all_cells test" << std::endl;
        simulate_SBML_for_cell((*all_cells)[i], (*all_cells)[i]->phenotype , 0.01);
    }
} 


std::vector<std::string> energy_coloring_function( Cell* pCell )
{
	// color 0: cytoplasm fill 
	// color 1: outer outline 
	// color 2: nuclear fill 
	// color 3: nuclear outline 
	
	// some bookkeeping 
	// static int energy_i = pCell->custom_data.find_variable_index( "energy" ); 
	
	std::vector< std::string > output( 4, "white" ); 

	if (pCell->ID == 0)
		//std::cout << "--- coloring fn: energy = " <<pCell->custom_data[energy_vi] << std::endl; 

	if (pCell->custom_data[energy_vi] > 1.8)
		output[0] = "rgb(0,255,0)";
	else if (pCell->custom_data[energy_vi] > 1.29)
		output[0] = "rgb(255,0,0)";
	else if (pCell->custom_data[energy_vi] > 0.9)
		output[0] = "rgb(0,0,0)";
	else if (pCell->custom_data[energy_vi] > 0.5)
		output[0] = "rgb(255,255,0)";
	else if (pCell->custom_data[energy_vi] > 0.1)
		output[0] = "rgb(0,255,255)";
	else 
		output[0] = "rgb(255,0,255)";

	// if (pCell->is_out_of_domain)
	// 	output[0] = "rgb(128,128,128)";
	// else if (!pCell->is_movable)
	// 	output[0] = "rgb(0,0,0)";
/*
	if( pCell->phenotype.death.dead == false )
	{
		int red = (int) round( 255.0 * pCell->custom_data[alpha_i] ); 
		int green = (int) round( 255.0 * pCell->custom_data[beta_i] ); 
		int blue = (int) round( 255.0 * pCell->custom_data[resistance_i] ); 
		int grey = (int) round( 255.0 * pCell->custom_data[energy_i] ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", red, green, blue );
		output[2].assign( szTempString ); // nucleus by alpha, beta, resistance "genes"
		sprintf( szTempString , "rgb(%u,%u,%u)", grey, grey, grey );
		output[0].assign( szTempString ); // cyto by energy 
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic ) 
		// Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown // switch to black???
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	*/
	
	return output; 
}