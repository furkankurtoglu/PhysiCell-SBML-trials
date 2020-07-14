#include "./internal_viral_dynamics.h" 

using namespace PhysiCell; 

Submodel_Information internal_viral_dynamics_info; 

void internal_virus_model_setup( void )
{
	internal_viral_dynamics_info.name = "internal viral dynamics"; 
	internal_viral_dynamics_info.version = "0.2.0";
	internal_viral_dynamics_info.main_function= internal_virus_model; 

	// what custom data do I need? 
	
	internal_viral_dynamics_info.cell_variables.push_back( "virion" ); // adhered, in process of endocytosis 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated virion" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral RNA" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral protein" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "assembled virion" ); 

	internal_viral_dynamics_info.cell_variables.push_back( "virion uncoating rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated to RNA rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "protein synthesis rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion assembly rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion export rate" ); 

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	internal_viral_dynamics_info.register_model();
	
	return; 
}

void internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled virion" ); 

	static int nUV = pCell->custom_data.find_variable_index( "uncoated virion" ); 
	static int nR  = pCell->custom_data.find_variable_index( "viral RNA" ); 
	static int nP  = pCell->custom_data.find_variable_index( "viral protein" ); 
	
	// copy virions from "internalized variables" to "custom variables"
/*	
	pCell->custom_data[nV_internal] = 
		phenotype.molecular.internalized_total_substrates[nV_external]; 
	// this transfer is now handled in receptor dynamics 
*/		
	pCell->custom_data[nA_internal] = 
		phenotype.molecular.internalized_total_substrates[nA_external]; 
		
	// actual model goes here 

	// uncoat endocytosed virus
	double dV = dt * pCell->custom_data["virion uncoating rate"] * pCell->custom_data[nV_internal] ;
	if( dV > pCell->custom_data[nV_internal] )
	{ dV = pCell->custom_data[nV_internal]; } 
	pCell->custom_data[nV_internal] -= dV; 
	pCell->custom_data[nUV] += dV; 

	// convert uncoated virus to usable mRNA 
	double dR = dt * pCell->custom_data["uncoated to RNA rate"] * pCell->custom_data[nUV]; 
	if( dR > pCell->custom_data[nUV] )
	{ dR = pCell->custom_data[nUV]; }
	pCell->custom_data[nUV] -= dR; 
	pCell->custom_data[nR] += dR; 
	
	// use mRNA to create viral protein 
	double dP = dt * pCell->custom_data["protein synthesis rate"] * pCell->custom_data[nR];
	pCell->custom_data[nP] += dP; 

	// degrade mRNA 
	

	// degrade protein 
	
	
	// assemble virus 
	double dA = dt * pCell->custom_data["virion assembly rate"] * pCell->custom_data[nP]; 
	if( dA > pCell->custom_data[nP] )
	{ dA = pCell->custom_data[nP]; } 
	pCell->custom_data[nP] -= dA; 
	pCell->custom_data[nA_internal] += dA; 
	
	// set export rate 
	
	phenotype.secretion.net_export_rates[nA_external] = 
		pCell->custom_data["virion export rate" ] * pCell->custom_data[nA_internal]; 
 
	// copy data from custom variables to "internalized variables" 
	
/*	
	phenotype.molecular.internalized_total_substrates[nV_external] = 
		pCell->custom_data[nV_internal];
*/		
	phenotype.molecular.internalized_total_substrates[nA_external] = 
		pCell->custom_data[nA_internal];

	
	return; 
}


void simulate_SBML_for_cell(Cell* pCell, Phenotype& phenotype , double dt)
{   
    // SBML indices
	static int SBML_virion = 0;
    static int SBML_uncoated_viron = 1;
	static int SBML_viral_RNA = 2;
    static int SBML_viral_protein = 3;
    static int SBML_assembled_viron = 4;
    
	rrc::RRVectorPtr vptr;
	rrc::RRCDataPtr result;


    // BioFVM indices
    static int i_virion = microenvironment.find_density_index( "virion" ); 
	static int i_assembled_viron = microenvironment.find_density_index( "assembled virion" );
	
    // Internal Amounts
    double internal_virion = phenotype.molecular.internalized_total_substrates[i_virion];
	double internal_assembled_viron = phenotype.molecular.internalized_total_substrates[i_assembled_viron];
    
    // Custom Data indices
    double i_virion_i = pCell->custom_data.find_variable_index( "virion" );
    double i_uncoated_viron_i = pCell->custom_data.find_variable_index( "uncoated virion" );
    double i_viral_RNA_i = pCell->custom_data.find_variable_index( "viral RNA" );
    double i_viral_protein_i = pCell->custom_data.find_variable_index( "viral_protein" );
    double i_assembled_viron_i = pCell->custom_data.find_variable_index( "assembled virion" );
    double cell_volume = phenotype.volume.total;

    
    // Calculating internal concentrations & Updating cell data
    pCell->custom_data[i_virion_i] = internal_virion / cell_volume;
    pCell->custom_data[i_assembled_viron_i] = internal_assembled_viron / cell_volume;
    
    
    //std::cout <<  "Internal Viral Amount: " << internal_virion  << std::endl;
    //std::cout <<  "Internal Viral Concentration: " << pCell->custom_data[i_virion_i]  << std::endl;

    // Geting molecular model structure
    vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	
    // Setting New Values to SBML
/*     vptr->Data[SBML_virion] = pCell->custom_data[i_virion_i];
    vptr->Data[SBML_uncoated_viron] = pCell->custom_data[i_uncoated_viron_i];
    vptr->Data[SBML_viral_RNA] = pCell->custom_data[i_viral_RNA_i];
    vptr->Data[SBML_viral_protein] = pCell->custom_data[i_viral_protein_i];
    vptr->Data[SBML_assembled_viron] = pCell->custom_data[i_assembled_viron_i]; */
    
    rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);
    
    // SBML Simulation
	//result = rrc::simulateEx (pCell->phenotype.molecular.model_rr, 0, 0.01, 2);  // start time, end time, and number of points

    //int idx = (result->RSize - 1) * result->CSize + 1;
    //std::cout << "Cell ID 0) Saving last energy value (cell custom var) = " << result->Data[idx] << std::endl;
/*     for (int idx=0; idx<20; idx++)
    {
        std::cout << idx << ", " << result->Data[8] << std::endl;
    } */

    // Result Indicing!!!!!

    //pCell->custom_data[i_virion_i]  = result->Data[6];
    //std::cout << result->Data[6] << std::endl;
    
/*     pCell->custom_data[i_uncoated_viron_i]  = result->Data[7];
	pCell->custom_data[i_viral_RNA_i]  = result->Data[8];
    pCell->custom_data[i_viral_protein_i] = result->Data[7];
    pCell->custom_data[i_assembled_viron_i] = result->Data[7]; */
    
    // Updating Internal Amount
    
    // phenotype.molecular.internalized_total_substrates[i_virion] = pCell->custom_data[i_virion_i]*cell_volume;
    // phenotype.molecular.internalized_total_substrates[i_assembled_viron] = pCell->custom_data[i_assembled_viron_i]*cell_volume;
}



void simulate_SBML_for_all_cells(void) 
{
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        //std::cout << "simulate_SBML_for_all_cells test" << std::endl;
        simulate_SBML_for_cell((*all_cells)[i], (*all_cells)[i]->phenotype , 0.01);
    }
} 

