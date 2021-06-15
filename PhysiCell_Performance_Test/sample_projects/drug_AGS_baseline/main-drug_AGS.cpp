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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 

#include "./custom_modules/custom.h" 
#include "./custom_modules/custom_main.h"

using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	// load and parse settings file(s)
	
	bool XML_status = false; 
	if( argc > 1 )
	{ XML_status = load_PhysiCell_config_file( argv[1] ); }
	else
	{ XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" ); }
	if( !XML_status )
	{ exit(-1); }

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// PNRG setup 
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	setup_microenvironment(); // modify this in the custom code 

	// User parameters
	double time_add_tnf = parameters.ints("time_add_tnf");
	double time_put_tnf = 0;
	double duration_add_tnf = parameters.ints("duration_add_tnf");
	double time_tnf_next = 0;
	double time_remove_tnf = parameters.ints("time_remove_tnf");
	double concentration_tnf = parameters.doubles("concentration_tnf") * microenvironment.voxels(0).volume * 0.000001;
	double time_add_akti = parameters.ints("time_add_akti");
	double time_put_akti = 0;
	double duration_add_akti = parameters.ints("duration_add_akti");
	double time_akti_next = 0;
	double time_remove_akti = parameters.ints("time_remove_akti");
	double concentration_akti = parameters.doubles("concentration_akti") * microenvironment.voxels(0).volume * 0.000001;
	double time_add_taki = parameters.ints("time_add_taki");
	double time_put_taki = 0;
	double duration_add_taki = parameters.ints("duration_add_taki");
	double time_taki_next = 0;
	double time_remove_taki = parameters.ints("time_remove_taki");
	double concentration_taki = parameters.doubles("concentration_taki") * microenvironment.voxels(0).volume * 0.000001;
	double membrane_lenght = parameters.ints("membrane_length");
	double time_add_bcati = parameters.ints("time_add_bcati");
	double time_put_bcati = 0;
	double duration_add_bcati = parameters.ints("duration_add_bcati");
	double time_bcati_next = 0;
	double time_remove_bcati = parameters.ints("time_remove_bcati");
	double concentration_bcati = parameters.doubles("concentration_bcati") * microenvironment.voxels(0).volume * 0.000001;
	// o2 no dirichlet :
	// double time_add_o2 = parameters.ints("time_add_o2");
	// double time_put_o2 = 0;
	// double duration_add_o2 = parameters.ints("duration_add_o2");
	// double time_o2_next = 0;
	// double time_remove_o2 = parameters.ints("time_remove_o2");
	// double concentration_o2 = parameters.doubles("concentration_o2") * microenvironment.voxels(0).volume * 0.000001;
	// o2 change dirichlet :
	double time_o2_change_dirich = parameters.ints("time_o2_change_dirich");
	double new_o2_dirich = parameters.doubles("new_o2_dirich");

	// do small diffusion steps alone to initialize densities
	int k = microenvironment.find_density_index("tnf");
	if ( k >= 0 ) 
		inject_density_sphere(k, concentration_tnf, membrane_lenght);
	for ( int i = 0; i < 25; i ++ )
		microenvironment.simulate_diffusion_decay( 5*diffusion_dt );
	
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	/* Users typically start modifying here. START USERMODS */ 
	
	create_cell_types();
	
	setup_tissue();

	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	display_citations(); 
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}
	
	// main loop 
	
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			/*
			  Custom add-ons could potentially go here. 
			*/			
			// TNF 
			// if ( PhysiCell_globals.current_time >= time_put_tnf )
			// {
			// 	time_tnf_next = PhysiCell_globals.current_time + duration_add_tnf;
			// 	time_put_tnf += time_add_tnf;
			// }

			// if ( PhysiCell_globals.current_time >= time_remove_tnf )
			// {
			// 	int k = microenvironment.find_density_index("tnf");
			// 	if ( k >= 0 )
			// 		remove_density(k);
			// 	time_remove_tnf += PhysiCell_settings.max_time;
			// }

			// if ( PhysiCell_globals.current_time <= time_tnf_next )
			// {
			// 	int k = microenvironment.find_density_index("tnf");
			// 	if ( k >= 0 ) 
			// 		inject_density_sphere(k, concentration_tnf, membrane_lenght);
			// }
				
			// o2 no dirichlet
			// if ( PhysiCell_globals.current_time >= time_put_o2 )
			// {
			// 	time_o2_next = PhysiCell_globals.current_time + duration_add_o2;
			// 	time_put_o2 += time_add_o2;
			// }

			// if ( PhysiCell_globals.current_time >= time_remove_o2 )
			// {
			// 	int k = microenvironment.find_density_index("oxygen");
			// 	if ( k >= 0 )
			// 		remove_density(k);
			// 	time_remove_o2 += PhysiCell_settings.max_time;
			// }

			// if ( PhysiCell_globals.current_time <= time_o2_next )
			// {
			// 	int k = microenvironment.find_density_index("oxygen");
			// 	if ( k >= 0 )
			// 		inject_density_sphere(k, concentration_o2, membrane_lenght);
			// }

			// o2 change dirichlet
			double time_o2_change_dirich_low = time_o2_change_dirich - 0.005;
			double time_o2_change_dirich_high = time_o2_change_dirich + 0.005;
			
			if ( ( PhysiCell_globals.current_time >= time_o2_change_dirich_low ) && ( PhysiCell_globals.current_time <= time_o2_change_dirich_high ) )
			{
				int k = microenvironment.find_density_index("oxygen");
				if ( k >= 0 )
					for( int n=0; n < microenvironment.number_of_voxels() ; n++ )
					{
						microenvironment.update_dirichlet_node(n, k, new_o2_dirich);
					}
			}


			//AKTi
			// if ( PhysiCell_globals.current_time >= time_put_akti )
			// {
			// 	time_akti_next = PhysiCell_globals.current_time + duration_add_akti;
			// 	time_put_akti += time_put_akti;
			// }

			// if ( PhysiCell_globals.current_time >= time_remove_akti )
			// {
			// 	int k = microenvironment.find_density_index("AKTi");
			// 	if ( k >= 0 )
			// 		remove_density(k);
			// 	time_remove_akti += PhysiCell_settings.max_time;
			// }

			// if ( PhysiCell_globals.current_time <= time_akti_next )
			// {
			// 	int k = microenvironment.find_density_index("AKTi");
			// 	if ( k >= 0 )
			// 		inject_density_sphere(k, concentration_akti, membrane_lenght);
			// }
			
			//TAKi
			// if ( PhysiCell_globals.current_time >= time_put_taki )
			// {
			// 	time_taki_next = PhysiCell_globals.current_time + duration_add_taki;
			// 	time_put_taki += time_add_taki;
			// }

			// if ( PhysiCell_globals.current_time >= time_remove_taki )
			// {
			// 	int k = microenvironment.find_density_index("TAKi");
			// 	if ( k >= 0 )
			// 		remove_density(k);
			// 	time_remove_taki += PhysiCell_settings.max_time;
			// }

			// if ( PhysiCell_globals.current_time <= time_taki_next )
			// {
			// 	int k = microenvironment.find_density_index("TAKi");
			// 	if ( k >= 0 ) 
			// 		inject_density_sphere(k, concentration_taki, membrane_lenght);
			// }
			
			//BCATi
			// if ( PhysiCell_globals.current_time >= time_put_bcati )
			// {
			// 	time_bcati_next = PhysiCell_globals.current_time + duration_add_bcati;
			// 	time_put_bcati += time_put_bcati;
			// }

			// if ( PhysiCell_globals.current_time >= time_remove_bcati )
			// {
			// 	int k = microenvironment.find_density_index("BCATi");
			// 	if ( k >= 0 )
			// 		remove_density(k);
			// 	time_remove_bcati += PhysiCell_settings.max_time;
			// }

			// if ( PhysiCell_globals.current_time <= time_bcati_next )
			// {
			// 	int k = microenvironment.find_density_index("BCATi");
			// 	if ( k >= 0 )
			// 		inject_density_sphere(k, concentration_bcati, membrane_lenght);
			// }
			

			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			PhysiCell_globals.current_time += diffusion_dt;
		}

		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}
