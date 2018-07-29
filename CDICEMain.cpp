//****************************************************************************
// CDICE
//     C/C++ VERSION OF DICE
//             Economic - Climate Projections
//             Based on GAMS version of DICE2013R
//				(DICE2013Rv2_102213_vanilla_v24b.gms)
//				
//            William D. Nordhaus, Yale University
//				The Climate Casino (2013), Yale University Press
//              http://www.econ.yale.edu/~nordhaus/homepage/index.html
//			 for previous versions, see also
//				DICE 2007:
//				William D. Nordhaus
//				A Question of Balance (2008), Yale University Press
//				DICE1999:
//				William D. Norhaus and Joseph Boyer
//				Warming the World (2000), The MIT Press
//				DICE 1994:
//				William D. Nordhaus,
//				Managing the Global Commons (1994), The MIT Press
//            
//           
//          Translated to C/C++ by Gregory Garner
//            V 1.0  2013 Penn State University
//
//*****************************************************************************/
/* Copyright (C) 2013 Gregory Garner, Klaus Keller, Patrick Reed, Martha Butler, and others

  CDICE is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at your
  option) any later version.
 
  CDICE is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.
 
  You should have received a copy of the GNU Lesser General Public License
  along with the CDICE code.  If not, see <http://www.gnu.org/licenses/>.
 */


//#ifdef CDICEMOD


#include "CDICE.h"

char dummychar;  //debug aid

DICE dice;
DICE *dicePtr = &dice;
DICE dice2;
DICE *dice2Ptr = &dice2;


//=================================================================================
//
// calc_CDICE(DICE *dice)
//
// Run the full model.  
//
// Return()
// - 'dice' structure variables are updated
//=================================================================================
void calc_CDICE(DICE *dice1, DICE *dice2)
{
	
  // Run the model over the dice object
  processModel(dice1);
  
  // Make temporary variables to hold the changes in utility
  double dutil_e;
  double dutil_c;
  int i;
  int t;
  int tsub;
  
  // Loop over the available time steps
  for(t=1; t<dice1->config.nPeriods; t++) {
    
    // Reset the temp dice object
    processModel(dice2);
    
    // Increment the consumption at this time step
    dice2->econ.c[t-1] += 1.0;
    
    // Do the post processing of this run
    postProcess(dice2);
    
    // Store the utility from this run
    dutil_c = dice2->econ.utility - dice1->econ.utility;
    
    // Reset the consumption
    dice2->econ.c[t-1] = dice1->econ.c[t-1];        
    
    // Increment the emissions at this time step
    dice2->econ.e[t-1] += 1.0;
    
    // Run the model on the temporary dice object
    for(tsub=t; tsub<dice1->config.nPeriods; tsub++) {
      calcCarb(dice2, tsub);
      if(dice2->ctrl.use_doeclim) {
				doeclim_dice_ts(dice2, tsub);
      }
      else {
				calcClim(dice2, tsub);
      }
      calcEcon(dice2, tsub);
    }
    
    // Do the post processing of this run
    postProcess(dice2);
    
    // Store the utility from this run
    dutil_e = dice2->econ.utility - dice1->econ.utility;
    
    // Calculate the SCC from these differences in utility
    dice1->econ.scc[t-1] = -1000.0 * (dutil_e / dutil_c);
  }
  
  // Delete the temporary DICE object
  //free_CDICE(dice2);
	
	return;
}

//=================================================================================
//
// readFile(filename, idx1, idx2)
//
// Opens CDICE_Input file from python wrapper and pulls the correct lines of model
// parameters (indicated by idx1 and idx2).
// 
// Returns a vector of doubles for input to CDICE.
//=================================================================================
vector< vector<double> > readFile(string filename, int idx1, int idx2)
{
	ifstream inFile;
	inFile.open(filename.c_str());
	
	vector< vector <double> > params;
	vector<double> param_hold;
	string line;
	double p;
	int i=0;
	int j=0;
	
	while(getline(inFile, line))
	{
		param_hold.clear();
		stringstream inFile(line);
		if (i>=idx1)
		{
			if (i<idx2)
			{
				while(inFile >> p)
				{
					param_hold.push_back(p);
				}
				params.push_back(param_hold);
			}
		}
		i++;
	}
	
	return params;
}


//=================================================================================
//
// main(filename, N)
//
// 
// Main function receives the name of the file containing rows of factor (parameter) sets,
// and N which is the number CDICE parameter sets to do.
//
// The code divides the total number of runs to do evenly among the workers, then runs through
// them in parallel.  The final results are saved in output files, marked by worker's rank.
//
//=================================================================================
int main(int argc, char* argv[])
{
	// Initialize MPI
	MPI::Init(argc,argv);
	
  // Obtain the number of iterations, from the input	
  int N = atof(argv[2]);
  
  // Obtain the CDICE input file name, from the input
  string filename = argv[1];
  
  // Get the number of processes, and the id of this process
  int rank = MPI::COMM_WORLD.Get_rank();
  int nprocs = MPI::COMM_WORLD.Get_size();

  // Determine the chunk which each processor will need to do.
  int count = N/nprocs;
  
  // We use the rank of this processor to determine which chunk of iterations it
  // will perform
  int start;
  int stop;
  int rem = N%nprocs;
  
  if(rank<rem) {
		start = rank*(count+1);
		stop = start+count+1;
  }
  else {
		start = rem*(count+1)+(rank-rem)*count;
		stop = start+count;
  }


/*	// DEBUG -----------------------------------------------------------------------------
	string filename = argv[1];
	string name = "./output/CDICE_DOECLIM_output.txt";
	int stop = 1;
	int start = 0;
	// -----------------------------------------------------------------------------------
*/	
	
  // Load-in the data
  vector< vector <double> > params = readFile(filename, start, stop);
  vector<double> params_hold;
  
  // Initialize the DICE structures
  // Note: Allocating here allows us to test for
  // DOECLIM.  We can then adjust the output filename
  // accordingly.
  allocateDICE(dicePtr);
  allocateDICE(dice2Ptr);
  
  //This code generates the file name based on processor rank
  char rnk [50];
  int l = sprintf(rnk,"%d",rank);
  string name = "CDICE_";
  name += "output_";
  name += rnk;
  name +=".txt";

	
  // This code opens the output file.
  ofstream myfile;
  myfile.open(name.c_str());
  
  // Loop over params
  int param_len = stop - start;
  for(int j=0; j<param_len; j++) {
	  params_hold = params[j];
	  
	  // Populate the parameter values with this current parameter sample
	  populateDICE(dicePtr, params_hold);
	  populateDICE(dice2Ptr, params_hold);
	  
	  // Initialize DICE for this set of parameters (if needed)
	  initTS(dicePtr);
	  initTS(dice2Ptr);
	  if(dicePtr->ctrl.use_doeclim == 1){
	  	doeclim_load_hind_forc(dicePtr);
	  	doeclim_load_hind_forc(dice2Ptr);
		  doeclim_DICE_init(dicePtr);
		  doeclim_DICE_init(dice2Ptr);	  
		}

	  // Run the model
	  calc_CDICE(dicePtr, dice2Ptr);
	  
	  //myfile <<-99<< " " <<start+j<< endl;
	  // Print the index of the factor set
	  
	  // PRINT OUT MODEL VARIABLES-------------------------------------------
	  // These time-series variables will be printed out to stdout forming
	  // space-delimited columns.  The output can be considered a matrix of
	  // size (dicePtr->config.nPeriods x N) where N is the number of
	  // variables being output.  Time increases with each row, so row 1 would
	  // be the initial time-step and row dicePtr->config.nPeriods would be
	  // the last.
  
	  // Here's currently what is being output
	  // Column: Variable
	  // 1			Utility of per-capita consumption for a given period
	  // 2			Population
	  // 3			Discount Factor
	  // 4			NPV Damages
	  // 5			NPV Abatement
	  // 6			Social cost of carbon
	  // 7			Gross Atmospheric temperature
	  // 8			Gross GWP
	  // 9			Cumulative Carbon
  
	  // Time series variables
	  for(int i=0; i<dicePtr->config.nPeriods; i++) {
	  	
	  	// Figure out which temperature variable to use
	  	double this_temp = 0.0;
	  	if(dicePtr->ctrl.use_doeclim == 1){
	  		this_temp = dicePtr->clim.temp[(dicePtr->config.hind_ns -1) + (i * (int)dicePtr->config.tstep)];
	  	}
	  	else {
	  		this_temp = dicePtr->clim.tatm[i];
	  	}
	
			// Print each time step's econ parameters to file
		  myfile << start+j << " "
		  << i << " "
		  << dicePtr->econ.periodu[i] << " "
		  << dicePtr->econ.l[i] << " "
		  << dicePtr->econ.rr[i] << " "
		  << dicePtr->econ.ri[i] << " "
		  << dicePtr->econ.pv_damages[i] << " "
		  << dicePtr->econ.pv_abatecost[i] << " "
		  << dicePtr->econ.abatecost[i] << " "
		  << dicePtr->econ.damages[i] << " "
		  << dicePtr->econ.c[i] << " "
		  << dicePtr->econ.scc[i] << " "
		  << this_temp << " "
		  << dicePtr->carb.mat[i] << " "
		  << dicePtr->carb.forc[i] << " "
		  << dicePtr->econ.cca[i] << " "
		  << dicePtr->lims.cca_up_step << " " << endl;
	  }	
  }

  myfile.close();
  MPI::Finalize();

  return 0;

}


//#endif //CDICEMOD
