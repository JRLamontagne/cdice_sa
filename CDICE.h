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

#ifndef __cdice_h
#define __cdice_h

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <vector>
#include <mpi.h>

using namespace std;
using std::ifstream;


// Structure: Config
// Contains configuration parameters to run an instance of CDICE.
struct Config
{
  int nPeriods;		// Number of time steps
  double tstep;		// Number of years per time step
  int startYear;	// Model start year
  double *dateSeries;	// Human-readable year
  
  int hind_ns;	// Number of hindcast timesteps
  int dt;				// Timestep interval for DOECLIM (years)
  int ns;				// Number of total timesteps for DOECLIM
  int i1900;		// Index of timestep corresponding to the year 1900
  int n_calpoints;		// Number of calibration points

};

// Structure: Econ
// Contains all economic model parameters and variables (both exogenous
// and endogenous). 
struct Econ
{

  // Preferences ---------------------------------
  double elasmu;	// Elasticity of marginal utility of consumption
  double prstp;		// Initial rate of social time preference (per year)


  // Population and Technology -------------------
  double gama;		// Capital elasticity in production function
  double pop0;		// Initial world population [Millions]
  double popadj;	// Growth rate to calibrate to 2050 population projection
  double popasym;	// Asymptotic world population [Millions]
  double dk;		// Depreciation rate on capital (per year)
  double q0;		// Initial world gross output [Trillions 2005 US $]
  double k0;		// Initial capital value [Trillions 2005 US $]
  double a0;		// Initial level of total factor productivity (TFP)
  double ga0;		// Initial growth rate for TFP (per 5 years)
  double dela;		// Decline rate of TFP (per 5 years)

  // Climate damage parameters ---------------------
  double a10;		// Initial damage intercept
  double a20;		// Initial damage quadratic term
  double a1;		// Damage intercept
  double a2;		// Damage quadratic term
  double a3;		// Damage exponent

  // Abatement cost --------------------------------
  double expcost2;	// Exponent of control cost function
  double pback;		// Cost of backstop [2005$ per tCO2 2010]
  double gback;		// Initial cost decline backstop [cost per period]
  double tnopol;	// Period before which no emissions controls base
  double cprice0;	// Initial base carbon price [2005$ per tCO2]
  double gcprice;	// Growth rate of base carbon price (per year)

  // Emissions control and decarbonization
  double gsigma1;	// Initial growth of sigma (per year)
  double dsig;		// Decline rate of decarbonization (per period)
  double sig0;		// Carbon Intensity 2010 [kgCO2 per output 2005 USD 2010]

  // Emissions parameters ------------------------
  double eland0;	// Carbon emissions from land 2010 [GtCO2 per year]
  double deland;	// Decline rate of land emissions (per period)
  double e0;		// Industrial emissions 2010 [GtCO2 per year]
 
  // Scaling and inessential parameters ------------
  // "Note that these are unnecessary for the calculations but are for convenience"
  // Quoted directly from comments in original GAMS code
  double scale1;	// Multiplicitive scaling coefficient
  double scale2;        // Additive scaling coefficient

  // Exogenous timeseries variables
  double *l;		// Level of population and labor (millions)
  double *al;		// Level of total factor productivity
  double *cost1;	// Adjusted cost for backstop
  double *partfract;	// Fraction of emissions in control regime
  //double *gfacpop;	// Growth factor population (Appears in GAMS version, but never used)
  double *pbacktime;	// Backstop price
  double *scc;		// Social cost of carbon
  double *cpricebase;	// Carbon price in base case
  double *photel;	// Carbon price under no damages (Hotelling rent condition)
  double *rr;		// Average utility social discount rate
  double *ga;		// Growth rate of productivity from (sic)
  double *gl;		// Growth rate of labor
  double *gcost;	// Growth of cost factor
  double *gsig;		// Change in sigma (cumulative improvement of energy efficiency)
  double *sigma;	// CO2-equivalent-emissions output ratio 
  double *etree;	// Emissions from deforestation
  double *cumetree;	// Cumulative emissions from land

  // Endogenous variables
  double *c;		// Consumption [Trillions 2005 US$ per year]
  double *k;		// Capital stock [Trillions 2005 US$]
  double *cpc;		// Per capita consumption [Thousands 2005 US$ per year]
  double *i;		// Investment [trillions 2005 US$ per year]
  double *ri;		// Real interest rate (per annum)
  double *y;		// Gross world product net of abatement and damages [Trillions 2005 US$ per year]
  double *ygross;	// Gross world product GROSS of abatement and damages [Trillions 2005 US$ per year]
  double *ynet;		// Output net of damages equation [Trillions of 2005 US$ per year]
  double *damages;	// Damages [Trillions 2005 US$ per year]
  double *damfrac;	// Damages as fraction of gross output
  double *abatecost;	// Cost of emissions reductions [Trillions 2005 US$ per year]
  double *mcabate;	// Marginal cost of abatement [2005 US$ per ton CO2]
  double *periodu;	// One period utility function
  double *cprice;	// Carbon price [2005 US$ per ton CO2]
  double *cemutotper;	// Period utility
  double utility;	// Welfare function (Sum of discounted utility of per capita consumption)
  double *ri_disc;	// Real interest rate (discounted) for present value calculations
  double *pv_damages;	// Present value of damages
  double *pv_abatecost;	// Present value of abatement costs
  double *totalcost;	// Total costs (abatement + damages)
  double *pv_totalcost;	// Present value of total costs (abatement + damages)
  double *utility_vec;	// Vector that contains the discounted utility from stochastic runs
  double *e;		// Total CO2 emissions [GtCO2 per year]
  double *eind;		// Industrial emissions [GtCO2 per year]
  double *cca;		// Cumulative industrial carbon emissions [GtC]
  double *ccatot;	// Total carbon emissions [GtC]
  double *maxcca_vec;	// Vector that contains the maximum cumulative emissions from stochastic runs

};							

// Structure: Carb
// Contains all carbon model parameters and variables (both exogenous
// and endogenous). 
struct Carb
{

  // Carbon cycle ---------------------------------
  // Initial conditions	
  double mat0;		// Initial concentration in atmosphere 2010 [GtC]
  double mu0;		// Initial concentration in upper strata [GtC]
  double ml0;		// Initial concentration in lower strata [GtC]
  double mateq;		// Equilibrium concentration in atmosphere [GtC]
  double mueq;		// Equilibrium concentration in upper strata [GtC]
  double mleq;		// Equilibrium concentration in lower strata [GtC]

  // Flow parameters (Carbon cycle transition matricies)
  double b12;
  double b23;
  double b11;
  double b21;
  double b22;
  double b32;
  double b33;

  // Forcing parameters
  double fco22x;	// Forcings of equilibrium CO2 doubling [Wm-2]
  double fex0;		// 2010 forcings of non-CO2 greenhouse gases (GHG) [Wm-2]
  double fex1;		// 2100 forcings of non-CO2 GHG [Wm-2]

  // Availability of fossil fuels ------------------
  double fosslim;	// Maximum cummulative extraction fossil fuels [GtC]

  // Exogenous timeseries variables
  double *forcoth;	// Exogenous forcing for other GHG

  // Endogenous variables
  double *forc;		// Increase in radiative forcing [Wm-2 from 1900]
  double *mat;		// Carbon concentration increase in atmosphere [GtC from 1750]
  double *mu;		// Carbon concentration increase in shallow oceans [GtC from 1750]
  double *ml;		// Carbon concentration increase in lower oceans [GtC from 1750]
  
  double *atfrac;		// Atmospheric fraction since 1850
  double *atfrac2010;	// Atmospheric fraction since 2010
  double *ppm;		// Atmospheric concentration (actually, mixing ratio)


};

// Structure: Clim
// Contains all climate parameters and variables (both exogenous
// and endogenous). 
struct Clim
{
  // Climate model parameters ---------------------
  double t2xco2;	// Equilibrium temperature impact [dC per doubling CO2]
  double tocean0;	// Initial lower stratum temperature change [dC from 1900]
  double tatm0;		// Initial atmospheric temperature change [dC from 1900]
  double c1;		// Climate equation coefficient for upper level
  double c3;		// Transfer coefficient upper to lower stratum
  double c4;		// Transfer coefficient for lower level
  double lam;		// Climate model parameter

  // DOEclim model parameters
  double t2co;
  double kappa;
  double alpha;
  double *forc;
  double flnd;
  double powtoheat;
  double cal;
  double cas;
  double taucfl;
  double taukls;
  double taucfs;
  double tauksl;
  double bsi;
  double *IB;
  double *A;
  double fso;
  double taudif;
  double *Ker;

  // Endogenous variables ==================================================
  double *tatm;		// Increase temperature of atmosphere [dC from 1900]
  double *tocean;	// Increase temperature of lower oceans [dC from 1900]
  double *maxtemp_vec;	// Vector that contains the maximum temperature increase from stochastic runs

  // DOEclim endogenous variables
  double *temp;
  double *temp_landair;
  double *temp_sst;
  double *heat_mixed;
  double *heat_interior;
  double *heatflux_mixed;
  double *heatflux_interior;
  
  // Variables to hold the calibration data
  double *cs_cal;		// Array to hold the climate sensitivity calibration data
  double *kv_cal;		// Array to hold the ocean diffusivity calibration data
  double *alpha_cal;		// Array to hold the aerosol forcing calibration data

};


// Structure: Dvars
// Contains the decision variables and parameters associated with such.
struct Dvars
{
  double limmiu;	// Upper limit on control rate after 2150
  double miu0;		// Initial emissions control rate for base case 2010

  double *miu;	// Emission control rate GHGs
  double *s;	// Gross savings rate as function of gross world production

};


// Structure: ctrl
// Contains control variables for this model instance
struct Ctrl
{
  int use_doeclim;              // Should the doeclim model be used?
	
};


// Structure: Lims
// Contains the bounds on variables needed for stability
struct Lims
{
	double optlrsav;	// Optimal long-run savings rate used for transversality
  double k_lo;
  double mat_lo;
  double mu_lo;
  double ml_lo;
  double c_lo;
  double tocean_up;
  double tocean_lo;
  double tatm_lo;
  double tatm_up;
  double cpc_lo;
  double y_lo;
  double ygross_lo;
  double i_lo;
  double cca_up;
  int cca_up_step;
  double *miu_up;

};


// Structure: DICE
// The overall structure that defines an instance of CDICE.
struct DICE
{
	
  Config config;
  Dvars dvars;
  Ctrl ctrl;
  Clim clim;
  Econ econ;
  Carb carb;
  Lims lims;

	
};	


// Function prototypes -----------------------
// Functions in CDICEMain:
void calc_CDICE(DICE *dice);		// Sets up a model execution
vector< vector<double> > readFile(string filename, int idx1, int idx2);	// Read input params for CDICE
//int main();

// Functions in CDICEInit:
void allocateCtrl(DICE *dice);	 // Allocate and initialize CTRL structure/contents
void allocateConfig(DICE *dice); // Allocate and initialize CONFIG structure/contents
void allocateCarb(DICE *dice);	 // Allocate and initialize CARB structure/contents
void allocateClim(DICE *dice);	 // Allocate and initialize CLIM structure/contents
void allocateEcon(DICE *dice);	 // Allocate and initialize ECON structure/contents
void allocateDvars(DICE *dice);	 // Allocate and initialize DVARS structure/contents
void allocateLims(DICE *dice);	 // Allocate and initialize LIMS structure/contents
void initTS(DICE *dice);	 // Initialize the time series variables in DICE 
void allocateDICE(DICE *dice);	// Allocates an instance of DICE
void populateDICE(DICE *dice, vector<double> parms);
void free_CDICE(DICE *dice);  // Frees up memory by completely deleting a DICE structure

// Functions in CDICE:
void processModel(DICE *dice);	// Performs one model execution
void calcExog(DICE *dice);		// Calculate exogenous time series
void calcCarb1(DICE *dice);		// Calculate carbon submodel time step1
void calcClim1(DICE *dice);		// Calculate climate submodel time step1
void calcEcon1(DICE *dice);		// Calculate economic submodel time step1
void calcCarb(DICE *dice, int t);	// Calculate carbon submodel time step
void calcClim(DICE *dice, int t);	// Calculate climate model time step
void calcEcon(DICE *dice, int t);	// Calculate economic model time step
void postProcess(DICE *dice);		// Calculate derived final model results

// Functions in CDICE_doeclim:
void count_lines(const char *filename, int *numlines);
void invert_1d_2x2_matrix(double * x, double * y);
void sum_1d_2x2_matrix(double * x, double * y, double * z);
void doeclim_ts(int tstep, int dt, int ns, double *forcing, double t2co, double kappa, 
      double flnd, double powtoheat, double cal, double cas,
      double taucfl, double taukls, double taucfs, double tauksl,
      double bsi, double *IB, double *A, double fso,
      double taudif, double *Ker,
      double *temp, double *temp_landair, double *temp_sst, 
      double *heat_mixed, double *heat_interior, 
      double *heatflux_mixed, double *heatflux_interior);
void doeclim_init(double t2co, double kappa, int dt, int ns, 
		  double *flnd, double *powtoheat, double *cal, double *cas, double *taucfl,
		  double *taukls, double *taucfs, double *tauksl, double *bsi, 
		  double *IB, double *A, double *fso, double *taudif, double *Ker);
void doeclim_DICE_init(DICE *dice);
void doeclim_hindcast(DICE *dice);
void doeclim_dice_ts(DICE *dice, int ts);
void doeclim_load_hind_forc(DICE *dice);
void doeclim_load_calibration(DICE *dice, int n);
void apply_calibration(DICE *dice, double this_cs);

#endif // __cdice_h

//#endif //CDICEMOD
