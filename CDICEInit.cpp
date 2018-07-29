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


//===========================================================================================
// allocateConfig(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Configuration structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateConfig(DICE *dice)
{
    // Allocate the configuration structure----------------------------------------------------

  dice->config.nPeriods = 100;	// Number of time steps
  dice->config.tstep = 5.0;	// Number of years per time step
  dice->config.startYear = 2015;	// Model start year	
  dice->config.dateSeries = new double[dice->config.nPeriods];	// Human-readable year

	if(dice->ctrl.use_doeclim == 1) {
	  dice->config.hind_ns = 116; // Number of hindcast timesteps (1900 - 2010)
  	dice->config.i1900 = 0;  // Index of year 1900
  	dice->config.dt = 1;
  	dice->config.ns = dice->config.hind_ns+(dice->config.tstep*(dice->config.nPeriods-1)); // 1-yr res.
  	dice->config.n_calpoints = 100;
  }

  return;
}


//===========================================================================================
// allocateCtrl(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Control structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateCtrl(DICE *dice)
{
  // Allocate the control variable structure-------------------------------------------------
  // DEFAULTS
  dice->ctrl.use_doeclim = 1;			// Use the DOEClim model (1 = DOEClim, 0 = DICE)

  return;
}


//===========================================================================================
// allocateCarb(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Carbon structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateCarb(DICE *dice)
{

  dice->carb.forcoth = new double[dice->config.nPeriods];	// Exogenous forcing for other GHG
  dice->carb.forc = new double[dice->config.nPeriods];		// Increase in radiative forcing [Wm-2 from 1900]
  dice->carb.mat = new double[dice->config.nPeriods];		// Carbon concentration increase in atmosphere [GtC from 1750]
  dice->carb.mu = new double[dice->config.nPeriods];		// Carbon concentration increase in shallow oceans [GtC from 1750]
  dice->carb.ml = new double[dice->config.nPeriods];		// Carbon concentration increase in lower oceans [GtC from 1750]
  
  dice->carb.atfrac = new double[dice->config.nPeriods];
  dice->carb.atfrac2010 = new double[dice->config.nPeriods];
  dice->carb.ppm = new double[dice->config.nPeriods];


  return;
}

//===========================================================================================
// allocateClim(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Climate structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateClim(DICE *dice)
{
  
  if(dice->ctrl.use_doeclim==1) {
    // DOEclim model variables
    dice->clim.forc = new double[dice->config.ns];
    dice->clim.IB = new double[4];
    dice->clim.A = new double[4];
    dice->clim.Ker = new double[dice->config.ns];
    
    dice->clim.temp = new double[dice->config.ns];
    dice->clim.temp_landair = new double[dice->config.ns];
    dice->clim.temp_sst = new double[dice->config.ns];
    dice->clim.heat_mixed = new double[dice->config.ns];
    dice->clim.heat_interior = new double[dice->config.ns];
    dice->clim.heatflux_mixed = new double[dice->config.ns];
    dice->clim.heatflux_interior = new double[dice->config.ns];
    
    dice->clim.cs_cal = new double [dice->config.n_calpoints];
    dice->clim.kv_cal = new double [dice->config.n_calpoints];
    dice->clim.alpha_cal = new double [dice->config.n_calpoints];
    
		//doeclim_DICE_init(dice);

  }

  else {
    // DICE Climate model variables
    dice->clim.tatm = new double[dice->config.nPeriods];		// Increase temperature of atmosphere [dC from 1900]
    dice->clim.tocean = new double[dice->config.nPeriods];	// Increase temperature of lower oceans [dC from 1900]
    
  }

  return;
}


//===========================================================================================
// allocateEcon(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Economics structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateEcon(DICE *dice)
{
  // Allocate the economic model structure---------------------------------------------
  dice->econ.l = new double[dice->config.nPeriods];		// Level of population and labor (millions)
  dice->econ.al = new double[dice->config.nPeriods];		// Level of total factor productivity
  dice->econ.rr  = new double [dice->config.nPeriods];		// Average utility social discount rate
  dice->econ.ga = new double[dice->config.nPeriods];		// Growth rate of productivity from (sic)
  dice->econ.gl = new double[dice->config.nPeriods];		// Growth rate of labor
  dice->econ.gcost = new double[dice->config.nPeriods];		// Growth of cost factor
  dice->econ.cost1 = new double[dice->config.nPeriods];		// Adjusted cost for backup
  dice->econ.partfract = new double[dice->config.nPeriods];	// Fraction of emissions in control regime
  //dice->econ.gfacpop = new double[dice->config.nPeriods];	// Growth factor population (Appears in GAMS version, but never used
  dice->econ.pbacktime = new double[dice->config.nPeriods];	// Backstop price
  dice->econ.scc = new double[dice->config.nPeriods];		// Social cost of carbon
  dice->econ.cpricebase = new double[dice->config.nPeriods];	// Carbon price in base case
  dice->econ.photel = new double[dice->config.nPeriods];	// Carbon price under no damages (Hotelling rent condition)
  dice->econ.gsig = new double[dice->config.nPeriods];		// Change in sigma (cumulative improvement of energy efficiency)
  dice->econ.sigma = new double[dice->config.nPeriods];		// CO2-equivalent-emissions output ratio 
  dice->econ.etree = new double[dice->config.nPeriods];		// Emissions from deforestation
  dice->econ.cumetree = new double[dice->config.nPeriods];

  dice->econ.c = new double[dice->config.nPeriods];		// Consumption [Trillions 2005 US$ per year]
  dice->econ.k = new double[dice->config.nPeriods];		// Capital stock [Trillions 2005 US$]
  dice->econ.cpc = new double[dice->config.nPeriods];		// Per capita consumption [Thousands 2005 US$ per year]
  dice->econ.i = new double[dice->config.nPeriods];		// Investment [trillions 2005 US$ per year]
  dice->econ.ri = new double[dice->config.nPeriods];		// Real interest rate (per annum)
  dice->econ.y = new double[dice->config.nPeriods];		// Gross world product net of abatement and damages [Trillions 2005 US$ per year]
  dice->econ.ygross = new double[dice->config.nPeriods];	// Gross world product GROSS of abatement and damages [Trillions 2005 US$ per year]
  dice->econ.ynet = new double[dice->config.nPeriods];		// Output net of damages equation [Trillions of 2005 US$ per year]
  dice->econ.damages = new double[dice->config.nPeriods];	// Damages [Trillions 2005 US$ per year]
  dice->econ.damfrac = new double[dice->config.nPeriods];	// Damages as fraction of gross output
  dice->econ.abatecost = new double[dice->config.nPeriods];	// Cost of emissions reductions [Trillions 2005 US$ per year]
  dice->econ.mcabate = new double[dice->config.nPeriods];	// Marginal cost of abatement [2005 US$ per ton CO2]
  dice->econ.periodu = new double[dice->config.nPeriods];	// One period utility function
  dice->econ.cprice = new double[dice->config.nPeriods];	// Carbon price [2005 US$ per ton CO2]
  dice->econ.cemutotper = new double[dice->config.nPeriods];	// Period utility
  dice->econ.ri_disc = new double[dice->config.nPeriods];	// Real interest rate (discounted) for present value calculations
  dice->econ.pv_damages = new double[dice->config.nPeriods];	// Present value of damages
  dice->econ.pv_abatecost = new double[dice->config.nPeriods];	// Present value of abatement costs
  dice->econ.totalcost = new double[dice->config.nPeriods];	// Total costs (abatement + damages)
  dice->econ.pv_totalcost = new double[dice->config.nPeriods];	// Present value of total costs (abatement + damages)
  dice->econ.e = new double[dice->config.nPeriods];		// Total CO2 emissions [GtCO2 per year]
  dice->econ.eind = new double[dice->config.nPeriods];		// Industrial emissions [GtCO2 per year]
  dice->econ.cca = new double[dice->config.nPeriods];		// Cumulative industrial carbon emissions [GtC]
  dice->econ.ccatot = new double[dice->config.nPeriods];

  return;
}

//===========================================================================================
// allocateDvars(DICE *dice) 
//
// Allocates space and initializes parameters/variables for the Decision Variable structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateDvars(DICE *dice)
{

  dice->dvars.miu = new double[dice->config.nPeriods];	// Emission control rate GHGs
  dice->dvars.s = new double[dice->config.nPeriods];	// Gross savings rate as function of gross world production
  

  return;
}


//===========================================================================================
// allocateLims(DICE *dice) 
//
// Allocates space and initializes the bounds in the Lims structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateLims(DICE *dice)
{
  dice->lims.k_lo = 1.0;
  dice->lims.mat_lo = 10.0;
  dice->lims.mu_lo = 100.0;
  dice->lims.ml_lo = 1000.0;
  dice->lims.c_lo = 2.0;
  dice->lims.tocean_up = 20.0;
  dice->lims.tocean_lo = -1.0;
  dice->lims.tatm_lo = 0.0;
  dice->lims.tatm_up = 40.0;
  dice->lims.cpc_lo = 0.01;
  dice->lims.y_lo = 0.0;
  dice->lims.ygross_lo = 0.0;
  dice->lims.i_lo = 0.0;
  dice->lims.miu_up = new double[dice->config.nPeriods];	// Upper limit of control policy for each time period

}


//===========================================================================================
// initTS(DICE *dice) 
//
// Initializes the time series variables in the DICE structure.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void initTS(DICE *dice)
{

  for (int i=0; i<dice->config.nPeriods; i++) {
    
    dice->config.dateSeries[i] = dice->config.startYear + (dice->config.tstep * i);
    dice->econ.l[i] = 0.0;
    dice->econ.al[i] = 0.0;		
    dice->econ.sigma[i] = 0.0;	 
    dice->econ.rr[i] = 0.0;		
    dice->econ.ga[i] = 0.0;		
    dice->carb.forcoth[i] = 0.0;
    dice->econ.gl[i] = 0.0;		
    dice->econ.gcost[i] = 0.0;
    dice->econ.gsig[i] = 0.0;		
    dice->econ.etree[i] = 0.0;
    dice->econ.cumetree[i] = 0.0;
    dice->econ.cost1[i] = 0.0;
    //dice->econ.gfacpop[i] = 0.0;	// Appears in GAMS version, but is never used
    dice->econ.pbacktime[i] = 0.0;
    dice->econ.scc[i] = 0.0;
    dice->econ.cpricebase[i] = 0.0;
    dice->econ.photel[i] = 0.0;

    //dice->dvars.miu[i] = dice->dvars.miu0;
    dice->carb.forc[i] = 0.0;
    dice->carb.mat[i] = 0.0;
    dice->carb.mu[i] = 0.0;
    dice->carb.ml[i] = 0.0;
    dice->carb.atfrac[i] = 0.0;
    dice->carb.atfrac2010[i] = 0.0;
    dice->carb.ppm[i] = 0.0;
    dice->econ.e[i] = 0.0;
    dice->econ.eind[i] = 0.0;
    dice->econ.c[i] = 0.0;
    dice->econ.k[i] = 0.0;
    dice->econ.cpc[i] = 0.0;
    dice->econ.i[i] = 0.0;
    dice->econ.ri[i] = 0.0;
    dice->econ.y[i] = 0.0;
    dice->econ.ygross[i] = 0.0;
    dice->econ.ynet[i] = 0.0;
    dice->econ.damages[i] = 0.0;
    dice->econ.damfrac[i] = 0.0;
    dice->econ.abatecost[i] = 0.0;
    dice->econ.mcabate[i] = 0.0;
    dice->econ.cca[i] = 0.0;
    dice->econ.ccatot[i] = 0.0;
    dice->econ.periodu[i] = 0.0;
    dice->econ.cprice[i] = 0.0;
    dice->econ.cemutotper[i] = 0.0;

    dice->econ.ri_disc[i] = 0.0;
    dice->econ.pv_damages[i] = 0.0;
    dice->econ.pv_abatecost[i] = 0.0;
    dice->econ.totalcost[i] = 0.0;
    dice->econ.pv_totalcost[i] = 0.0;

    // If not using DOEclim, then initialize DICE clim time series
    if(dice->ctrl.use_doeclim == 0) {
      dice->clim.tatm[i] = 0.0;
      dice->clim.tocean[i] = 0.0;
    }
    
  }
	
	// If using DOEclim...
	if(dice->ctrl.use_doeclim == 1) {

		// Initialize the appropriate time series
		for(int i=0; i<dice->config.ns; i++) {
  	  dice->clim.temp[i] = 0.0;
    	dice->clim.temp_landair[i] = 0.0;
	    dice->clim.temp_sst[i] = 0.0;
  	  dice->clim.heat_mixed[i] = 0.0;
    	dice->clim.heat_interior[i] = 0.0;
	    dice->clim.heatflux_mixed[i] = 0.0;
 	   dice->clim.heatflux_interior[i] = 0.0;
  	}

  
  }

  // Calculate values for the exogenous variables -------------------------
  // dice->econ.partfract is used to determine the upper limit of the
  // control policy 'miu'
  calcExog(dice);
  for(int i=0; i<dice->config.nPeriods; i++) {
    if (i < 29) {
      dice->lims.miu_up[i] = 1.0;
    }
    else	{
      dice->lims.miu_up[i] = dice->dvars.limmiu;
    }
  }
  
  return;

}


//=====================================================
// populateDice
//
// Applies the parameter set to this instance
// of DICE
//=====================================================

void populateDICE(DICE *dice, vector<double> parms)
{
	// Parameter index variable
	int ipar = 0;

	// Preferences ---------------------------------
	dice->econ.elasmu = (parms[ipar++]);		// Elasticity of marginal utility of consumption
  dice->econ.prstp = (parms[ipar++]);		// Initial rate of social time preference (per year)
  
  // Population and Technology -------------------
  dice->econ.gama = (parms[ipar++]);		// Capital elasticity in production function
  dice->econ.pop0 = (parms[ipar++]);		// Initial world population [Millions]
  dice->econ.popadj = (parms[ipar++]);		// Growth rate to calibrate to 2050 population projection
  dice->econ.popasym = (parms[ipar++]);	// Asymptotic world population [Millions]
  dice->econ.dk = (parms[ipar++]);			// Depreciation rate on capital (per year)
  dice->econ.q0 = (parms[ipar++]);			// Initial world gross output [Trillions 2005 US $]
  dice->econ.k0 = (parms[ipar++]);			// Initial capital value [Trillions 2005 US $]
  dice->econ.a0 = (parms[ipar++]);			// Initial level of total factor productivity (TFP)
  dice->econ.ga0 = (parms[ipar++]);		// Initial growth rate for TFP (per 5 years)
  dice->econ.dela = (parms[ipar++]);		// Decline rate of TFP (per 5 years)
	dice->lims.optlrsav = (dice->econ.dk + 0.004) / (dice->econ.dk + 0.004 * dice->econ.elasmu + dice->econ.prstp) * dice->econ.gama;
	
	// Emissions parameters ------------------------
	dice->econ.gsigma1 = (parms[ipar++]);		// Initial growth of sigma (per year)
  dice->econ.dsig = (parms[ipar++]);		// Decline rate of decarbonization (per period)
  dice->econ.eland0 = (parms[ipar++]);		// Carbon emissions from land 2010 [GtCO2 per year]
  dice->econ.deland = (parms[ipar++]);		// Decline rate of land emissions (per period)
  dice->econ.e0 = (parms[ipar++]);		// Industrial emissions 2010 [GtCO2 per year]
  dice->dvars.miu0 = (parms[ipar++]);		// Initial emissions control rate for base case 2010
  //dice->econ.sig0 = dice->carb.e0 / (dice->econ.q0 * (1 - dice->econ.miu0));
						// Carbon Intensity 2010 [kgCO2 per output 2005 USD 2010]

	// Carbon cycle ---------------------------------
	// Initial conditions	
	dice->carb.mat0 = (parms[ipar++]);		// Initial concentration in atmosphere 2010 [GtC]
  dice->carb.mu0 = (parms[ipar++]);		// Initial concentration in upper strata [GtC]
  dice->carb.ml0 = (parms[ipar++]);		// Initial concentration in lower strata [GtC]
  dice->carb.mateq = (parms[ipar++]);		// Equilibrium concentration in atmosphere [GtC]
  dice->carb.mueq = (parms[ipar++]);		// Equilibrium concentration in upper strata [GtC]
  dice->carb.mleq = (parms[ipar++]);		// Equilibrium concentration in lower strata [GtC]
  
  // Flow parameters (Carbon cycle transition matricies)
  dice->carb.b12 = (parms[ipar++]);
  dice->carb.b23 = (parms[ipar++]);
  dice->carb.b11 = 1 - dice->carb.b12;
  dice->carb.b21 = dice->carb.b12 * dice->carb.mateq / dice->carb.mueq;
  dice->carb.b22 = 1 - dice->carb.b21 - dice->carb.b23;
  dice->carb.b32 = dice->carb.b23 * dice->carb.mueq / dice->carb.mleq;
  dice->carb.b33 = 1 - dice->carb.b32;

	// Climate model parameters ---------------------
	if(dice->ctrl.use_doeclim == 1) {
	  dice->clim.t2co = (parms[ipar++]);
  	apply_calibration(dice, dice->clim.t2co);
  }
  else {
  	dice->clim.t2xco2 = (parms[ipar++]);	// Equilibrium temperature impact [dC per doubling CO2]
  }
  dice->carb.fex0 = (parms[ipar++]);	// 2010 forcings of non-CO2 greenhouse gases (GHG) [Wm-2]
  dice->carb.fex1 = (parms[ipar++]);	// 2100 forcings of non-CO2 GHG [Wm-2]
	if(dice->ctrl.use_doeclim == 0) {
    dice->clim.tocean0 = (parms[ipar++]);	// Initial lower stratum temperature change [dC from 1900]
    dice->clim.tatm0 = (parms[ipar++]);	// Initial atmospheric temperature change [dC from 1900]
    dice->clim.c1 = (parms[ipar++]);	// Climate equation coefficient for upper level
    dice->clim.c3 = (parms[ipar++]);	// Transfer coefficient upper to lower stratum
    dice->clim.c4 = (parms[ipar++]);	// Transfer coefficient for lower level
  }
  dice->carb.fco22x = (parms[ipar++]);	// Forcings of equilibrium CO2 doubling [Wm-2]
  if(dice->ctrl.use_doeclim == 0) {
  	dice->clim.lam = dice->carb.fco22x / dice->clim.t2xco2;	// Climate model parameter
  }
  
  // Climate damage parameters ---------------------
  dice->econ.a10 = (parms[ipar++]);			// Initial damage intercept
  dice->econ.a20 = (parms[ipar++]);		// Initial damage quadratic term
  dice->econ.a1 = (parms[ipar++]);			// Damage intercept
  dice->econ.a2 = (parms[ipar++]);		// Damage quadratic term
  dice->econ.a3 = (parms[ipar++]);			// Damage exponent
  
  // Abatement cost --------------------------------
  dice->econ.expcost2 = (parms[ipar++]);		// Exponent of control cost function
  dice->econ.pback = (parms[ipar++]);		// Cost of backstop [2005$ per tCO2 2010]
  dice->econ.gback = (parms[ipar++]);		// Initial cost decline backstop [cost per period]
  dice->dvars.limmiu = (parms[ipar++]);		// Upper limit on control rate after 2150
  dice->econ.tnopol = (parms[ipar++]);		// Period before which no emissions controls base
  dice->econ.cprice0 = (parms[ipar++]);		// Initial base carbon price [2005$ per tCO2]
  dice->econ.gcprice = (parms[ipar++]);		// Growth rate of base carbon price (per year)
  
  // Availability of fossil fuels ------------------
  dice->carb.fosslim = (parms[ipar++]);	// Maximum cummulative extraction fossil fuels [GtC]
  dice->lims.cca_up = dice->carb.fosslim;
  
  
  // Scaling and inessential parameters ------------
	// "Note that these are unnecessary for the calculations but are for convenience"
	// Quoted directly from comments in original GAMS code
  dice->econ.scale1 = 0.0302455265681763;	// Multiplicitive scaling coefficient
  dice->econ.scale2 = -10993.704;     // Additive scaling coefficient
  
  // Objective
  dice->econ.utility = 0.0;		// Welfare function (Sum of discounted utility of per capita consumption)


	// Decision variables
 	for (int i=0; i<dice->config.nPeriods; i++) {
		if(i > 0) {
		  dice->dvars.miu[i] = parms[ipar++];
		} else {
			dice->dvars.miu[i] = dice->dvars.miu0;
		}
	}
	for (int i=0; i<dice->config.nPeriods; i++) {
		if (i < dice->config.nPeriods - 10) {
		  dice->dvars.s[i] = parms[ipar++];
		}	else {
			dice->dvars.s[i] = dice->lims.optlrsav;
		}
	}

	
	return;
}

//===========================================================================================
// allocateDICE(DICE *dice) 
//
// Defines and allocates the components of the DICE-DOECLIM model.
//
// Space is allocated for each model component. The hindcast forcing data are and
// calibration data are loaded. All the timeseries vectors are initialized.
//
// Return ()
// - 'dice' structure variables are updated
//===========================================================================================
void allocateDICE(DICE *dice)
{
	allocateCtrl(dice);
	allocateConfig(dice);
	allocateCarb(dice);
	allocateEcon(dice);
	allocateDvars(dice);
	allocateLims(dice);
	allocateClim(dice);
	
	if(dice->ctrl.use_doeclim == 1) {	
  	doeclim_load_calibration(dice, dice->config.n_calpoints);
  }

	return;
}

//=================================================================================
//
// void free_CDICE(DICE *dice)
//
// Frees up memory by destroying a DICE object
//
//=================================================================================
void free_CDICE(DICE *dice)
{
	// Delete the allocated arrays
	delete dice->config.dateSeries;	
	delete dice->carb.forcoth;	
  delete dice->carb.forc;		
  delete dice->carb.mat;		
  delete dice->carb.mu;		
  delete dice->carb.ml;		
  delete dice->econ.l;		
  delete dice->econ.al;		
  delete dice->econ.rr ;		
  delete dice->econ.ga;		
  delete dice->econ.gl;		
  delete dice->econ.gcost;		
  delete dice->econ.cost1;		
  delete dice->econ.partfract;	
  delete dice->econ.pbacktime;	
  delete dice->econ.scc;		
  delete dice->econ.cpricebase;	
  delete dice->econ.photel;	
  delete dice->econ.gsig;		
  delete dice->econ.sigma;		
  delete dice->econ.etree;		
  delete dice->econ.c;		
  delete dice->econ.k;		
  delete dice->econ.cpc;		
  delete dice->econ.i;		
  delete dice->econ.ri;		
  delete dice->econ.y;		
  delete dice->econ.ygross;	
  delete dice->econ.ynet;		
  delete dice->econ.damages;	
  delete dice->econ.damfrac;	
  delete dice->econ.abatecost;	
  delete dice->econ.mcabate;	
  delete dice->econ.periodu;	
  delete dice->econ.cprice;	
  delete dice->econ.cemutotper;	
  delete dice->econ.ri_disc;	
  delete dice->econ.pv_damages;	
  delete dice->econ.pv_abatecost;	
  delete dice->econ.totalcost;	
  delete dice->econ.pv_totalcost;	
  delete dice->econ.e;		
  delete dice->econ.eind;		
  delete dice->econ.cca;		
  delete dice->econ.maxcca_vec;
  delete dice->econ.utility_vec;
  delete dice->dvars.miu;	
  delete dice->dvars.s;	
  delete dice->lims.miu_up;	

	if(dice->ctrl.use_doeclim == 1) {	
		delete dice->clim.forc;
    delete dice->clim.IB;
    delete dice->clim.A;
    delete dice->clim.Ker;
    delete dice->clim.temp;
    delete dice->clim.temp_landair;
    delete dice->clim.temp_sst;
    delete dice->clim.heat_mixed;
    delete dice->clim.heat_interior;
    delete dice->clim.heatflux_mixed;
    delete dice->clim.heatflux_interior;
    delete dice->clim.cs_cal;
    delete dice->clim.kv_cal;
    delete dice->clim.alpha_cal;
  } else {
    delete dice->clim.tatm;		
    delete dice->clim.tocean;	
  }
  
/*  // Delete the structures in the DICE object
  delete dice->config;
  delete dice->dvars;
  delete dice->ctrl;
  delete dice->clim;
  delete dice->econ;
  delete dice->carb;
  delete dice->lims;
*/  
  // Delete the parent structure of the DICE object
//  delete dice;
  
}

//#endif  //CDICEMOD
