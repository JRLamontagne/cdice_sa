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

extern DICE dice;
extern DICE* dicePtr;


//============================================================================
// calcExog(DICE *dice)
// 
// Calculate exogenous time series in the Config section of the dice structure
//
// Return ()
// - 'dice' structure config time series are updated
//============================================================================
void calcExog(DICE *dice)
{
  Config *config = &dice->config;
  Econ *econ = &dice->econ;
  Carb *carb = &dice->carb;
  Clim *clim = &dice->clim;
  Dvars *dvars = &dice->dvars;
  Ctrl *ctrl = &dice->ctrl;
  int tlim = config->nPeriods;
  
  
  // Initial timestep calculations and settings ----------------------------
  econ->l[0] = econ->pop0;			// Population
  econ->ga[0] = econ->ga0;			// Growth rate for TFP
  econ->al[0] = econ->a0;			// Level of TFP
  econ->gsig[0] = econ->gsigma1;		// Change in sigma (cumulative improvement in energy efficiency)
  econ->sigma[0] = econ->sig0;		// CO2-equivalent-emissions output ratio
  econ->pbacktime[0] = econ->pback;	// Cost of backstop
  econ->cost1[0] = econ->pbacktime[0]*econ->sigma[0]/econ->expcost2/1000;

  // Adjusted cost for backstop
  econ->etree[0] = econ->eland0;		// Emissions from deforestation
  econ->cumetree[0] = 100.0;
  econ->rr[0] = 1;					// Average utility social discount rate
  carb->forcoth[0] = carb->fex0;		// Forcing from other GHG

  // Fraction of emissions under control regime
  econ->cpricebase[0] = econ->cprice0; // Carbon price in base case
  econ->sig0 = econ->e0 / (econ->q0 * (1 - dvars->miu0));
                                       // Carbon Intensity 2010 [kgCO2 per output 2005 USD 2010]
  
  // Iterate over the full time series
  for (int t=1; t<tlim; t++)
    {
      // Population
      econ->l[t] = econ->l[t-1]*(pow(econ->popasym/econ->l[t-1], econ->popadj));
      
      // Growth rate for TFP
      econ->ga[t] = econ->ga0*exp(-1*econ->dela*config->tstep*double(t));
      
      // Level for TFP
      econ->al[t] = econ->al[t-1]/(1-econ->ga[t-1]);
      
      // Change in sigma (cumulative improvement in energy efficiency)
      econ->gsig[t] = econ->gsig[t-1]*(pow(1+econ->dsig, config->tstep));
      
      // CO2-equivalent-emissions output ratio
      econ->sigma[t] = econ->sigma[t-1]*exp(econ->gsig[t-1]*config->tstep);
      
      // Cost of backstop
      econ->pbacktime[t] = econ->pback*pow(1-econ->gback, double(t));
      
      // Adjusted cost for backstop
      econ->cost1[t] = econ->pbacktime[t]*econ->sigma[t]/econ->expcost2/1000;
      
      // Emissions from deforestation
      econ->etree[t] = econ->eland0*pow(1-econ->deland, double(t));
      
      // Cumulative emissions from land
      econ->cumetree[t] = econ->cumetree[t-1]+econ->etree[t-1]*(5/3.666);
      
      // Average utility social discount rate	
      econ->rr[t] = 1/pow(1+econ->prstp, config->tstep*double(t));
      
      // Forcings from other GHG
      if (t < 18) {
	carb->forcoth[t] = carb->fex0+(1.0/17.0)*(carb->fex1-carb->fex0)*double(t);
      }
      else {
	carb->forcoth[t] = carb->fex0+(carb->fex1-carb->fex0);
      }
      
      // Carbon price in base case
      econ->cpricebase[t] = econ->cprice0*pow(1+econ->gcprice, config->tstep*double(t));
      
    }
  
  return;
}

//=============================================================================
// calcCarb1(DICE *dice)
//
// Calculate the first time step of the carbon-related variables
//
// Return ()
// - 'dice' carb structure variables are updated
//=============================================================================
void calcCarb1(DICE *dice)
{

  Carb *carb = &dice->carb;
  Lims *lims = &dice->lims;
  
  // Carbon pools
  carb->mat[0] = carb->mat0;
  if(carb->mat[0] < lims->mat_lo) {
    carb->mat[0] = lims->mat_lo;
  }
  
  carb->mu[0] = carb->mu0;
  if(carb->mu[0] < lims->mu_lo) {
    carb->mu[0] = lims->mu_lo;
  }

  carb->ml[0] = carb->ml0;
  if(carb->ml[0] < lims->ml_lo) {
    carb->ml[0] = lims->ml_lo;
  }

  // Radiative forcing
  carb->forc[0] = carb->fco22x*(log(carb->mat[0]/588.000)/log(2.0))+carb->forcoth[0];
  
  return;
}


//=============================================================================
// calcClim1(DICE *dice)
//
// Calculate the first time step of the climate-related variables
//
// Return ()
// - 'dice' clim structure variables are updated
//=============================================================================
void calcClim1(DICE *dice)
{
  Clim *clim = &dice->clim;
  Lims *lims = &dice->lims;

  // Atmospheric temperature
  clim->tatm[0] = clim->tatm0;
  if(clim->tatm[0] < lims->tatm_lo) {
    clim->tatm[0] = lims->tatm_lo;
  }
  if(clim->tatm[0] > lims->tatm_up) {
    clim->tatm[0] = lims->tatm_up;
  }
  
  // Oceanic temperature
  clim->tocean[0] = clim->tocean0;
  if(clim->tocean[0] < lims->tocean_lo) {
    clim->tocean[0] = lims->tocean_lo;
  }
  if(clim->tocean[0] > lims->tocean_up) {
    clim->tocean[0] = lims->tocean_up;
  }
  
  return;
}

//=============================================================================
//	calcEcon1(DICE *dice)
//
// Calculate the first time step for the economic variables
//
// Return ()
// - 'dice' endog structure variables updated
//=============================================================================
void calcEcon1(DICE *dice)
{

  Econ *econ  = &dice->econ;
  Clim *clim = &dice->clim;
  Dvars *dvars = &dice->dvars;
  Ctrl *ctrl = &dice->ctrl;
  Lims *lims = &dice->lims;
  Config *config = &dice->config;

  // Variable to hold the appropriate temperature
  double temp = 0.0;
  
  // Where do we get the temperature from?
  if(ctrl->use_doeclim == 1) {
    temp = clim->temp[config->hind_ns-1] - clim->temp[config->i1900];
    //temp = clim->temp[config->hind_ns-1];
  }
  else {
    temp = clim->tatm[0];
  }
  
  // Marginal cost of abatement
  econ->mcabate[0] = econ->pbacktime[0]*pow(dvars->miu[0],(econ->expcost2-1));
  
  // Damages fraction
  econ->damfrac[0] = econ->a1*temp+econ->a2*pow(temp,econ->a3);
  
  // Capital stock
  econ->k[0] = econ->k0;
  if(econ->k[0] < lims->k_lo) {
    econ->k[0] = lims->k_lo;
  }
  
  // Carbon price
  econ->cprice[0] = econ->pbacktime[0]*pow(dvars->miu[0],(econ->expcost2-1));
  
  // Gross world product (gross abatement and damages)
  econ->ygross[0] = econ->al[0]*pow((econ->l[0]/1000.0),(1.0-econ->gama))*pow(econ->k[0],econ->gama);
  if(econ->ygross[0] < lims->ygross_lo) {
    econ->ygross[0] = lims->ygross_lo;
  }
  
  // Industrial emissions
  econ->eind[0] = econ->sigma[0]*econ->ygross[0]*(1.0-dvars->miu[0]);
  
  // Total emissions	
  econ->e[0] = econ->eind[0]+econ->etree[0];
  
  // Cumulative emissions
  econ->cca[0] = 90.0;
  if(econ->cca[0] > lims->cca_up) {
    econ->cca[0] = lims->cca_up + 1.0;
  }
  
  // Cumulative total emissions
  econ->ccatot[0] = econ->cca[0] + econ->cumetree[0];

  // Net output damages
  econ->ynet[0] = econ->ygross[0]*(1.0-econ->damfrac[0]);
  
  // Abatement costs
  econ->abatecost[0] = econ->ygross[0]*econ->cost1[0]*pow(dvars->miu[0],econ->expcost2);
  
  // Damages
  econ->damages[0] = econ->ygross[0]*econ->damfrac[0];
  
  // Gross world product (net of abatement and damages)
  econ->y[0] = econ->ynet[0]-econ->abatecost[0];
  if(econ->y[0] < lims->y_lo) {
    econ->y[0] = lims->y_lo;
  }
  
  // Investment
  econ->i[0] = dvars->s[0]*econ->y[0];
  if(econ->i[0] < lims->i_lo) {
    econ->i[0] = lims->i_lo;
  }
  
  // Consumption
  econ->c[0] = econ->y[0]-econ->i[0];
  if(econ->c[0] < lims->c_lo) {
    econ->c[0] = lims->c_lo;
  }

  // Per capita consumption
  econ->cpc[0] = 1000.0*econ->c[0]/econ->l[0];
  if(econ->cpc[0] < lims->cpc_lo) {
    econ->cpc[0] = lims->cpc_lo;
  }

  return;
}


//=============================================================================
// calcCarb(DICE *dice, int t)
//
// Calculate the carbon model component for time step 't'
//
// Return()
// - 'dice' endogenous structure variables updated
//=============================================================================
void calcCarb(DICE *dice, int t)
{

  Econ *econ = &dice->econ;
  Carb *carb = &dice->carb;
  Lims *lims = &dice->lims;
  
  // Carbon pools (carbon cycle box model)
  carb->mat[t] = (econ->e[t-1]*(5.0/3.666))+carb->b11*carb->mat[t-1]+carb->b21*carb->mu[t-1];
  if(carb->mat[t] < lims->mat_lo) {
    carb->mat[t] = lims->mat_lo;
  }

  carb->mu[t] = carb->b12*carb->mat[t-1]+carb->b22*carb->mu[t-1]+carb->b32*carb->ml[t-1];
  if(carb->mu[t] < lims->mu_lo) {
    carb->mu[t] = lims->mu_lo;
  }
  carb->ml[t] = carb->b33*carb->ml[t-1]+carb->b23*carb->mu[t-1];
  if(carb->ml[t] < lims->ml_lo) {
    carb->ml[t] = lims->ml_lo;
  }
  
  // Radiatice forcing
  carb->forc[t] = carb->fco22x*(log(carb->mat[t]/588.000)/log(2.0))+carb->forcoth[t];
  
  return;
}


//=============================================================================
// calcClim(DICE *dice, int t)
//
// Calculate the climate model component for time step 't'
//
// Return()
// - 'dice' endogenous structure variables updated
//=============================================================================
void calcClim(DICE *dice, int t)
{

  Config *config = &dice->config;
  Econ *econ = &dice->econ;
  Carb *carb = &dice->carb;
  Clim *clim = &dice->clim;
  Dvars *dvars = &dice->dvars;
  Ctrl *ctrl = &dice->ctrl;
  Lims *lims = &dice->lims;

  clim->tatm[t] = clim->tatm[t-1]+clim->c1*((carb->forc[t]-((carb->fco22x/clim->t2xco2)*clim->tatm[t-1]))-(clim->c3*(clim->tatm[t-1]-clim->tocean[t-1])));
  if(clim->tatm[t] < lims->tatm_lo) {
    clim->tatm[t] = lims->tatm_lo;
  }
  if(clim->tatm[t] > lims->tatm_up) {
    clim->tatm[t] = lims->tatm_up;
  }

  clim->tocean[t] = clim->tocean[t-1]+clim->c4*(clim->tatm[t-1]-clim->tocean[t-1]);
  if(clim->tocean[t] < lims->tocean_lo) {
    clim->tocean[t] = lims->tocean_lo;
  }
  if(clim->tocean[t] > lims->tocean_up) {
    clim->tocean[t] = lims->tocean_up;
  }

  return;
}

//=============================================================================
// calcEcon(DICE *dice, int t)
//
// Calculate the economic model component for the time step 't'.
//
// Return()
// - 'dice' endogenous structure variables updated
//=============================================================================
void calcEcon(DICE *dice, int t)
{

  Config *config = &dice->config;
  Econ *econ = &dice->econ;
  Carb *carb = &dice->carb;
  Clim *clim = &dice->clim;
  Dvars *dvars = &dice->dvars;
  Ctrl *ctrl = &dice->ctrl;
  Lims *lims = &dice->lims;
	
  // Variable to hold the appropriate temperature
  double temp = 0.0;
  
  // Where do we get the temperature from?
  if(ctrl->use_doeclim == 1) {
    int this_index = config->hind_ns+(t*config->tstep)-1;
    temp = clim->temp[this_index] - clim->temp[config->i1900];
    //temp = clim->temp[this_index];
  }
  else {
    temp = clim->tatm[t];
  }
  
  // Apply the limit test to this temperature (for numeric
  // stability in the damage calculation
  if(temp < lims->tatm_lo) {
    temp = lims->tatm_lo;
  }
  if(temp > lims->tatm_up) {
    temp = lims->tatm_up;
  }  

  // Marginal cost of abatement
  econ->mcabate[t] = econ->pbacktime[t]*pow(dvars->miu[t],(econ->expcost2-1));
  
  // Damages fraction
  econ->damfrac[t] = econ->a1*temp+econ->a2*pow(temp,econ->a3);
  
  // Capital stock
  econ->k[t] = pow((1-econ->dk),config->tstep)*econ->k[t-1]+config->tstep*econ->i[t-1];
  if(econ->k[t] < lims->k_lo) {
    econ->k[t] = lims->k_lo;
  }
  
  // Carbon price
  econ->cprice[t] = econ->pbacktime[t]*pow(dvars->miu[t],(econ->expcost2-1));
  
  // Gross world product (gross abatement and damages)
  econ->ygross[t] = econ->al[t]*pow((econ->l[t]/1000.0),(1.0-econ->gama))*pow(econ->k[t],econ->gama);
  if(econ->ygross[t] < lims->ygross_lo) {
    econ->ygross[t] = lims->ygross_lo;
  }
  
  // Net output damages
  econ->ynet[t] = econ->ygross[t]*(1.0-econ->damfrac[t]);
  
  // Abatement costs
  econ->abatecost[t] = econ->ygross[t]*econ->cost1[t]*pow(dvars->miu[t],econ->expcost2);
  
  // Damages
  econ->damages[t] = econ->ygross[t]*econ->damfrac[t];
  
  // Industrial emissions
  econ->eind[t] = econ->sigma[t]*econ->ygross[t]*(1.0-dvars->miu[t]);
  
  // Total emissions	
  econ->e[t] = econ->eind[t]+econ->etree[t];
  
  // Cumulative emissions
  econ->cca[t] = econ->cca[t-1]+econ->eind[t-1]*5.0/3.666;
  if(econ->cca[t] > lims->cca_up) {
    econ->cca[t] = lims->cca_up + 1.0;
    if(lims->cca_up_step == 0) {
    	lims->cca_up_step = t;
    }
  }
  
  // Cumulative total emissions
  econ->ccatot[t] = econ->cca[t] + econ->cumetree[t];

  // Gross world product (net of abatement and damages)
  econ->y[t] = econ->ynet[t]-econ->abatecost[t];
  if(econ->y[t] < lims->y_lo) {
    econ->y[t] = lims->y_lo;
  }
  
  // Investment
  econ->i[t] = dvars->s[t]*econ->y[t];
  if(econ->i[t] < lims->i_lo) {
    econ->i[t] = lims->i_lo;
  }

  // Consumption
  econ->c[t] = econ->y[t]-econ->i[t];
  if(econ->c[t] < lims->c_lo) {
    econ->c[t] = lims->c_lo;
  }
  // Per capita consumption
  econ->cpc[t] = 1000.0*econ->c[t]/econ->l[t];
  if(econ->cpc[t] < lims->cpc_lo) {
    econ->cpc[t] = lims->cpc_lo;
  }

  // Real interest rate
  econ->ri[t-1] = (1.0+econ->prstp)*pow((econ->cpc[t]/econ->cpc[t-1]),(econ->elasmu/config->tstep))-1.0;
  
  return;
}


//=============================================================================
// postProcess(DICE *dice)
//
// Calculate objective function and any additional quantities added by the user
//
// Return()
// - 'dice' endogenous structure variables are updated
//=============================================================================
void postProcess(DICE *dice)
{
  Config *config = &dice->config;
  Econ *econ = &dice->econ;
  Carb *carb = &dice->carb;
  Clim *clim = &dice->clim;
  Dvars *dvars = &dice->dvars;
  Ctrl *ctrl = &dice->ctrl;


  // Calculate additional model outputs ------------------------------------------
  for (int t=0; t<config->nPeriods; t++)
    {
      
      // One period utility function
      econ->periodu[t] = (pow((econ->c[t]*1000.0/econ->l[t]),(1.0-econ->elasmu))-1.0)/(1.0-econ->elasmu)-1.0;
      
      // Period utility
      econ->cemutotper[t] = econ->periodu[t]*econ->l[t]*econ->rr[t];
      
      // Total cost of abatement and damages
      econ->totalcost[t] = econ->damages[t] + econ->abatecost[t];
      
      // Atmospheric fraction 1850
		carb->atfrac[t] = (carb->mat[t]-588)/(econ->ccatot[t]+0.000001);
		
		// Atmospheric fraction 2010
		carb->atfrac2010[t] = (carb->mat[t] - carb->mat0)/(0.00001+econ->ccatot[t]-econ->ccatot[0]);
		
		// Atmospheric concentration (mixing ratio)
		carb->ppm[t] = carb->mat[t]/2.13;
      
    }
  
  
  // Discounted utility of consumption (objective function) -----------------------
  double tempsum;
  tempsum = 0.0;
  
  for (int t=0; t<config->nPeriods; t++)
    {
      tempsum += econ->cemutotper[t];
    }
  
  econ->utility = (config->tstep*econ->scale1*tempsum)+econ->scale2;
  
  
  // Present Value Calculations for damages and abatement costs --------------------
  econ->pv_damages[0] = econ->damages[0];
  econ->pv_abatecost[0] = econ->abatecost[0];
  econ->pv_totalcost[0] = econ->totalcost[0];
  
  for (int t=1; t<config->nPeriods; t++)
    {
      if(econ->ri[t] <= 0.0) {
	econ->ri_disc[t] = 0.0;
      }
      else {
	econ->ri_disc[t] = 1.0/pow(1.0 + econ->ri[t],config->tstep*double(t));
      }
      
      econ->pv_damages[t] = econ->pv_damages[t-1]+econ->ri_disc[t]*econ->damages[t];
      econ->pv_abatecost[t] = econ->pv_abatecost[t-1]+econ->ri_disc[t]*econ->abatecost[t];
      econ->pv_totalcost[t] = econ->pv_totalcost[t-1]+econ->ri_disc[t]*econ->totalcost[t];
    }
  
  return;
}


//=============================================================================
// processModel(DICE *dice)
//
// Perform a full iteration of the model
// - calculate the exogenous time series
// - do time step 1 of both climate and econ model components
// - loop the remaining time steps (climate and econ)
// - calculate the objective function and other model results
//
// Return()
// - Full 'dice' structure is calculated.
//=============================================================================
void processModel(DICE *dice)
{
	
  // Calculate exogenous time series
  calcExog(dice);
  
   // Calculate the first time step of carbon model
  calcCarb1(dice);

  // Calculate the first time step of climate model
  if(dice->ctrl.use_doeclim == 1) {
    doeclim_hindcast(dice);
  }
  else {
    calcClim1(dice);
  }
  
  // Calculate the first time step of economic model
  calcEcon1(dice);
  
  // Calculate the rest of the time series
  for (int tIndex=1; tIndex<dice->config.nPeriods; tIndex++)
    {
      calcCarb(dice, tIndex);

      if(dice->ctrl.use_doeclim == 1) {
				doeclim_dice_ts(dice, tIndex);
    	}
      else {
				calcClim(dice, tIndex);
      }

      calcEcon(dice, tIndex);
    }
  
  // Do the post-processing
  postProcess(dice);
  
  return;
}


//#endif    //CDICEMOD

