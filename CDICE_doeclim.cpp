/* Copyright (C) 2013 Gregory Garner, Klaus Keller, 
   Patrick Reed, Martha Butler, and others

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

#include "CDICE.h"

extern DICE dice;
extern DICE* dicePtr;

//=========================================================================
//
// void count_lines(const char *filename, int *numlines)
//
// Counts the number of lines in a file
// 
// filename = name of the files
//
// numlines = number of lines in the file
//
//=========================================================================
void count_lines(const char *filename, int *numlines)
{
	int ch;	
	*numlines = 0;
	FILE * fp = fopen(filename, "r");
	
	do {
		ch = fgetc(fp);
		if(ch == '\n') { *numlines = *numlines + 1; }
	} while (ch != EOF);

//	if(ch != '\n' && *numlines != 0) 
//    *numlines = *numlines + 1;

	fclose(fp);
	
	return;
}


//=========================================================================
//
// void invert_1d_2x2_matrix(double * x, double * y)
//
// Calculates the inverse of 'x' and stores it in 'y'
//
// Assumes that the matrix is set up as such:
//
//      x = [a,b,c,d]  -->  x = |a, b|
//                              |c, d|
// 
// x = Input 1-d matrix
//
// y = Inverted 1-d matrix
//
//=========================================================================

void invert_1d_2x2_matrix(double * x, double * y)
{
  double temp_d = (x[0]*x[3] - x[1]*x[2]);
  if(temp_d == 0) {
    fprintf(stderr, "Matrix inversion divide by zero. Aborting...\n");
    exit(EXIT_FAILURE);
  }
  double temp = 1/temp_d;
  y[0] = temp * x[3];
  y[1] = temp * -1 * x[1];
  y[2] = temp * -1 * x[2];
  y[3] = temp * x[0];

  return;
}

//=========================================================================
//
// void sum_1d_2x2_matrix(double * x, double * y, double * z)
//
// Calculates the element-wise sum of matrix 'x' and 'y'
// storing output in 'z'.
//
// Assumes that the matricies are set up as such:
//
//      x = [a,b,c,d]  -->  x = |a, b|
//                              |c, d|
// 
// x = Input 1-d matrix
// y = Input 1-d matrix
//
// z = Element-wise sum of matricies 'x' and 'y'
//
//=========================================================================

void sum_1d_2x2_matrix(double * x, double * y, double * z)
{
  for(int i=0; i < 4; i++) {
    z[i] = x[i] + y[i];
  }

  return;
}

//=========================================================================
//
// doeclim_ts(int tstep, int dt, int ns, double forcing, 
//            double flnd, double powtoheat, double cal, double cas,
//            double taucfl, double taukls, double taucfs, double tauksl,
//            double bsi, double *IB, double *A, double fso,
//            double taudif, double *Ker,
//            double *temp, double *temp_landair, double *temp_sst, 
//            double *heat_mixed, double *heat_interior, 
//            double *heatflux_mixed, double *heatflux_interior)
//
// Calculates an iteration of the DOEclim model for the provided tstep.
//
// 
// tstep = Time step (year of interest; 0 = DOEclim base year)
// dt = Number of years in a single time step
// ns = Number of time steps in DOEclim
// forcing = Radiative forcing (W/m2)
// flnd, powtoheat, cal, cas, taucfl, taukls, taucfs, tauksl,
// bsi, IB, A, fso, taudif, Ker = Model parameters
// 
// temp:               global mean temperature anomaly (K), preindustrial
// temp_landair:       land air temperature anomaly (K)
// temp_sst:           sea surface temperature anomaly (K)
// heat_mixed:         mixed layer heat anomaly (10^22 J)
// heat_interior:      interior ocean heat anomaly (10^22 J)
// heatflux_mixed:     heat uptake of the mixed layer (W/m^2)
// heatflux_interior:  heat uptake of the interior ocean (W/m^2)
//
// Assumptions: land surface temperature = land air temperature
//              mixed layer temperature  = sea surface temperatures 
//                                       = marine air temperature / bsi 
//
//=========================================================================

void doeclim_ts(int tstep, int dt, int ns, double *forcing, 
            double flnd, double powtoheat, double cal, double cas,
            double taucfl, double taukls, double taucfs, double tauksl,
            double bsi, double *IB, double *A, double fso,
            double taudif, double *Ker,
            double *temp, double *temp_landair, double *temp_sst, 
            double *heat_mixed, double *heat_interior, 
            double *heatflux_mixed, double *heatflux_interior)
{
  // Initialize variables for time-stepping through the model
  double DQ1 = 0.0;
  double DQ2 = 0.0;
  double QC1 = 0.0;
  double QC2 = 0.0;
  double DelQL = 0.0;
  double DelQO = 0.0;
  double DPAST1 = 0.0;
  double DPAST2 = 0.0;
  double DTEAUX1 = 0.0;
  double DTEAUX2 = 0.0;

	// Reset the endogenous varibales for this time step
	temp[tstep] = 0.0;
	temp_landair[tstep] = 0.0;
	temp_sst[tstep] = 0.0;
	heat_mixed[tstep] = 0.0;
	heat_interior[tstep] = 0.0;
	heatflux_mixed[tstep] = 0.0;
	heatflux_interior[tstep] = 0.0;

  // Assume land and ocean forcings are equal to global forcing
  double *QL = forcing;
  double *QO = forcing;
  
  if (tstep > 0) {
    
    DelQL = QL[tstep] - QL[tstep - 1];
    DelQO = QO[tstep] - QO[tstep - 1];
    
    // Assume linear forcing change between tstep and tstep+1
    QC1 = (DelQL/cal*(1.0/taucfl+1.0/taukls)-bsi*DelQO/cas/taukls);
    QC2 = (DelQO/cas*(1.0/taucfs+bsi/tauksl)-DelQL/cal/tauksl);
    QC1 = QC1 * pow(double(dt), 2.0)/12.0;
    QC2 = QC2 * pow(double(dt), 2.0)/12.0;

    // ----------------- Initial Conditions --------------------
    // Initialization of temperature and forcing vector:
    // Factor 1/2 in front of Q in Equation A.27, EK05, and Equation 2.3.27, TK07 is a typo! 
    // Assumption: linear forcing change between n and n+1
    DQ1 = 0.5*double(dt)/cal*(QL[tstep]+QL[tstep-1]);
    DQ2 = 0.5*double(dt)/cas*(QO[tstep]+QO[tstep-1]);
    DQ1 = DQ1 + QC1;
    DQ2 = DQ2 + QC2;

    // ---------- SOLVE MODEL ------------------
    // Calculate temperatures
    for(int i = 0; i <= tstep; i++) {
      DPAST2 = DPAST2 + temp_sst[i] * Ker[ns-tstep+i-1];
    }
    DPAST2 = DPAST2 * fso * pow((double(dt)/taudif), 0.5);
    
    DTEAUX1 = A[0] * temp_landair[tstep-1] + A[1] * temp_sst[tstep-1];
    DTEAUX2 = A[2] * temp_landair[tstep-1] + A[3] * temp_sst[tstep-1];

    temp_landair[tstep] = IB[0] * (DQ1 + DPAST1 + DTEAUX1) + IB[1] * (DQ2 + DPAST2 + DTEAUX2);
    temp_sst[tstep] = IB[2] * (DQ1 + DPAST1 + DTEAUX1) + IB[3] * (DQ2 + DPAST2 + DTEAUX2);    
    
  }
  
  else {  // Handle the initial conditions
    
    temp_landair[0] = 0.0;
    temp_sst[0] = 0.0;

  }

  temp[tstep] = flnd * temp_landair[tstep] + (1.0 - flnd) * bsi * temp_sst[tstep];

  // Calculate ocean heat uptake [W/m^2]
  // heatflux[tstep] captures in the heat flux in the period between tstep-1 and tstep.
  // Numerical implementation of Equation 2.7, EK05, or Equation 2.3.13, TK07)
  // ------------------------------------------------------------------------
  
  if (tstep > 0) {
    
    heatflux_mixed[tstep] = cas*(temp_sst[tstep] - temp_sst[tstep-1]);
    
    for (int i=0; i < tstep; i++) { 
      heatflux_interior[tstep] = heatflux_interior[tstep] + temp_sst[i]*Ker[ns-tstep+i];
    }
    
    heatflux_interior[tstep] = cas*fso/pow((taudif*dt), 0.5)*(2.0*temp_sst[tstep] - heatflux_interior[tstep]);
    
    heat_mixed[tstep] = heat_mixed[tstep-1] + heatflux_mixed[tstep] * (powtoheat*dt);
    
    heat_interior[tstep] = heat_interior[tstep-1] + heatflux_interior[tstep] * (fso*powtoheat*dt);
    
  }
  
  else {   // Handle the initial conditions
    
    heatflux_mixed[0] = 0.0;
    heatflux_interior[0] = 0.0;
    heat_mixed[0] = 0.0;
    heat_interior[0] = 0.0;
    
  }


  return;
}


//=========================================================================
//
// void doeclim_init(double t2co, double kappa, int dt, int ns, 
//		  double *flnd, double *powtoheat, double *cal, double *cas, double *taucfl,
//		  double *taukls, double *taucfs, double *tauksl, double *bsi, 
//		  double *IB, double *A, double *fso, double *taudif, double *Ker)
//
// Initializes the doeclim climate model parameters.
//
// Input parameters:
//   t2co = Climate sensitivity (K)
//   kappa = Ocean heat diffusivity (cm^2 s^-1)
//   dt = Number of years per time step
//   ns = Number of time steps
//
//=========================================================================

void doeclim_init(double t2co, double kappa, int dt, int ns, 
		  double *flnd, double *powtoheat, double *cal, double *cas, double *taucfl,
		  double *taukls, double *taucfs, double *tauksl, double *bsi, 
		  double *IB, double *A, double *fso, double *taudif, double *Ker)
{
  
  // Define and initialize parameters
  double pi = 3.141592653589793;
  double B[4] = {0.0, 0.0, 0.0, 0.0};
  double C[4] = {0.0, 0.0, 0.0, 0.0};
  double KT0[ns];
  double KTA1[ns];
  double KTB1[ns];
  double KTA2[ns];
  double KTB2[ns];
  double KTA3[ns];
  double KTB3[ns];

  for(int i=0; i<ns; i++) {
    KT0[i] = 0.0;
    KTA1[i] = 0.0;
    KTB1[i] = 0.0;
    KTA2[i] = 0.0;
    KTB2[i] = 0.0;
    KTA3[i] = 0.0;
    KTB3[i] = 0.0;
  }

  // Define the doeclim parameters
  double ak = 0.31;
  double bk = 1.59;
  double csw = 0.13;
  double earth_area = 5100656 * pow(10.0, 8);
  double kcon = 3155.0;
  double q2co = 3.7; 
  double rlam = 1.43;
  double secs_per_Year = 31556926.0;
  double zbot = 4000.0;
  *bsi = 1.3;
  *cal = 0.52;
  *cas = 7.80;
  *flnd = 0.29;
  *fso = 0.95;
  
  // Dependent Model Parameters
  double ocean_area = (1.0 - *flnd) * earth_area;
  double cnum = rlam * *flnd + *bsi * (1.0 - *flnd);
  double cden = rlam * *flnd - ak * (rlam - *bsi);
  double cfl = *flnd * cnum / cden * q2co / t2co - bk * (rlam - *bsi) / cden;
  double cfs = (rlam * *flnd - ak / (1.0 - *flnd) * (rlam - *bsi)) * cnum / cden * q2co / t2co + rlam * *flnd / (1.0 - *flnd) * bk * (rlam - *bsi) / cden;
  double kls = bk * rlam * *flnd / cden - ak * *flnd * cnum / cden * q2co / t2co; 
  double keff = kcon * kappa;
  double taubot = pow(zbot,2) / keff;
  *powtoheat = ocean_area * secs_per_Year / pow(10.0,22);
  *taucfs = *cas / cfs;
  *taucfl = *cal / cfl;
  *taudif = pow(*cas,2) / pow(csw,2) * pi / keff;
  *tauksl  = (1.0 - *flnd) * *cas / kls;
  *taukls  = *flnd * *cal / kls;
  
  // First order
  KT0[ns-1] = 4.0 - 2.0 * pow(2.0, 0.5);
  KTA1[ns-1] = -8.0 * exp(-taubot / double(dt)) + 4.0 * pow(2.0, 0.5) * exp(-0.5 * taubot / double(dt));
  KTB1[ns-1] = 4.0 * pow((pi * taubot / double(dt)), 0.5) * (1.0 + erf(pow(0.5 * taubot / double(dt), 0.5)) - 2.0 * erf(pow(taubot / double(dt), 0.5)));

  // Second order
  KTA2[ns-1] =  8.0 * exp(-4.0 * taubot / double(dt)) - 4.0 * pow(2.0, 0.5) * exp(-2.0 * taubot / double(dt));
  KTB2[ns-1] = -8.0 * pow((pi * taubot / double(dt)), 0.5) * (1.0 + erf(pow((2.0 * taubot / double(dt)), 0.5)) - 2.0 * erf(2.0 * pow((taubot / double(dt)), 0.5)) );

  // Third order
  KTA3[ns-1] = -8.0 * exp(-9.0 * taubot / double(dt)) + 4.0 * pow(2.0, 0.5) * exp(-4.5 * taubot / double(dt));
  KTB3[ns-1] = 12.0 * pow((pi * taubot / double(dt)), 0.5) * (1.0 + erf(pow((4.5 * taubot / double(dt)), 0.5)) - 2.0 * erf(3.0 * pow((taubot / double(dt)), 0.5)) );

  // Calculate the kernel component vectors
  for(int i=0; i<(ns-1); i++) {
    
    // First order
    KT0[i] = 4.0 * pow((double(ns-i)), 0.5) - 2.0 * pow((double(ns+1-i)), 0.5) - 2.0 * pow(double(ns-1-i), 0.5);      
    KTA1[i] = -8.0 * pow(double(ns-i), 0.5) * exp(-taubot / double(dt) / double(ns-i)) + 4.0 * pow(double(ns+1-i), 0.5) * exp(-taubot / double(dt) / double(ns+1-i)) + 4.0 * pow(double(ns-1-i), 0.5) * exp(-taubot/double(dt) / double(ns-1-i));
    KTB1[i] =  4.0 * pow((pi * taubot / double(dt)), 0.5) * ( erf(pow((taubot / double(dt) / double(ns-1-i)), 0.5)) + erf(pow((taubot / double(dt) / double(ns+1-i)), 0.5)) - 2.0 * erf(pow((taubot / double(dt) / double(ns-i)), 0.5)) );

    // Second order
    KTA2[i] =  8.0 * pow(double(ns-i), 0.5) * exp(-4.0 * taubot / double(dt) / double(ns-i)) - 4.0 * pow(double(ns+1-i), 0.5) * exp(-4.0 * taubot / double(dt) / double(ns+1-i)) - 4.0 * pow(double(ns-1-i), 0.5) * exp(-4.0 * taubot / double(dt) / double(ns-1-i));
    KTB2[i] = -8.0 * pow((pi * taubot / double(dt)), 0.5) * ( erf(2.0 * pow((taubot / double(dt) / double(ns-1-i)), 0.5)) + erf(2.0 * pow((taubot / double(dt) / double(ns+1-i)), 0.5)) - 2.0 * erf(2.0 * pow((taubot / double(dt) / double(ns-i)), 0.5)) );

    // Third order
    KTA3[i] = -8.0 * pow(double(ns-i), 0.5) * exp(-9.0 * taubot / double(dt) / double(ns-i)) + 4.0 * pow(double(ns+1-i), 0.5) * exp(-9.0 * taubot / double(dt) / double(ns+1-i)) + 4.0 * pow(double(ns-1-i), 0.5) * exp(-9.0 * taubot / double(dt) / double(ns-1-i)); 
    KTB3[i] = 12.0 * pow((pi * taubot / double(dt)), 0.5) * ( erf(3.0 * pow((taubot / double(dt) / double(ns-1-i)), 0.5)) + erf(3.0 * pow((taubot / double(dt) / double(ns+1-i)), 0.5)) - 2.0 * erf(3.0 * pow((taubot / double(dt) / double(ns-i)), 0.5)) );

  }

  // Sum up the kernel components
  for(int i=0; i<ns; i++) {
    
    Ker[i] = KT0[i] + KTA1[i] + KTB1[i] + KTA2[i] + KTB2[i] + KTA3[i] + KTB3[i];

  }

  // Switched on (To switch off, comment out lines below)
  C[0] = 1.0 / pow(*taucfl, 2.0) + 1.0 / pow(*taukls, 2.0) + 2.0 / *taucfl / *taukls + *bsi / *taukls / *tauksl;
  C[1] = -1 * *bsi / pow(*taukls, 2.0) - *bsi / *taucfl / *taukls - *bsi / *taucfs / *taukls - pow(*bsi, 2.0) / *taukls / *tauksl;
  C[2] = -1 * *bsi / pow(*tauksl, 2.0) - 1.0 / *taucfs / *tauksl - 1.0 / *taucfl / *tauksl -1.0 / *taukls / *tauksl;
  C[3] = 1.0 / pow(*taucfs, 2.0) + pow(*bsi, 2.0) / pow(*tauksl, 2.0) + 2.0 * *bsi / *taucfs / *tauksl + *bsi / *taukls / *tauksl;
  for(int i=0; i<4; i++) {
    C[i] = C[i] * (pow(double(dt), 2.0) / 12.0);
  }

  //------------------------------------------------------------------
  // Matrices of difference equation system B*T(i+1) = Q(i) + A*T(i)
  // T = (TL,TO)
  // (Equation A.27, EK05, or Equations 2.3.24 and 2.3.27, TK07)
  B[0] = 1.0 + double(dt) / (2.0 * *taucfl) + double(dt) / (2.0 * *taukls);
  B[1] = double(-dt) / (2.0 * *taukls) * *bsi;
  B[2] = double(-dt) / (2.0 * *tauksl);
  B[3] = 1.0 + double(dt) / (2.0 * *taucfs) + double(dt) / (2.0 * *tauksl) * *bsi + 2.0 * *fso * pow((double(dt) / *taudif), 0.5);

  A[0] = 1.0 - double(dt) / (2.0 * *taucfl) - double(dt) / (2.0 * *taukls);
  A[1] = double(dt) / (2.0 * *taukls) * *bsi;
  A[2] = double(dt) / (2.0 * *tauksl);
  A[3] = 1.0 - double(dt) / (2.0 * *taucfs) - double(dt) / (2.0 * *tauksl) * *bsi + Ker[ns-1] * *fso * pow((double(dt) / *taudif), 0.5);

  for (int i=0; i<4; i++) {
    B[i] = B[i] + C[i];
    A[i] = A[i] + C[i];
  }

  // Calculate the inverse of B
  invert_1d_2x2_matrix(B, IB);

  return;
}


//=========================================================================
//
// void doeclim_DICE_init(DICE *dice)
//
// Initializes the doeclim model for use in CDICE by initializing the first
// time step and running the hindcast up to 2010.
//
// Input parameters:
//   *dice = pointer to DICE object
//
//=========================================================================

void doeclim_DICE_init(DICE *dice)
{
  Clim *clim = &dice->clim;
  Config *config = &dice->config;
/*
  for(int i=0; i<clim->ns; i++) {
    clim->Ker[i] = 0.0;
    clim->temp[i] = 0.0;
    clim->temp_landair[i] = 0.0;
    clim->temp_sst[i] = 0.0;
    clim->heat_mixed[i] = 0.0;
    clim->heat_interior[i] = 0.0;
    clim->heatflux_mixed[i] = 0.0;
    clim->heatflux_interior[i] = 0.0;
    clim->forc[i] = 0.0;
  }
  for(int i=0; i<4; i++) {
		clim->IB[i] = 0.0;
    clim->A[i] = 0.0;
  }
*/
  // Initialize the doeclim model
  doeclim_init(clim->t2co, clim->kappa, config->dt, config->ns,
		 &clim->flnd, &clim->powtoheat, &clim->cal, &clim->cas,
		 &clim->taucfl, &clim->taukls, &clim->taucfs, &clim->tauksl,
		 &clim->bsi, clim->IB, clim->A, &clim->fso, &clim->taudif,
		 clim->Ker);
  
  // Run the hindcast
  //doeclim_hindcast(dice);

}


//=========================================================================
//
// void doeclim_hindcast(DICE *dice)
//
// Runs a hindcast of DOEclim using the current setup of the DICE object
//
// Input parameters:
//   *dice = pointer to DICE object
//
//=========================================================================

void doeclim_hindcast(DICE *dice)
{
  Clim *clim = &dice->clim;
  Config *config = &dice->config;
  
  for(int t=0; t<config->hind_ns; t++) {
    
    doeclim_ts(t, config->dt, config->ns, clim->forc, 
            clim->flnd, clim->powtoheat, clim->cal, clim->cas,
            clim->taucfl, clim->taukls, clim->taucfs, clim->tauksl,
            clim->bsi, clim->IB, clim->A, clim->fso,
            clim->taudif, clim->Ker, clim->temp, clim->temp_landair, clim->temp_sst, 
            clim->heat_mixed, clim->heat_interior, 
            clim->heatflux_mixed, clim->heatflux_interior);
  }
		
  // Index (config->hind_ns - 1) should now contain the initial values
  // of temperature and for running CDICE.

}


//=========================================================================
//
// void doeclim_dice_ts(DICE *dice, int ts)
//
// Runs a projection time step of the doeclim model with a DICE object.
// Uses the endogenous forcing from the DICE carbon model.
//
// Input parameters:
//   *dice = pointer to DICE object
//
//=========================================================================

void doeclim_dice_ts(DICE *dice, int ts)
{

  Clim *clim = &dice->clim;
  Carb *carb = &dice->carb;
  Config *config = &dice->config;

  // Get the range of doeclim time steps associated with this
  // DICE time step.
  int start = config->hind_ns + ((ts-1) * config->tstep);
  int finish = config->hind_ns + (ts * config->tstep);

  // Get a linear interpolation of the forcing for this time step
  double dforc = (carb->forc[ts] - carb->forc[ts-1]) / double(config->tstep);

  for(int i=start; i<finish; i++) {
    
    // Calculate the forcing needed for this doeclim time step
    clim->forc[i] = carb->forc[ts-1] + (double(i-start+1) * dforc);
    
    // Run doeclim for this doeclim time step using the interpolated forcing.
    doeclim_ts(i, config->dt, config->ns, clim->forc, 
	       clim->flnd, clim->powtoheat, clim->cal, clim->cas,
	       clim->taucfl, clim->taukls, clim->taucfs, clim->tauksl,
	       clim->bsi, clim->IB, clim->A, clim->fso,
	       clim->taudif, clim->Ker, clim->temp, clim->temp_landair, clim->temp_sst, 
	       clim->heat_mixed, clim->heat_interior, 
	       clim->heatflux_mixed, clim->heatflux_interior);
  }

}


//=========================================================================
//
// void doeclim_load_hind_forc(DICE *dice)
//
// Loads the historic forcing file into the DICE object
//
// Input parameters:
//   *dice = pointer to DICE object
//
//=========================================================================

void doeclim_load_hind_forc(DICE *dice)
{

  FILE * fp;
  fp = fopen("./giss_forc_component_1900_2015_extrap.txt", "r");

  int counter = 0;
  float temp_forc;
	float yr, ghg, o3, sh2o, refa, aie, bc, snow, stra, solar, land;
  while(!feof(fp)) {
    if (fscanf(fp, "%g %g %g %g %g %g %g %g %g %g %g", &yr, &ghg, &o3, &sh2o, &refa, &aie, &bc, &snow, &stra, &solar, &land) != 11) { break; }

    temp_forc = ghg + o3 + sh2o + stra + solar + land + (dice->clim.alpha * (refa + aie + bc + snow));

    dice->clim.forc[counter] = double(temp_forc);
    counter++;
    if (counter >= dice->config.hind_ns) {
      break;
    }
  }
  fclose(fp);
}

//=========================================================================
//
// void doeclim_load_calibration(DICE *dice, int n)
//
// Loads a file containing the calibration results. File should have 3
// columns (cs, kv, alpha).
//
// Input parameters:
//   *dice = pointer to DICE object
//   n = number of calibration points in file
//
//=========================================================================
void doeclim_load_calibration(DICE *dice, int n)
{
	// Calibration file name
	const char calfile[] = "calibration_results.txt";
	
	// Open the file and read in the calibration data
	FILE * fp;
  fp = fopen(calfile, "r");

  int counter = 0;
  float this_cs, this_kv, this_alpha;
  while(!feof(fp)) {
    if (fscanf(fp, "%g %g %g", &this_cs, &this_kv, &this_alpha) != 3) {
      break;
    }
    dice->clim.cs_cal[counter] = this_cs;
    dice->clim.kv_cal[counter] = this_kv;
    dice->clim.alpha_cal[counter] = this_alpha;
    counter++;
    if (counter >= n) {
      break;
    }
  }
  fclose(fp);
  
}


//=========================================================================
//
// void apply_calibration(DICE *dice, double this_cs)
//
// Applies the calibration for doeclim using the climate sensitivity
//
// cs, kv, and alpha are the arrays containing the calibration data. n
// is the number of elements in one of the calibration arrays (assuming
// all the calibration arrays are of the same size). this_kv and this_alpha
// are the calibrated values for the given climate sensitivity this_cs.
//
//=========================================================================
void apply_calibration(DICE *dice, double this_cs)
{
	// Loop iterator
	int i;
	
	// Setup variables to find the bounds around a particular cs
	int ihigh = 0;
	int ilow = 0;
	double diffh = 999.9;
	double diffl = -999.9;
	double this_diff = 0.0;
	double cs_range = 0.0;
	double cs_diff = 0.0;
		
	// Loop through the calibration data
	for(i=0; i<dice->config.n_calpoints; i++) {

		// Find the difference between this_cs and the calibration point
		this_diff = dice->clim.cs_cal[i] - this_cs;
		
		// If this is an exact match, set the kv and alpha values and return
		if(this_diff == 0.0) {
			ilow = i;
			ihigh = i;
			break;
		}
		
		// Otherwise, test this index as a bound
		if(this_diff >= 0.0 && this_diff < diffh) {
			ihigh = i;
			diffh = this_diff;
			}
		else if(this_diff <= 0.0 && this_diff > diffl) {
			ilow = i;
			diffl = this_diff;
		}	
	}
	
	// Now that we have the bounds, calculate the kv and alpha values
	// through a simple linear-interpolation between the bound points
	if(ilow == ihigh) {
		dice->clim.kappa = dice->clim.kv_cal[ilow];
		dice->clim.alpha = dice->clim.alpha_cal[ilow];
	}
	else {
		cs_range = dice->clim.cs_cal[ihigh] - dice->clim.cs_cal[ilow];
		cs_diff = this_cs - dice->clim.cs_cal[ilow];
		dice->clim.kappa = ((cs_diff/cs_range) * (dice->clim.kv_cal[ihigh] - dice->clim.kv_cal[ilow])) + dice->clim.kv_cal[ilow];
		dice->clim.alpha = ((cs_diff/cs_range) * (dice->clim.alpha_cal[ihigh] - dice->clim.alpha_cal[ilow])) + dice->clim.alpha_cal[ilow];
	}

	return;
}
