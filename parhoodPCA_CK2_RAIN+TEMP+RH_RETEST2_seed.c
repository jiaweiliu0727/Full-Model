#include "headW2.h"
gsl_rng *r;
//////////////////Begin Dot Product////////////////////////////////////////

float DotProduct (int Length, double *Holder, double *PCA)
{
//	printf("here\n");
//	printf("%f %f\n", Coefficients[0], PC[0]);
  double answer = 0;

  int i;
  for(i=0;i<Length;i++)
    answer += PCA[i]*Holder[i];

  return(answer);


}

//////////////////End dot product/////////////////////////////////////////



int main(void)
{
//int test = atoi(argv[1]);                  // test is 99 for checking, 0 otherwise

int test = 0; //CK// Second test mode.  Like the full program but less runs and less MISER calls

//int test = 99;
//int test = 66; //CK// Second test mode.  Like the full program but less runs and less MISER calls

////////////
///NOTE:  Going to try to change nuF to a site-specific parameter and make R(0) a general parameter
////////// Just going to try making a new parameter for average R(0).  Make it number 23.
////////// Comment out param #3 (general nuF) but leave structure in place if we go back to it
////////// Pipe param #23 into R for all populations
////////// Make it fit nuF, rather than R(0), for each population.  Sounds simple enough....

STRUCTURE Params;

int pro = 1;//atoi(argv[1]);						// pro and argv[1] are the inputs (argv[i] is the i^th input)
printf("Profile Parameter is %d\n",pro);	fflush(stdout);
// ------------------------------------- Adustable accuracy vs. speed ------------------------------------------------ //
int num_runs	 = 20;
double parm_inc, host_inc, initR_inc;	//int inc_gamma_box= 1;

//if (pro==1)	{	parm_inc=200.0;		host_inc=100.0;	initR_inc=100.0;	}
if (pro==1)	{	parm_inc=32.0;		host_inc=10.0;	initR_inc=18.0;	}
//if (pro==1)	{	parm_inc=34.0;		host_inc=18.0;	initR_inc=20.0;	}
else		{	parm_inc=15.0;		host_inc=10.0;	initR_inc=10.0;	}

//if (test==66)	{	printf("for checking CK MODE2!!!\n");        num_runs=5;} //CK// New test mode!

if (test==66)	{	printf("for checking!!!\n");        num_runs=5;	parm_inc=16.0;	host_inc=8.0;	initR_inc=8.0;}

//if (test==66)	{	printf("for checking!!!\n");        num_runs=5;	parm_inc=6.0;	host_inc=4.0;	initR_inc=4.0;	}
printf("runs=%d\t incs: parm=%2.0f\t S_0=%2.0f\t R_0=%2.0f\n",num_runs,parm_inc,host_inc,initR_inc);

//---------------//CK// Best fit params and initial conditions from previous run //CK//-------------------------------------------------------//

//double bestparams[19]={1.000000e+02,	1.000000e+00,	8.254042e-01,	2.154435e-05,	1.000000e-02,	3.981072e-01,	5.623413e-01,	0.000000e+00,	5.000000e+00,	5.000000e+00,	1.000000e-01,	3.000000e-01,	3.000000e-01,	1.000000e-02,	2.618589e+00,	1.191667e+01,	2.040000e+00,	2.083333e-02,	1.000000e-01};
//double initialS[6]={1.2, 33.2, 4.8, 61.4, 45.4, 154.6};	 //CK// S(0)
//double initialR[6]={1.500000e-01,	1.000001e-02,	3.750000e-01,	6.250009e-03,	7.312500e-02,	6.250009e-03};	//CK// R(0)

double initialR[8+1];

//---------------//CK// Initial conditions for S from field observations //CK//-------------------------------------------------------//

//double initialS[7]={1.92, 53.1, 7.68, 98.2, 72.6, 6.72, 5.76};	 //CK// my observed average egg density per circle at each site S(0)
//double initialS[8]={1.92*10, 98.2, 7.68*5, 72.6, 53.1, 6.72, 5.76, 247.36};	 //CK// my observed average egg density per circle at each site S(0)
//double initialS[8]={14.22, 242.43125, 18.96, 179.23125, 131.090625, 5.53, 4.74, 610.67};	 //CK// my observed average egg density per circle at each
double initialS[8]={14.22*2, 242.43125, 18.96*3, 179.23125, 131.090625, 5.53, 4.74, 610.67};	 //CK// my observed average egg density per circle at each

// ------------------------------------------------------------------------------------------------------------------ //
int i=0; int j;int ii; int jj; int k;
int run;	            int changer;	    double index, tot_index;

int num_adj_pars=29;			// number of adjustable parameters

double inner_parm;				//double outer_parm;
double inner_parm2;				//double outer_parm;
int pop;
double log_pop;
//Params.th_id=0;
// -------------------------------------------- MISTER STUFF --------------------------------------------------------- //
inputdata(&Params);				// gets Params.DATA[j][i][0-2] and Params.MAXT[i] from inputdata.h
//int calls=2;					//CK// turned down to run without stochasticity
int calls=400;					// number of stochastic simulations for each parameter and IC set
if (test==99)	calls=10;
if (test==66)	calls=2; //CK// second test mode
//if (test==66)	calls=5; //CK// second test mode
size_t dim;
size_t dim2;
// --------------------------------------- Name for Output Files ----------------------------------------------------- //
char strFileName[99];					// from filenames.h
GetString(pro,0,strFileName,98);	printf("file names:\t   %s\n",strFileName);		fflush(stdout);		//getc(stdin);
FILE *fp_results;
// ---------------------------------------- Random Number Stuff ------------------------------------------------------ //
gsl_rng *r_seed;
r_seed=random_setup();
//printf("Random Seed: %f\n", r_seed); getc(stdin);
// -------------------- parameter high/low values and increments and fixed parameter values ------------------------- //
global_fixed_parms(&Params);  // gets Params.PARS[i] for fixed parameters from bounds.h
parm_range_inc(&Params,parm_inc,host_inc,initR_inc,num_adj_pars); // gets Params.parm_set,low,high,R_END from bounds.h
// ------------------------------------ Declare Likelihood Quanitites ----------------------------------------------- //
double pop_lhood, pop_lhood2, pop_err,post_hood;	// population lhood (and posterior lhood) calculated for each initS and initR
double pop_best_lhood;					// likelihood and error for best initS and initR
double total_lhood;						// sum of pop_best_lhood over all patches
double best_post_hood;	double best_lhood=0;		// best post_hood and lhood
double prior[num_adj_pars];


/////////////////////////////Begin reading in principal component analysis results//////////////////////

	int NumberOfParams=28;			// 32 parameters

	int Realizations=25000;

	double RandomNumber;

	double LogJumpToNew;
	double LogJumpToOld;
	double ProbOfAcceptance;

	double LogOldPosterior;
	double LogNewPosterior;
	double LogNewPrior;

	const gsl_rng_type * T;

	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r, time(NULL)*(int)getpid());
	srand(time(NULL));

	double delta=0;

	double SDpca[NumberOfParams];
	double Coefficients[NumberOfParams*NumberOfParams];
	double Center[NumberOfParams];
	double Scale[NumberOfParams];
	double PC[NumberOfParams];
	double Old_PC[NumberOfParams];
	double Holder[NumberOfParams];
	double PCAparams[NumberOfParams];
	double Old_Params[NumberOfParams];

	//double AcceptedVect[NumberOfParams]={0.0};
	//double LoopVect[NumberOfParams]={0.0};

	int a;
	int b;
	int ticker;
	int ticker2;

	int Case;

	int Accepted=0;
	signed int LoopNumber=-1;

	int ParCnt2 = NumberOfParams;

	double SigmaInflation=1.20;  //CK// Dave had 4 as his sigma inflation factor - FORMERLY WAS 0.5, WHICH IS A BAD THING

	double sigma[NumberOfParams];

	run=1;	changer=1;	best_post_hood=-10000000000;

/////////////////////////////Begin reading in principal component analysis results//////////////////////

	FILE *file;

	file=fopen("PCAsdCK2.txt", "r");
	for (a=0; a<(NumberOfParams); a++)
	{
		fscanf(file, "%lf\n", &SDpca[a]);
		//printf("%lf\n", SDpca[a]);
	}
	fclose(file);


	file=fopen("PCArotationsCK2.txt", "r");
	for (a=0; a<(NumberOfParams*NumberOfParams); a++)
	{
		fscanf(file, "%lf\n", &Coefficients[a]);
		//printf("%lf\n", Coefficients[a]);
	}
	fclose(file);

	file=fopen("PCAscaleCK2.txt", "r");
	for (a=0; a<NumberOfParams; a++)
	{
		fscanf(file, "%lf\n", &Scale[a]);
		//printf("%lf\n", Scale[a]);
	}
	fclose(file);

	file=fopen("PCAcenterCK2.txt", "r");

	for (a=0; a<NumberOfParams; a++)
	{
		fscanf(file, "%lf\n", &Center[a]);
		//printf("%lf\n", Center[a]);
	}
	fclose(file);

/////////////////////////////DONE reading in principal component analysis results//////////////////////


//////////////////////////////Start MCMC/////////////////////////////////////////////////////////////////////////////

///Generate first set of PC values////

	for(a=0; a<NumberOfParams; a++){
		sigma[a]=SigmaInflation*SDpca[a];
		//printf("%lf\n", sigma[a]);
	}


	for (a=0; a<NumberOfParams; a++){
		PC[a]=gsl_ran_gaussian (r, sigma[a]);    //CK//  PC contains current set of PC values
		//printf("%lf\n", PC[a]);
	}
	//getc(stdin);

///Begin Loop that goes through and tweaks PCs one at a time////

//while (1==1) {     //CK// INFINITE LOOP START!!!!
while (LoopNumber<=Realizations) {     //CK// BOUND LOOP START!!!!

		LoopNumber=LoopNumber+1;
//if (LoopNumber % NumberOfParams == 1)
//		{
		Case=LoopNumber%NumberOfParams;					//Determines which PC to change

		for (a=0; a<NumberOfParams; a++)
		{
			Old_PC[a]=PC[a];							//Store old PC values:
		}

		PC[Case]=gsl_ran_gaussian (r, sigma[Case]);		//Draw 1 new PC

		for (a=0;a<NumberOfParams; a++)								//Back transform PC's into model parameters
		{
			for (b=0; b<NumberOfParams; b++){
				Holder[b]=Coefficients[a*NumberOfParams+b];
			}


			PCAparams[a]=exp(DotProduct(ParCnt2, Holder, PC)*Scale[a]+Center[a]);
			Old_Params[a]=exp(DotProduct(ParCnt2, Holder, Old_PC)*Scale[a]+Center[a]);

		    //printf("%f\n", PC[a]);
		    //printf("%f\n", exp(PCAparams[a]));
		    //printf("PCAparams: %e\n", PCAparams[a]);  //getc(stdin);

		}

		LogJumpToNew= -log(gsl_ran_gaussian_pdf(PC[Case], sigma[Case]));		//Metropolis sampling step: LATER YOU WILL USE THIS TO CORRECT FOR PROPOSAL
		LogJumpToOld= -log(gsl_ran_gaussian_pdf(Old_PC[Case], sigma[Case]));

		//printf("New: %f Old: %f\n", LogJumpToNew, LogJumpToOld);

		//---------------------LINE UP THE NEW PCAparams WITH THE CURRENT ORDER OF PARAMS --------------------------//

		ticker=0;

		for (k=0;k<=num_adj_pars;k++)	{

			if (k==2||k==4||k==5||k==8||k==10||k==15)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 0.0;
			}
			//else if (k==6)	{	//CK// set conidia parameters to 0. NO CONIDIA
			//	Params.PARS[k] = 0.0;
			//}
			else if (k==0||k==1)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 1.0;
			}

			else if (k==7)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 10.0;
			}
			else if (k==9)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 50.0;
			}

			else{Params.PARS[k] = PCAparams[ticker];
				ticker++;}

			//printf("NEW Params.PARS: %e\n",Params.PARS[k]);
		}
//getc(stdin);

		//if(Params.PARS[11]>1.0){Params.PARS[11] = 1.0;}
		//if(Params.PARS[12]>1.0){Params.PARS[12] = 1.0;}
		//if(Params.PARS[13]>10.0){Params.PARS[13] = 10.0;}
		//if(Params.PARS[14]>10.0){Params.PARS[14] = 10.0;}
		//if(Params.PARS[21]>0.10){Params.PARS[21] = 0.10;}
		//if(Params.PARS[22]>0.251){Params.PARS[22] = 0.251;}
		//if(Params.PARS[23]>0.251){Params.PARS[23] = 0.251;}
		//if(Params.PARS[26]>0.251){Params.PARS[26] = 0.251;}


//		exit(0);

		for (j=1;j<=DATA_SETS;j++)	{

			ticker2= NumberOfParams - DATA_SETS + j - 1;

			initialR[j] = PCAparams[ticker2];
			//printf("ssNuF: %e\n",initialR[j]);
		}

		//printf("New: %f Old: %f\n", LogJumpToNew, LogJumpToOld);

/*			// ------------------ Loop over gamma for each global param set ------------------------------------------- //
				for (inner_parm2=Params.parm_low[14];inner_parm2<=Params.parm_high[14];inner_parm2+=Params.parm_step[14])	{
					Params.PARS[14] = pow(10,inner_parm2);
					//printf("LINE SEARCH FOR PARAMETER:%d\t value=%e\n",14,Params.PARS[14]); //getc(stdin);
*/
				total_lhood=0;

				// ----------------------- loop over patch numbers -------------------------------------------- //
				for (Params.pop=1;Params.pop<=DATA_SETS;Params.pop++)	{
					pop=Params.pop;
					pop_best_lhood = -1e9;

					int MAXT3=(Params.EXPDATA[pop][Params.MAXT2[pop]][2]+1)*7;
					if(pop==6){MAXT3=77;}  ///can't figure out why plot 6 (UMBS 2012) is fucked up.  says it only has 2 weeks of data, but that's not true at all...

					dim = 2*MAXT3;     //CK//  Changed dim to accomodate the longer EXP data sets... hopefully that works

					//dim = 2*Params.MAXT[pop];		// need to let this vary by patch (double because of nuV and nuF)
					gsl_monte_function G = { &Hood_Pops, dim, &Params };	// declares function calling Hood_Pops.h
					double xl[dim];	double xu[dim];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim;jj++)	{
						xl[jj]=0;
						xu[jj]=1;
					}

					Params.PARS[30+pop]=initialS[pop-1];  //CK//  I think these should just be initialS[pop], not initialS[pop-1].  Keep everything where pop starts at 1
					Params.PARS[50+pop]=initialR[pop];  //CK// Needs to be [pop-1].  initialS is where the conditions are read in, so it starts at 0 and needs to be adjusted for
					//Params.PARS[50+pop]=initialR[pop-1];  //CK// Needs to be [pop-1].  initialS is where the conditions are read in, so it starts at 0 and needs to be adjusted for

					//printf("S(0): %e\n",Params.PARS[30+pop]);
					//printf("ssNuF: %e\n",Params.PARS[50+pop]);

						// ----------- Use MISER to call function pop_lhood --------------------------------- //
						gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim);
						gsl_monte_miser_integrate (&G,xl,xu,dim,calls,r_seed,s,&pop_lhood,&pop_err);
						gsl_monte_miser_free(s);
							//printf("pop:%d lhood=%f log_error=%f error=%f best_post_hood=%f\n",
								//pop,pop_lhood,log(pop_err),pop_err,pop_best_lhood);
							// ------------------- check to see if these ICs are better  ------------------------ //
						pop_lhood2=log(pop_lhood)-700.0;  //CK// Converting back to log likelihoods for MCMC
							//printf("pop_lhood:%e\t pop_lhood2:%f\n",pop_lhood,pop_lhood2); //getc(stdin);

							//printf("pop_lhood:%e\t pop_lhood2:%f\t pop_best_lhood=%f\n",pop_lhood,pop_lhood2,pop_best_lhood); getc(stdin);
						total_lhood += pop_lhood2;
					//printf("total_lhood=%f\n",total_lhood);	getc(stdin);
				} //CK// end of going through patches

				//prior[5]=log(prior_dist(5,log10(Params.PARS[5])));	//printf("prior[5]=%f\n",prior[5]);
				//prior[6]=log(prior_dist(6,log10(Params.PARS[6])));	//printf("prior[6]=%f\n",prior[6]);
				//post_hood=total_lhood+prior[5]+prior[6];
				LogNewPosterior = 0.0;
				LogNewPosterior=total_lhood;
				//printf("lhood=%f\t param(5)=%f\t prior(5)=%e\t params(6)=%f\t prior(6)=%e\n",LogNewPosterior,log10(Params.PARS[5]),prior[5],log10(Params.PARS[6]),prior[6]);
				//printf("\t parm:%d posthood=%f\t best_post_hood=%f\t prior5=%f\t prior6=%f\n",i,LogNewPosterior,best_post_hood,prior[5],prior[6]);	//getc(stdin);

		//---------------------Done calculating likelihood of NEW PARAM set --------------------------//

		//---------------------LINE UP THE OLD Params WITH THE CURRENT ORDER OF PARAMS --------------------------//

		ticker=0;

		for (k=0;k<=num_adj_pars;k++)	{

			if (k==2||k==4||k==5||k==8||k==10||k==15)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 0.0;
			}
			//else if (k==6)	{	//CK// set conidia parameters to 0. NO VIRUS
			//	Params.PARS[k] = 0.0;
			//}
			else if (k==0||k==1)	{	//CK// set virus parameters to 0. NO CONIDIA
				Params.PARS[k] = 1.0;
			}

			else if (k==7)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 10.0;
			}
			else if (k==9)	{	//CK// set virus parameters to 0. NO VIRUS
				Params.PARS[k] = 50.0;
			}
			else{Params.PARS[k] = Old_Params[ticker];                //CK// PUT IN OLD_PARAMS HERE!
				ticker++;}

			//printf("Params.PARS: %e\n",Params.PARS[k]);
		}

		//if(Params.PARS[11]>1.0){Params.PARS[11] = 1.0;}
		//if(Params.PARS[12]>1.0){Params.PARS[12] = 1.0;}
		//if(Params.PARS[13]>10.0){Params.PARS[13] = 10.0;}
		//if(Params.PARS[14]>10.0){Params.PARS[14] = 10.0;}
		//if(Params.PARS[21]>0.10){Params.PARS[21] = 0.10;}
		//if(Params.PARS[22]>0.251){Params.PARS[22] = 0.251;}
		//if(Params.PARS[23]>0.251){Params.PARS[23] = 0.251;}

		for (j=1;j<=DATA_SETS;j++)	{

			ticker2= NumberOfParams - DATA_SETS + j - 1;

			initialR[j] = Old_Params[ticker2];       //CK// PUT IN OLD_PARAMS HERE!
			//printf("ssNuF: %e\n",initialR[j]);
		}

		//printf("New: %f Old: %f\n", LogJumpToNew, LogJumpToOld);

/*			// ------------------ Loop over gamma for each global param set ------------------------------------------- //
				for (inner_parm2=Params.parm_low[14];inner_parm2<=Params.parm_high[14];inner_parm2+=Params.parm_step[14])	{
					Params.PARS[14] = pow(10,inner_parm2);
					//printf("LINE SEARCH FOR PARAMETER:%d\t value=%e\n",14,Params.PARS[14]); //getc(stdin);
*/
				total_lhood=0;

				// ----------------------- loop over patch numbers -------------------------------------------- //
				for (Params.pop=1;Params.pop<=DATA_SETS;Params.pop++)	{
					pop=Params.pop;
					pop_best_lhood = -1e9;

					int MAXT3=(Params.EXPDATA[pop][Params.MAXT2[pop]][2]+1)*7;
					if(pop==6){MAXT3=77;}  ///can't figure out why plot 6 (UMBS 2012) is fucked up.  says it only has 2 weeks of data, but that's not true at all...

					dim = 2*MAXT3;     //CK//  Changed dim to accomodate the longer EXP data sets... hopefully that works

					//dim = 2*Params.MAXT[pop];		// need to let this vary by patch (double because of nuV and nuF)
					gsl_monte_function G = { &Hood_Pops, dim, &Params };	// declares function calling Hood_Pops.h
					double xl[dim];	double xu[dim];	// need to redeclare xl and xu since the size changes
					for (jj=0;jj<=dim;jj++)	{
						xl[jj]=0;
						xu[jj]=1;
					}

					Params.PARS[30+pop]=initialS[pop-1];  //CK//  I think these should just be initialS[pop], not initialS[pop-1].  Keep everything where pop starts at 1
					Params.PARS[50+pop]=initialR[pop];  //CK// Needs to be [pop-1].  initialS is where the conditions are read in, so it starts at 0 and needs to be adjusted for

					//printf("S(0): %e\n",Params.PARS[30+pop]);
					//printf("ssNuF: %e\n",Params.PARS[50+pop]);

						// ----------- Use MISER to call function pop_lhood --------------------------------- //
						gsl_monte_miser_state *s = gsl_monte_miser_alloc(dim);
						gsl_monte_miser_integrate (&G,xl,xu,dim,calls,r_seed,s,&pop_lhood,&pop_err);
						gsl_monte_miser_free(s);
							//printf("pop:%d lhood=%f log_error=%f error=%f best_post_hood=%f\n",
								//pop,pop_lhood,log(pop_err),pop_err,pop_best_lhood);
							// ------------------- check to see if these ICs are better  ------------------------ //
						pop_lhood2=log(pop_lhood)-700.0;  //CK// Converting back to log likelihoods for MCMC
							//printf("pop_lhood:%e\t pop_lhood2:%f\n",pop_lhood,pop_lhood2); //getc(stdin);

							//printf("pop_lhood:%e\t pop_lhood2:%f\t pop_best_lhood=%f\n",pop_lhood,pop_lhood2,pop_best_lhood); getc(stdin);
						total_lhood += pop_lhood2;
					//printf("total_lhood=%f\n",total_lhood);	getc(stdin);
				} //CK// end of going through patches

				//prior[5]=log(prior_dist(5,log10(Params.PARS[5])));	//printf("prior[5]=%f\n",prior[5]);
				prior[6]=log(prior_dist(6,log10(Params.PARS[6])));	//printf("prior[6]=%f\n",prior[6]);
				//post_hood=total_lhood+prior[5]+prior[6];
				LogOldPosterior=0.0;
				LogOldPosterior=total_lhood;
				//printf("lhood=%f\t param(5)=%f\t prior(5)=%e\t params(6)=%f\t prior(6)=%e\n",total_lhood,log10(Params.PARS[5]),prior[5],log10(Params.PARS[6]),prior[6]);
				//printf("\t parm:%d posthood=%f\t best_post_hood=%f\t prior5=%f\t prior6=%f\n",i,post_hood,best_post_hood,prior[5],prior[6]);	//getc(stdin);

		//---------------------Done calculating likelihood of OLD PARAM set --------------------------//

		// ----------------------Compare New Param set with Old Param Set ----------------------------------------- //

		//ProbOfAcceptance=exp(-LogNewPosterior-LogJumpToOld + LogOldPosterior+LogJumpToNew);    //Probability of accepting the new PC
		ProbOfAcceptance=exp(LogNewPosterior+LogJumpToOld - LogOldPosterior-LogJumpToNew);    //Probability of accepting the new PC
//		printf("Reasons to stay=%.0f\n", LogNewPosterior+LogJumpToOld);
//		printf("Reasons to leave=%.0f\n",LogOldPosterior+LogJumpToNew);
//		printf("ProbOfAcceptance=%f\n", ProbOfAcceptance);
//The larger the value, the more likely to accept:  -LogNewPosterior	-LogJumpToOld 	+ LogOldPosterior	+LogJumpToNew


		//printf("LogOldPosterior: %f\t LogNewPosterior: %f\t ProbOfAcceptance: %f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance); getc(stdin);

		Params.LoopVect[Case] = Params.LoopVect[Case] + 1;

		if (ProbOfAcceptance>1 || gsl_rng_uniform_pos (r) < ProbOfAcceptance)   //MH-MCMC algorithm
		{
			LogOldPosterior=LogNewPosterior;
			//Accepted=Accepted+1;
			//printf("Accepted\n");
			Params.AcceptedVect[Case] = Params.AcceptedVect[Case] + 1;
		}

		else
		{
			for (a=0; a<NumberOfParams; a++)
			{
				PC[a]=Old_PC[a];
			}

			//printf("Rejected\n");

		}

		//getc(stdin);

		//exit(0);

		//printf("%f\t%f\t%f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance); getc(stdin);

/*		if (LoopNumber % 100 ==0)   	//thin interval, set to 1 for unthinned output
		{
			//printf("%f\t%f\t%f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance);
			//printf("%f %f %f %f %f %f %f\n", PC[0], PC[1], PC[2], PC[3], PC[4], (double)Accepted/((double)LoopNumber+1.0), LogOldPosterior-LogNewPrior);  	//print to screen
			fflush(stdout);
		}
*/
		if (LoopNumber % 10 == 0)   	//CK// output results every 5000 loops (probably not best plan but we'll see)
		{
		// ------------------------------------------ output results to file  --------------------------------------- //
			//Params.PARS[0] = Accepted;
			//Params.PARS[1] = LoopNumber;

			//printf("%f\t%f\t%f\n", LogOldPosterior, LogNewPosterior, ProbOfAcceptance);

			fp_results = fopen(strFileName,"a");
			output_file(&Params,fp_results,LogOldPosterior,num_adj_pars,pro); // prints to output file (filenames.h)
			fclose(fp_results);

			fflush(stdout);
		// ------------------------------------------------------------------------------------------------------- //
		}


}   //end of the infinite while loop

/*
				if (post_hood > best_post_hood)				{
					best_post_hood = post_hood;
					best_lhood= total_lhood;				// max lhood associated with this posterior
					Params.MLE[i] = Params.PARS[i];
					Params.MLE[14] = Params.PARS[14];

					//printf("\t switch to new MLEs!\t parm:%d has new MLE value=%e with best post_hood=%f\n",i,Params.MLE[i],best_post_hood);
					//printf("\t switch to new MLEs!\t parm:%d has new MLE value=%e with best post_hood=%f\n",14,Params.MLE[14],best_post_hood);

					//if (test==66){printf("\t switch to new MLEs!\t parm:%d has new MLE value=%e with best post_hood=%f\n",i,Params.MLE[i],best_post_hood);}
					// ------------- update MLE values for initS and initR ------------------------------ //
					for (pop=1;pop<=DATA_SETS;pop++)	{
						Params.MLE_host[pop] = Params.best_initS[pop];
						Params.MLE_initR[pop]= Params.best_initR[pop];
						//printf("pop:%d\t S=%e\t R=%e\n",pop,Params.MLE_host[pop],Params.MLE_initR[pop]);
					}//getc(stdin);
					changer++;
				}
				// ---------------------------------------------------------------------------------------- //
*/
//	}	//CK// END OF LINE SEARCH FOR GAMMA



//}		//CK// END OF LINE SEARCH FOR SELECTED PARAMETER
			//printf("parameter: %d\t has MLE value=%f\n",i,Params.MLE[i]);
	//		Params.PARS[i]=Params.MLE[i];     //////CK///////!!!!!!!!  Found the spot where it puts the best value back into PARS for continued checking
	//		Params.PARS[14]=Params.MLE[14];     //////CK///////!!!!!!!!  Found the spot where it puts the best value back into PARS for continued checking
	//	}		//CK/  END OF IF STATEMENT FOR SPECIFIC PARAM NUMBERS
	//}		//CK// END OF LOOP FOR TICKING THROUGH PARAM NUMBERS FOR LINE SEARCHES
	//printf("profile=%d: value=%f\t lhood=%f\t post_hood=%f\t changer=%d\n",pro,Params.PARS[pro],best_lhood,best_post_hood,changer);
//}		//CK// END OF INNER LOOP!!

index=(Params.PARS[pro]-Params.parm_low[pro])/Params.parm_step[pro];
tot_index=(Params.parm_high[pro]-Params.parm_low[pro])/Params.parm_step[pro];
printf("INDEX:%3.0f of %3.0f\t VALUE:%f\t LHOOD:%f\n",index+1,tot_index,Params.PARS[pro],LogOldPosterior);

// ------------------------------------------ output results to file  --------------------------------------- //
//fp_results = fopen(strFileName,"a");
//output_file(&Params,fp_results,LogOldPosterior,num_adj_pars,pro); // prints to output file (filenames.h)
//fclose(fp_results);
// ------------------------------------------------------------------------------------------------------- //
//CK//}   //end of the infinite while loop
free_i3tensor(Params.DATA,0,DATA_SETS,0,MAX_WEEKS,0,3);
free_i3tensor(Params.EXPDATA,0,DATA_SETS,0,MAX_WEEKS,0,3);
//free_i3tensor(Params.WDATA,0,DATA_SETS,0,MAX_WEEKS2,0,3);
printf("DONE!!!\n");

return 0;
}
