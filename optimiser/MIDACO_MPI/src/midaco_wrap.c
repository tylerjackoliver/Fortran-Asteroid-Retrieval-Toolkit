#include "licenseKey.h"

int firstRun = 0;

/* @brief Wraps the (C) MIDACO function in C++

	For the list of parameters, see the MIDACO documentation.
*/
int midaco_wrap(long int *p, long int *o, long int *n, long int *ni, 
	long int *m, long int *me, double *x, double *f, double *g, double *xl, double *xu, 
	long int *iflag, long int *istop, double *param, double *rw, long int *lrw, long int *iw, 
	long int *liw, double *pf, long int *lpf, long int *save2file, long int *maxeval, 
	long int *maxtime, long int *printeval){

	int returnCode = 0;

	if (firstRun == 0)
	{
		returnCode = midaco_print(1 , *printeval, *save2file, iflag, 
								istop, f, g, x, xl, 
								xu, *o, *n, *ni, *m, 
								*me, rw, pf, *maxeval, 
								*maxtime, param, *p, licenseKeyJack);
		firstRun++;
	}

	returnCode = midaco(p, o, n, ni, 
						m, me, x, f, g, 
						xl, xu, iflag, istop, 
						param, rw, lrw, iw, 
						liw, pf, lpf, licenseKeyJack);
	returnCode = midaco_print(2, *printeval, *save2file, iflag, 
							  istop, f, g, x, xl, 
							  xu, *o, *n, *ni, *m, 
							  *me, rw, pf, *maxeval, 
							  *maxtime, param, *p, licenseKeyJack);
	return returnCode;
}

/* @brief Reset whether or not we've already run MIDACO - useful for spawning successive optimisations.
*/
int reset_run()
{
	firstRun = 0;
	return 0;
}
