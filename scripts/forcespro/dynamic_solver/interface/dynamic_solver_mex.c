/*
dynamic_solver : A fast customized optimization solver.

Copyright (C) 2013-2020 EMBOTECH AG [info@embotech.com]. All rights reserved.


This software is intended for simulation and testing purposes only. 
Use of this software for any commercial purpose is prohibited.

This program is distributed in the hope that it will be useful.
EMBOTECH makes NO WARRANTIES with respect to the use of the software 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. 

EMBOTECH shall not have any liability for any damage arising from the use
of the software.

This Agreement shall exclusively be governed by and interpreted in 
accordance with the laws of Switzerland, excluding its principles
of conflict of laws. The Courts of Zurich-City shall have exclusive 
jurisdiction in case of any dispute.

*/

#include "mex.h"
#include "math.h"
#include "../include/dynamic_solver.h"
#ifndef SOLVER_STDIO_H
#define SOLVER_STDIO_H
#include <stdio.h>
#endif



/* copy functions */
void copyCArrayToM_double(double *src, double *dest, solver_int32_default dim) 
{
    solver_int32_default i;
    for( i = 0; i < dim; i++ ) 
    {
        *dest++ = (double)*src++;
    }
}
void copyMArrayToC_double(double *src, double *dest, solver_int32_default dim) 
{
    solver_int32_default i;
    for( i = 0; i < dim; i++ ) 
    {
        *dest++ = (double) (*src++) ;
    }
}
void copyMValueToC_double(double * src, double * dest)
{
	*dest = (double) *src;
}



extern void (dynamic_solver_float *x, dynamic_solver_float *y, dynamic_solver_float *l, dynamic_solver_float *p, dynamic_solver_float *f, dynamic_solver_float *nabla_f, dynamic_solver_float *c, dynamic_solver_float *nabla_c, dynamic_solver_float *h, dynamic_solver_float *nabla_h, dynamic_solver_float *hess, solver_int32_default stage, solver_int32_default iteration);
dynamic_solver_extfunc pt2function_dynamic_solver = &;


/* Some memory for mex-function */
static dynamic_solver_params params;
static dynamic_solver_output output;
static dynamic_solver_info info;

/* THE mex-function */
void mexFunction( solver_int32_default nlhs, mxArray *plhs[], solver_int32_default nrhs, const mxArray *prhs[] )  
{
	/* file pointer for printing */
	FILE *fp = NULL;

	/* define variables */	
	mxArray *par;
	mxArray *outvar;
	const mxArray *PARAMS = prhs[0];
	double *pvalue;
	solver_int32_default i;
	solver_int32_default exitflag;
	const solver_int8_default *fname;
	const solver_int8_default *outputnames[100] = {"x001","x002","x003","x004","x005","x006","x007","x008","x009","x010","x011","x012","x013","x014","x015","x016","x017","x018","x019","x020","x021","x022","x023","x024","x025","x026","x027","x028","x029","x030","x031","x032","x033","x034","x035","x036","x037","x038","x039","x040","x041","x042","x043","x044","x045","x046","x047","x048","x049","x050","x051","x052","x053","x054","x055","x056","x057","x058","x059","x060","x061","x062","x063","x064","x065","x066","x067","x068","x069","x070","x071","x072","x073","x074","x075","x076","x077","x078","x079","x080","x081","x082","x083","x084","x085","x086","x087","x088","x089","x090","x091","x092","x093","x094","x095","x096","x097","x098","x099","x100"};
	const solver_int8_default *infofields[10] = { "it", "it2opt", "res_eq", "res_ineq",  "rsnorm",  "rcompnorm",  "pobj",  "mu",  "solvetime",  "fevalstime"};
	
	/* Check for proper number of arguments */
    if (nrhs != 1) 
	{
        mexErrMsgTxt("This function requires exactly 1 input: PARAMS struct.\nType 'help dynamic_solver_mex' for details.");
    }    
	if (nlhs > 3) 
	{
        mexErrMsgTxt("This function returns at most 3 outputs.\nType 'help dynamic_solver_mex' for details.");
    }

	/* Check whether params is actually a structure */
	if( !mxIsStruct(PARAMS) ) 
	{
		mexErrMsgTxt("PARAMS must be a structure.");
	}

	/* copy parameters into the right location */
	par = mxGetField(PARAMS, 0, "x0");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	
	{
        mexErrMsgTxt("PARAMS.x0 not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.x0 must be a double.");
    }
    if( mxGetM(par) != 1200 || mxGetN(par) != 1 ) 
	{
    mexErrMsgTxt("PARAMS.x0 must be of size [1200 x 1]");
    }
#endif	 
	if ( (mxGetN(par) != 0) && (mxGetM(par) != 0) )
	{
		copyMArrayToC_double(mxGetPr(par), params.x0,1200);

	}
	par = mxGetField(PARAMS, 0, "xinit");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	
	{
        mexErrMsgTxt("PARAMS.xinit not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.xinit must be a double.");
    }
    if( mxGetM(par) != 9 || mxGetN(par) != 1 ) 
	{
    mexErrMsgTxt("PARAMS.xinit must be of size [9 x 1]");
    }
#endif	 
	if ( (mxGetN(par) != 0) && (mxGetM(par) != 0) )
	{
		copyMArrayToC_double(mxGetPr(par), params.xinit,9);

	}
	par = mxGetField(PARAMS, 0, "all_parameters");
#ifdef MEXARGMUENTCHECKS
    if( par == NULL )	
	{
        mexErrMsgTxt("PARAMS.all_parameters not found");
    }
    if( !mxIsDouble(par) )
    {
    mexErrMsgTxt("PARAMS.all_parameters must be a double.");
    }
    if( mxGetM(par) != 1200 || mxGetN(par) != 1 ) 
	{
    mexErrMsgTxt("PARAMS.all_parameters must be of size [1200 x 1]");
    }
#endif	 
	if ( (mxGetN(par) != 0) && (mxGetM(par) != 0) )
	{
		copyMArrayToC_double(mxGetPr(par), params.all_parameters,1200);

	}


	#if SET_PRINTLEVEL_dynamic_solver > 0
		/* Prepare file for printfs */
        fp = fopen("stdout_temp","w+");
		if( fp == NULL ) 
		{
			mexErrMsgTxt("freopen of stdout did not work.");
		}
		rewind(fp);
	#endif

	/* call solver */
	exitflag = dynamic_solver_solve(&params, &output, &info, fp, pt2function_dynamic_solver);
	
	#if SET_PRINTLEVEL_dynamic_solver > 0
		/* Read contents of printfs printed to file */
		rewind(fp);
		while( (i = fgetc(fp)) != EOF ) 
		{
			mexPrintf("%c",i);
		}
		fclose(fp);
	#endif

	/* copy output to matlab arrays */
	plhs[0] = mxCreateStructMatrix(1, 1, 100, outputnames);
		outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x001, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x001", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x002, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x002", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x003, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x003", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x004, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x004", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x005, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x005", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x006, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x006", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x007, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x007", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x008, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x008", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x009, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x009", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x010, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x010", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x011, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x011", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x012, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x012", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x013, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x013", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x014, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x014", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x015, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x015", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x016, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x016", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x017, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x017", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x018, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x018", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x019, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x019", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x020, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x020", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x021, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x021", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x022, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x022", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x023, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x023", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x024, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x024", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x025, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x025", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x026, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x026", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x027, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x027", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x028, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x028", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x029, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x029", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x030, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x030", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x031, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x031", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x032, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x032", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x033, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x033", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x034, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x034", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x035, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x035", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x036, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x036", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x037, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x037", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x038, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x038", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x039, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x039", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x040, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x040", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x041, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x041", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x042, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x042", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x043, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x043", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x044, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x044", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x045, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x045", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x046, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x046", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x047, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x047", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x048, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x048", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x049, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x049", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x050, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x050", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x051, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x051", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x052, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x052", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x053, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x053", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x054, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x054", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x055, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x055", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x056, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x056", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x057, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x057", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x058, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x058", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x059, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x059", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x060, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x060", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x061, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x061", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x062, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x062", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x063, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x063", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x064, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x064", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x065, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x065", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x066, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x066", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x067, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x067", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x068, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x068", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x069, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x069", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x070, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x070", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x071, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x071", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x072, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x072", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x073, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x073", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x074, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x074", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x075, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x075", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x076, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x076", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x077, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x077", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x078, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x078", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x079, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x079", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x080, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x080", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x081, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x081", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x082, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x082", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x083, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x083", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x084, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x084", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x085, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x085", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x086, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x086", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x087, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x087", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x088, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x088", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x089, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x089", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x090, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x090", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x091, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x091", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x092, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x092", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x093, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x093", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x094, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x094", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x095, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x095", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x096, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x096", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x097, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x097", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x098, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x098", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x099, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x099", outvar);

	outvar = mxCreateDoubleMatrix(12, 1, mxREAL);
	copyCArrayToM_double( output.x100, mxGetPr(outvar), 12);
	mxSetField(plhs[0], 0, "x100", outvar);

	/* copy exitflag */
	if( nlhs > 1 )
	{
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(plhs[1]) = (double)exitflag;
	}

	/* copy info struct */
	if( nlhs > 2 )
	{
	        plhs[2] = mxCreateStructMatrix(1, 1, 10, infofields);
         
		
		/* iterations */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)info.it;
		mxSetField(plhs[2], 0, "it", outvar);

		/* iterations to optimality (branch and bound) */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = (double)info.it2opt;
		mxSetField(plhs[2], 0, "it2opt", outvar);
		
		/* res_eq */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.res_eq;
		mxSetField(plhs[2], 0, "res_eq", outvar);

		/* res_ineq */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.res_ineq;
		mxSetField(plhs[2], 0, "res_ineq", outvar);

		/* rsnorm */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.rsnorm;
		mxSetField(plhs[2], 0, "rsnorm", outvar);

		/* rcompnorm */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.rcompnorm;
		mxSetField(plhs[2], 0, "rcompnorm", outvar);
		
		/* pobj */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.pobj;
		mxSetField(plhs[2], 0, "pobj", outvar);

		/* mu */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.mu;
		mxSetField(plhs[2], 0, "mu", outvar);

		/* solver time */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.solvetime;
		mxSetField(plhs[2], 0, "solvetime", outvar);

		/* solver time */
		outvar = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(outvar) = info.fevalstime;
		mxSetField(plhs[2], 0, "fevalstime", outvar);
	}
}