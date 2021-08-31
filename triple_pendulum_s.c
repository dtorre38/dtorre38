/*
Generated 09-Aug-2008 11:36:26 by SD/FAST, Kane's formulation
(sdfast B.2.8 #30123) on machine ID unknown
Copyright (c) 1990-1997 Symbolic Dynamics, Inc.
Copyright (c) 1990-1997 Parametric Technology Corp.
RESTRICTED RIGHTS LEGEND: Use, duplication, or disclosure by the U.S.
Government is subject to restrictions as set forth in subparagraph
(c)(1)(ii) of the Rights in Technical Data and Computer Software
clause at DFARS 52.227-7013 and similar clauses in the FAR and NASA
FAR Supplement.  Symbolic Dynamics, Inc., Mountain View, CA 94041
*/
#include <math.h>

/* These routines are passed to sdroot. */

void sdposfunc(double vars[3],
    double param[1],
    double resid[1])
{
    int i;
    double pos[3],vel[3];

    for (i = 0; i < 3; i++) {
        vel[i] = 0.;
    }
    sdang2st(vars,pos);
    sdstate(param[0],pos,vel);
    sdperr(resid);
}

void sdvelfunc(double vars[3],
    double param[4],
    double resid[1])
{

    sdstate(param[3],param,vars);
    sdverr(resid);
}

void sdstatfunc(double vars[3],
    double param[4],
    double resid[3])
{
    double pos[3],qdotdum[3];

    sdang2st(vars,pos);
    sdstate(param[3],pos,param);
    sduforce(param[3],pos,param);
    sdperr(resid);
    sdderiv(qdotdum,&resid[0]);
}

void sdstdyfunc(double vars[6],
    double param[1],
    double resid[3])
{
    double pos[3],qdotdum[3];

    sdang2st(vars,pos);
    sdstate(param[0],pos,&vars[3]);
    sduforce(param[0],pos,&vars[3]);
    sdperr(resid);
    sdverr(&resid[0]);
    sdderiv(qdotdum,&resid[0]);
}

/* This routine is passed to the integrator. */

void sdmotfunc(double time,
    double state[6],
    double dstate[6],
    double param[1],
    int *status)
{

    sdstate(time,state,&state[3]);
    sduforce(time,state,&state[3]);
    sdderiv(dstate,&dstate[3]);
    *status = 0;
}

/* This routine performs assembly analysis. */

void sdassemble(double time,
    double state[6],
    int lock[3],
    double tol,
    int maxevals,
    int *fcnt,
    int *err)
{
    double perrs[1],param[1];
    int i;

    sdgentime(&i);
    if (i != 113626) {
        sdseterr(50,42);
    }
    param[0] = time;
    *err = 0;
    *fcnt = 0;
    sdposfunc(state,param,perrs);
    *fcnt = *fcnt+1;
}

/* This routine performs initial velocity analysis. */

void sdinitvel(double time,
    double state[6],
    int lock[3],
    double tol,
    int maxevals,
    int *fcnt,
    int *err)
{
    double verrs[1],param[4];
    int i;

    sdgentime(&i);
    if (i != 113626) {
        sdseterr(51,42);
    }
    for (i = 0; i < 3; i++) {
        param[i] = state[i];
    }
    param[3] = time;
    *err = 0;
    *fcnt = 0;
    sdvelfunc(&state[3],param,verrs);
    *fcnt = *fcnt+1;
}

/* This routine performs static analysis. */

void sdstatic(double time,
    double state[6],
    int lock[3],
    double ctol,
    double tol,
    int maxevals,
    int *fcnt,
    int *err)
{
    double resid[3],param[4],jw[9],dw[72],rw[48];
    int iw[24],rooterr,i;

    sdgentime(&i);
    if (i != 113626) {
        sdseterr(52,42);
    }
    for (i = 0; i < 3; i++) {
        param[i] = state[3+i];
    }
    param[3] = time;
    sdroot(sdstatfunc,state,param,3,3,3,lock,
      ctol,tol,maxevals,jw,dw,rw,iw,resid,fcnt,&rooterr);
    sdstatfunc(state,param,resid);
    *fcnt = *fcnt+1;
    if (rooterr == 0) {
        *err = 0;
    } else {
        if (*fcnt >= maxevals) {
            *err = 2;
        } else {
            *err = 1;
        }
    }
}

/* This routine performs steady motion analysis. */

void sdsteady(double time,
    double state[6],
    int lock[6],
    double ctol,
    double tol,
    int maxevals,
    int *fcnt,
    int *err)
{
    double resid[3],param[1];
    double jw[18],dw[162],rw[75];
    int iw[36],rooterr,i;

    sdgentime(&i);
    if (i != 113626) {
        sdseterr(53,42);
    }
    param[0] = time;
    sdroot(sdstdyfunc,state,param,3,6,3,lock,
      ctol,tol,maxevals,jw,dw,rw,iw,resid,fcnt,&rooterr);
    sdstdyfunc(state,param,resid);
    *fcnt = *fcnt+1;
    if (rooterr == 0) {
        *err = 0;
    } else {
        if (*fcnt >= maxevals) {
            *err = 2;
        } else {
            *err = 1;
        }
    }
}

/* This routine performs state integration. */

void sdmotion(double *time,
    double state[6],
    double dstate[6],
    double dt,
    double ctol,
    double tol,
    int *flag,
    int *err)
{
    static double step;
    double work[36],ttime,param[1];
    int vintgerr,which,ferr,i;

    sdgentime(&i);
    if (i != 113626) {
        sdseterr(54,42);
    }
    param[0] = ctol;
    ttime = *time;
    if (*flag != 0) {
        sdmotfunc(ttime,state,dstate,param,&ferr);
        step = dt;
        *flag = 0;
    }
    if (step <= 0.) {
        step = dt;
    }
    sdvinteg(sdmotfunc,&ttime,state,dstate,param,dt,&step,6,tol,work,&vintgerr,&
      which);
    *time = ttime;
    *err = vintgerr;
}

/* This routine performs state integration with a fixed-step integrator. */

void sdfmotion(double *time,
    double state[6],
    double dstate[6],
    double dt,
    double ctol,
    int *flag,
    double *errest,
    int *err)
{
    double work[24],ttime,param[1];
    int ferr,i;

    sdgentime(&i);
    if (i != 113626) {
        sdseterr(55,42);
    }
    param[0] = ctol;
    *err = 0;
    ttime = *time;
    if (*flag != 0) {
        sdmotfunc(ttime,state,dstate,param,&ferr);
        *flag = 0;
    }
    sdfinteg(sdmotfunc,&ttime,state,dstate,param,dt,6,work,errest,&ferr);
    if (ferr != 0) {
        *err = 1;
    }
    *time = ttime;
}
