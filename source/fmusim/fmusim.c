#include "fmusim.h"
#include "fmuio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef _MSC_VER
#define TRUE 1
#define FALSE 0
#define min(a,b) (a>b ? b : a)
#endif

// simulate the given FMU using the forward euler method.
// time events are processed by reducing step size to exactly hit tNext.
// state events are checked and fired only at the end of an Euler step. 
// the simulator may therefore miss state events and fires state events typically too late.
int fmuSimulate(FMU* fmu, double tEnd, double h, fmiBoolean loggingOn, char separator, const char* resultFileName) {
    int i;
    double dt, tPre;
    fmiBoolean timeEvent, stateEvent, stepEvent;
    double simtime;
    int nx;                          // number of state variables
    int nz;                          // number of state event indicators
    double *x;                       // continuous states
    double *xdot;                    // the crresponding derivatives in same order
    double *z = NULL;                // state event indicators
    double *prez = NULL;             // previous values of state event indicators
    fmiEventInfo eventInfo;          // updated by calls to initialize and eventUpdate
    ModelDescription* md;            // handle to the parsed XML file        
    const char* guid;                // global unique id of the fmu
    fmiCallbackFunctions callbacks;  // called by the model during simulation
    fmiComponent c;                  // instance of the fmu 
    fmiStatus fmiFlag;               // return code of the fmu functions
    fmiReal t0 = 0;                  // start time
    fmiBoolean toleranceControlled = fmiFalse;
    int nSteps = 0;
    int nTimeEvents = 0;
    int nStepEvents = 0;
    int nStateEvents = 0;
    FILE* file = 0;
    clock_t timing;

    // instantiate the fmu
    md = fmu->modelDescription;
    guid = getString(md, att_guid);
    callbacks.logger = fmuLogger;
    callbacks.allocateMemory = calloc;
    callbacks.freeMemory = free;
    c = fmu->instantiateModel(getModelIdentifier(md), guid, callbacks, loggingOn);
    if (!c) return fmuError("could not instantiate model");
    
    // allocate memory 
    nx = getNumberOfStates(md);
    nz = getNumberOfEventIndicators(md);
    x    = (double *) calloc(nx, sizeof(double));
    xdot = (double *) calloc(nx, sizeof(double));
    if (nz>0) {
        z    =  (double *) calloc(nz, sizeof(double));
        prez =  (double *) calloc(nz, sizeof(double));
    }
    if (!x || !xdot || ((nz>0) && (!z || !prez))) return fmuError("out of memory");

    // open result file
    if(resultFileName != 0) {
        if((strlen(resultFileName)==1) && (resultFileName[0]=='-') )
            file = stdout;
        else {
            if (!(file=fopen(resultFileName, "w"))) {
                fprintf(stderr,"could not write %s\n", resultFileName);
                return 0; // failure
            }
        }
    }

    timing = clock();
    fprintf(stderr,"timing start %ld\n", timing);

    // set the start time and initialize
    simtime = t0;
    fmiFlag =  fmu->setTime(c, t0);
    if (fmiFlag > fmiWarning) return fmuError("could not set time");
    fmiFlag =  fmu->initialize(c, toleranceControlled, t0, &eventInfo);
    if (fmiFlag > fmiWarning)  fmuError("could not initialize model");
    if (eventInfo.terminateSimulation) {
        fprintf(stderr,"model requested termination at init");
        tEnd = simtime;
    }
  
    // output solution for simtime t0
    if(file) {
        outputRow(fmu, c, t0, file, separator, TRUE);  // output column names
        outputRow(fmu, c, t0, file, separator, FALSE); // output values
    }

    // enter the simulation loop
    while (simtime < tEnd) {
     // get current state and derivatives
     fmiFlag = fmu->getContinuousStates(c, x, nx);
     if (fmiFlag > fmiWarning) return fmuError("could not retrieve states");
     fmiFlag = fmu->getDerivatives(c, xdot, nx);
     if (fmiFlag > fmiWarning) return fmuError("could not retrieve derivatives");

     // advance simtime
     tPre = simtime;
     simtime = min(simtime+h, tEnd);
     timeEvent = eventInfo.upcomingTimeEvent && eventInfo.nextEventTime < simtime;
     if (timeEvent) simtime = eventInfo.nextEventTime;
     dt = simtime - tPre;
     fmiFlag = fmu->setTime(c, simtime);
     if (fmiFlag > fmiWarning) fmuError("could not set time");

     // perform one step
     for (i=0; i<nx; i++) x[i] += dt*xdot[i]; // forward Euler method
     fmiFlag = fmu->setContinuousStates(c, x, nx);
     if (fmiFlag > fmiWarning) return fmuError("could not set states");
     if (loggingOn) fprintf(stderr,"Step %d to t=%.16g\n", nSteps, simtime);
    
     // Check for step event, e.g. dynamic state selection
     fmiFlag = fmu->completedIntegratorStep(c, &stepEvent);
     if (fmiFlag > fmiWarning) return fmuError("could not complete intgrator step");

     // Check for state event
     for (i=0; i<nz; i++) prez[i] = z[i]; 
     fmiFlag = fmu->getEventIndicators(c, z, nz);
     if (fmiFlag > fmiWarning) return fmuError("could not retrieve event indicators");
     stateEvent = FALSE;
     for (i=0; i<nz; i++) 
         stateEvent = stateEvent || (prez[i] * z[i] < 0);  
     
     // handle events
     if (timeEvent || stateEvent || stepEvent) {
        
        if (timeEvent) {
            nTimeEvents++;
            if (loggingOn) fprintf(stderr,"time event at t=%.16g\n", simtime);
        }
        if (stateEvent) {
            nStateEvents++;
            if (loggingOn) for (i=0; i<nz; i++)
                fprintf(stderr,"state event %s z[%d] at t=%.16g\n",
                        (prez[i]>0 && z[i]<0) ? "-\\-" : "-/-", i, simtime);
        }
        if (stepEvent) {
            nStepEvents++;
            if (loggingOn) fprintf(stderr,"step event at t=%.16g\n", simtime);
        }

        // event iteration in one step, ignoring intermediate results
        fmiFlag = fmu->eventUpdate(c, fmiFalse, &eventInfo);
        if (fmiFlag > fmiWarning) return fmuError("could not perform event update");
        
        // terminate simulation, if requested by the model
        if (eventInfo.terminateSimulation) {
            fprintf(stderr,"model requested termination at t=%.16g\n", simtime);
            break; // success
        }

        // check for change of value of states
        if (eventInfo.stateValuesChanged && loggingOn) {
            fprintf(stderr,"state values changed at t=%.16g\n", simtime);
        }
        
        // check for selection of new state variables
        if (eventInfo.stateValueReferencesChanged && loggingOn) {
            fprintf(stderr,"new state variables selected at t=%.16g\n", simtime);
        }
       
     } // if event
     if(file)
        outputRow(fmu, c, simtime, file, separator, FALSE); // output values for this step
     nSteps++;
  } // while  

  timing = clock() - timing;

  // cleanup
  fmu->freeModelInstance(c);
  if(file)
    fclose(file);

  if (x!=NULL) free(x);
  if (xdot!= NULL) free(xdot);
  if (z!= NULL) free(z);
  if (prez!= NULL) free(prez);

  // print simulation summary 
  fprintf(stderr,"Simulation from %g to %g terminated successful\n", t0, tEnd);
  fprintf(stderr,"  steps ............ %d\n", nSteps);
  fprintf(stderr,"  fixed step size .. %g\n", h);
  fprintf(stderr,"  time events ...... %d\n", nTimeEvents);
  fprintf(stderr,"  state events ..... %d\n", nStateEvents);
  fprintf(stderr,"  step events ...... %d\n", nStepEvents);
  if(resultFileName && ((strlen(resultFileName)!=1) || (resultFileName[0] != '-')))
        fprintf(stderr,"CSV file '%s' written.\n", resultFileName);
  fprintf(stderr,"  simultation time.. %g seconds (resolution %g sec) \n",
          ((double)timing * 1.0) / CLOCKS_PER_SEC, 1.0/CLOCKS_PER_SEC);

  return 1; // success
}
