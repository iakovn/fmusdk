/* ------------------------------------------------------------------------- 
 * fmusim.h
 * Code for simulating models
 * Copyright 2010 QTronic GmbH. All rights reserved. 
 * -------------------------------------------------------------------------
 */

#ifndef fmusim_h
#define fmusim_h

#include "main.h"

int fmuSimulate(FMU* fmu, double tEnd, double h,
                fmiBoolean loggingOn, char separator, const char* resultFileName);

#endif // fmusim_h
