/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: ArgusBG.h,v 1.13 2007/07/12 20:30:49 wouter Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ArgusBG_H
#define ArgusBG_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class ArgusBG : public RooAbsPdf {
 public:
    ArgusBG(){ };
    ArgusBG(const char *name, const char *title,
            RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _c);
    ArgusBG(const char *name, const char *title,
            RooAbsReal& _m, RooAbsReal& _m0,
            RooAbsReal& _c, RooAbsReal& _p);
    ArgusBG(const ArgusBG& other, const char* name = 0);
    virtual TObject* clone(const char* newname) const { return new ArgusBG(*this, newname); }
    inline virtual ~ArgusBG() { }

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
            const char* rangeName = 0) const;
    Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

 protected:
        RooRealProxy m;
        RooRealProxy m0;
        RooRealProxy c;
        RooRealProxy p;

        Double_t evaluate() const;

 private:
        ClassDef(ArgusBG, 1)  // Argus background shape PDF
};

#endif
