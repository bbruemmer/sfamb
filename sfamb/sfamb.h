#ifndef SFA_INCLUDED
#define SFA_INCLUDED

#include <database.h>
#include <modelbase.h>
#include <maximize.h>
#include <oxdraw.h>
#include <quadpack.h>

//enum
//{ U_VAR = 3, Z_VAR = 4
//};


enum
{    M_MAXLIK
};

enum {POOLED, WT, LSDV, CFE, TFE}; // 0,1,2,3,4 for 'SetMethod()'

/*-------------------------- Sfa : Modelbase-----------------------------*/
class Sfa : Modelbase
{
    decl m_iModel;                  /* TRUE if model successfully estimated */
    decl m_fPrintDetails;        /* TRUE if output desired after estimation */
    decl m_vStart;
    decl m_iTl;                        /*tl=0: linear
                                ** tl=1:squares & crossprod on all variables
                                ** tl=n:squares & crossprod on first n vars */
    decl m_iResult;                             /* return code from MaxBFGS */
    decl m_cZ;                  /* number of heteroscedasticity variables   */
    decl m_cU;                  /* number of determinants of TE */
    decl m_var;                                   /*  data y~x [cT][1+m_cX] */
    decl m_mU;                  /* determinants of TE */
    decl m_mZ;                     /* Heteroscedasticity, if any [cT][m_cZ] */
    decl m_par;                                       /* parameters: [p][1] */
//    decl m_dLoglik;                                     /* log-likelihood */
    decl m_dLoglikOLS;                                /* log-likelihood OLS */
    decl m_mCovP;                           /*Estimated Variance Covariance */
    decl m_dAlpha;        /*confidence level for bounds on TE, default 0.05 */
    decl m_mScore;   /* [p][cT] Matrix of first derivates evaluated at m_par*/
    decl m_bUseRobustStdErr;
    decl m_bUse_maxSQP;             /* if nonlinear restrictions or parameter
			              bounds are to be considered, the use of
			              MaxSQP is necessary. Algorithm is
			              automatically switched if any of the
			              restriction functions is called      */
    decl  m_bUse_maxSQPF;                 /*alternate "feasible" algorithm 
				 not for admissible for equality constraints*/
    decl m_vUpper;                             /* Upper bounds on parameters*/
    decl m_vLower;                             /* Lower bounds on parameters*/
    decl m_fcon_ge0; /*0: no inequality cons., otherwise func to eval const.*/
    decl m_fcon_eq0;   /*0: no equality cons., otherwise func to eval const.*/
	decl m_bCost;  /*True: Use cost fct type of model,
					i.e., u + v instead of u - v (default) */
					
/*---Additional data members---*/
	decl m_vID;//[#obs][1], firms
	decl m_ID;//[#firms][1], firms
	decl m_vPer;//[#obs][1], period t
	decl m_mMeans;//[#firms][csel], mean values of selected vars
	decl m_mWithins;//[#obs][csel], w-transf. vars
	decl m_Ti;//[#firms][1], individual T-i
	decl m_vTi;//[#obs][1], holds number of individual T-i
	decl m_lnLW;//[#firms][1]vector
	decl m_mTldata;//data for TL form
	decl m_mVars;//only for eq(31)/inactive
	decl m_iDist;//for specification (WT)
	decl m_cParOLS;//for test stat. one-sided err
	
/*---Related to integration (CFE)---*/	
	fKotzetal(const x);
	decl m_lam;
	decl m_rat;
	decl m_vgreps;
	decl m_gti;	

/*---Additional or adjusted functions---*/	
	Ident(const vID, const vPer);
	PrepData(const mSel, iNorm);
	TE();					   //cf. 'TEint()'
	Ineff();
	AiHat();
	Elast(const sXname);
	IDandPer();
	GetLLFi();
	GetTldata();
	GetMeans();
	GetWithins();

	SetConstant();/*used with 'Deterministic' to omit "Constant" in case of the panel models*/
	Select(const iGroup, const aSel);/*to override the conventional 'Select' function
									and deal with "Constant"*/
	SetPrintSfa(const iPri);/*used with 'SetPrint' to print panel information*/
	DropGroupIf(const mifr);
	SetRobustStdErr(const bool);

		//OxPack-related
    SendVarStatus();
    SendSpecials();
//    SendFunctions();
    SendMethods();
    SendDialog(const sDialog);
    SendMenu(const sMenu);
//    SendResults(const sType);
//    ReceiveData();
    ReceiveModel();
    ReceiveDialog(const sDialog, const asOptions, const aValues);
    Sfa();                                              /* constructor */
    GetPackageName();
    GetPackageVersion();
    SetStart(const vStart);
    SetTranslog(const iTl);
    SetConfidenceLevel(const alpha);
    InitData();                                /* extract data from database */
    GridSearch(const beta, const s);                     /* find best gamma */
    DoEstimation(vPar);   /* do the estimation; */
    Covar();
    GetParNames();                      /* overridden for correct labels */
    Output();                   /* calls parent, prints additional info */
    TestGraphicAnalysis();
    GetResults(const par, const eff, const fct, const v); /* get results
                   ** par: coefficients~std.err.~prob
                   ** eff: TE~ lower bound~upper bound (95%)
                   ** fct: likelihood-related ~ variance related (see source)
                   ** v: Variance-Covariance-Matrix m_mCovP */
	GetResiduals(); 			    /* overridden to return the composed error w = v -u  */
	SetPrintDetails(const bool);
//    TE(const dAlpha);          /*Efficiency estimates, cf. GREENE(1992) p82
                       			/* and confidence intervals cf. Horrace&Schmidt(1995)*/
	TEint(const dAlpha);//changed name, cf. 'TE()'
    fSfa(const vP, const adFunc, const avScore, const amHessian);/* LL */

    /* Call one or more of the following functions to estimate subject
       to restrictions or bounds */
    fIneqRest(const fIn); /* Function for inequality restrictions*/
    fEqRest(const fIn);  /*Function for equality restrictions*/
    SetUpperBounds(const vBounds);
    SetLowerBounds(const vBounds);
    SetUse_maxSQPF(const bool);
	SetCost(const bCost);
public:
	/** Types of variables in the model. */
	enum
	{ Y_VAR, X_VAR, U_VAR, Z_VAR
	};
	GetGroupLabels();
	FindGroup(const theGroup); // only required for Ox versions < 8

};
/*------------------------ END Sfa : Database -------------------------*/
#endif /* SFA_INCLUDED */
