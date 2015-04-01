#include <oxstd.h>
#include <oxfloat.h>
#import <modelbase>
#import <maximize>
#import <maxsqp>
#include "sfamb.h"

/*-------------------------- Sfa : Database ---------------------------*/
Sfa::Sfa()
 {
    this->Modelbase();           // intialize base class
    m_iModel = 0;               // 1 when estimated
    m_bUse_maxSQP = m_bUse_maxSQPF = m_fPrintDetails = FALSE;
    m_iResult = -1;
    m_dLoglikOLS = 0;
    m_par = 0;
    m_dAlpha = 0.05;
    m_vStart = 0;
    m_iTl = 0;
    m_bUseRobustStdErr = TRUE;
    m_vUpper = m_vLower = <>;
    m_fcon_ge0 = m_fcon_eq0 = 0;
    if (m_fPrintDetails)
        print("Sfa class, object created on ",date(), ".\n\n");
}

Sfa::Ident(const vID, const vPer)
{switch (m_iMethod){
case 4:
case 3:
case 2:
case 1:
//----------------- Firms ID -----------------------
	decl vtid, i, vid;
	vtid = zeros(sizer(vID),1);//empty matrix
	for (i = 1; i < sizer(vID); ++i)
	{
		vtid[i][] = vID[i][] - vID[i-1][];//to find the first obs. of each firm	(note: value of the first firm is zero!)
	}									  
	vid = selectifr(vtid, vtid .!=0);//get rid of the zeros
	decl vpid = zeros(sizer(vid),1);
	for (i = 0; i < sizer(vid); ++i)
	{
		vpid[i][] = i+2;//"real" ID, i.e. the identifiers 1,2,3...  
	}
	m_ID = 1 | vpid;//[#firms][1]
//--------- ID vector for all observations (expansion)
	decl vFirm = vID[0][] | selectifr(vID, vtid .!= 0);//vector of original firm variable (FirmA, FirmB...)
	decl mData = vFirm ~ m_ID;
	decl vEmp, vRiden = zeros(sizer(vID),1);
	for (i=0; i < sizer(vID); i++)//LoopL[#obs]
	{
	vEmp = selectifr(mData, mData[][0] .== vID[i][]);
	if(sizer(vEmp)>1) {println("\nPanel structure incorrectly specified!",
							   "\nOrder of groups");exit(1);}
	vRiden[i][] = vEmp[][1];
	}
	m_vID = vRiden;//[#obs][1]
	m_vPer = vPer;//orig. data
//---------get individual T-i and expand-------------
	decl mFirm, mTi = zeros(maxc(m_vID), 1);
	for (i = 0; i < maxc(m_vID); i++)//LoopL[#firms] 
	{
	mFirm = selectifr(m_vPer, m_vID .== i+1);
	mTi[i][] = sizer(mFirm);//individual T-i
	if(mTi[i][]<2){
		println("\nPanel structure incorrectly specified!",
				"\nOnly one observation for group ", i+1);exit(1);}
	}
	m_Ti = mTi;
	decl vEm, vAllTi = zeros(sizer(m_vID),1);//expand it to m_vTi
	decl mData2 = m_ID ~ m_Ti;//#firms
	for (i=0; i < sizer(m_vID); i++)//LoopL[#obs], expansion
	{
	vEm = selectifr(mData2, mData2[][0] .== m_vID[i][]);
	vAllTi[i][] = vEm[][1];
	}
	m_vTi = vAllTi;//T-i as a [#obs][1] vector
break;
default:
	return FALSE;
break;
}}

Sfa::DropGroupIf(const mifr)
{
	decl mSel = selectifr(IDandPer(), mifr);

	decl i, vDiff, vDrop = zeros(sizer(mSel),1);

	for (i = 1; i < sizer(mSel); i++)
			{vDiff = mSel[i][0] - mSel[i-1][0];

			if(vDiff .!= 0)
			{vDrop[i][] = mSel[i][0];}}//find IDs

	vDrop = mSel[0][0] | vDrop[1:][];//append first ID

	vDrop = deleteifr(vDrop, vDrop .== 0);//delete zeros
	
	//expand:
	decl vIndex = zeros(sizer(m_vID), 1);
	
	for (i = 0; i < sizer(m_vID); i++){
		vIndex[i][] = selectifr(vDrop, vDrop .== m_vID[i][]);}

	//shorten the vectors of identification
	decl mvidnew = deleteifr(m_vID, vIndex .!= 0);
	decl mvpernew = deleteifr(m_vPer, vIndex .!= 0);

	//shorten the database
	m_mData = deleteifr(m_mData, vIndex .!= 0);

	println("\nRestricted Sample:");
	
	Ident(mvidnew, mvpernew);//new call to 'Ident()' to renew the corresp. data members
}

Sfa::SetPrintSfa(const iPri)
{switch (m_iMethod){
case 4:
case 3:
case 2:
case 1:
	print("\n\t  #groups:   #periods(max):  avg.T-i:");
	print(maxc(m_vID) ~ maxc(m_vPer) - minc(m_vPer) + 1 ~ sumc(m_Ti)/maxc(m_ID));

	SetPrint(iPri);
break;

default:
	SetPrint(iPri);
break;
}}

Sfa::IDandPer()
{switch (m_iMethod){
case 4:
case 3:
case 2:
case 1:
	return m_vID ~ m_vTi;//Firm number ~ individual T-i; [#obs][2]
break;
default:
	return FALSE;
break;
}}

Sfa::PrepData(const mSel, iNorm)
{
	decl mVar;
	
	if (iNorm == TRUE){
		mVar = log(mSel/meanc(mSel));}
	
	else {mVar = log(mSel);}

	return mVar;
}

Sfa::SendVarStatus(){
    return {
        {"&Y variable", 'Y', STATUS_GROUP + STATUS_ENDOGENOUS, Y_VAR},
        {"&X variable", 'X', STATUS_GROUP, X_VAR},
        {"&Technical effect", 'T', STATUS_GROUP2, U_VAR},
        {"&Heteroscedastic", 'H', STATUS_GROUP2, Z_VAR}
    };
}

Sfa::SendSpecials(){
    return {"Constant", "Trend"};
}

Sfa::SendMethods(){
    return {{"Maximum Likelihood", M_MAXLIK, FALSE, 0}};
}

// Sfa::SendDialog(const sDialog){
//     if (sDialog == "OP_SETTINGS")
//     {
//         return
//         {
//             { "Translog (0:no, 1:all X-vars, n>1:first n X-Vars)" },
//             { "value:", CTL_INT, m_iTl, "tl"},
//             { "Confidence level for efficiency bounds" },
//             { "significance:", CTL_DOUBLE,   m_dAlpha,   "alpha"},
//             { "Estimation Options"},
//             { "Use user supplied start values", CTL_STRMAT, m_vStart, "start"},
//             { "Use robust standard errors", CTL_CHECK, m_bUseRobustStdErr, "robstderr"},
//             { "Output Options"},
//             { "Detailed Output", CTL_CHECK, m_fPrintDetails, "printdetails"}
//         };
//     }
//     // allow base class to process unhandled cases
//     else return Modelbase::SendDialog(sDialog);
// }

Sfa::SendMenu(const sMenu){
    if (sMenu == "Test")
    {   return
            {{ "&Graphic Analysis", "OP_TEST_GRAPHICS"},
             0,
             { "&Exclusion Restrictions...", "OP_TEST_SUBSET"},
             { "&Linear Restrictions...",    "OP_TEST_LINRES"},
             0,
             { "Store in D&atabase...", "OP_TEST_STORE"}
             };
    }
}

Sfa::ReceiveModel(){
    decl as, i;

    DeSelect();

    // get selection of database variables
    Select(Y_VAR, "OxPackGetData"("SelGroup", Y_VAR));
    Select(X_VAR, "OxPackGetData"("SelGroup", X_VAR));
    Select(U_VAR, "OxPackGetData"("SelGroup", U_VAR));
    Select(Z_VAR, "OxPackGetData"("SelGroup", Z_VAR));

    decl freq, year1, period1, year2, period2;
    [year1, period1, year2, period2] = "OxPackGetData"("SelSample");
    ForceSelSample(year1, period1, year2, period2);

    // get functions
    // get method
    decl imethod, itforc, brecursive;
    [imethod, itforc, brecursive] = "OxPackGetData"("Method");
    SetMethod(imethod);
}

// Sfa::ReceiveDialog(const sDialog, const asOptions, const aValues){
//     if (sDialog == "OP_SETTINGS")
//     {
//             m_iTl       = aValues[0];
//             m_dAlpha    = aValues[1];
//             m_vStart    = aValues[2];
//             m_bUseRobustStdErr = aValues[3];
//             m_fPrintDetails = aValues[4];
// 			return 1;
//     }
//     else if (sDialog == "OP_TEST_GRAPHICS")
//         TestGraphicAnalysis();
//     else if (sDialog == "OP_TEST_STORE")
//     {
//          decl asdlg, asoptions, avalues, mEffs;
//          if ("OxPackDialog"("Store in Database",
//                   {       { "Store in database" },
//                           { "TE (point estimates)", CTL_CHECK, 1, "effs" },
//                           { "TE (lower bound)", CTL_CHECK, 1, "lower" },
//                           { "TE (upper bound)", CTL_CHECK, 1, "upper" }
//                   },
//                   &asoptions, &avalues)
//             )
//          {

//              mEffs = TE(m_dAlpha);
//              if (avalues[0] > 0)
//                      "OxPackStore"(mEffs[][0], 0, GetSize() - 1, "effs");
//              if (avalues[1] > 0)
//                      "OxPackStore"(mEffs[][1], 0, GetSize() - 1, "lower");
//              if (avalues[2] > 0)
//                      "OxPackStore"(mEffs[][2], 0, GetSize() - 1, "upper");
//          }
//     }
//     else
//     {       // allow base class to process unhandled cases
//             Modelbase::ReceiveDialog(sDialog, asOptions, aValues);
//     }
// }

Sfa::GetPackageName(){return "Sfa";}

Sfa::GetPackageVersion(){return "1.0";}

Sfa::SetConstant()
{switch (m_iMethod){
case 4:
case 3:
case 2:
case 1:
	return FALSE;
break;

default:
	Deterministic(-1);// (-1) == no seasonals
break;
}}

Sfa::Select(const iGroup, const aSel)
{switch (m_iMethod){
case 4:
case 3:
case 2:
case 1:
	decl aSelnew;

	if (iGroup != Y_VAR)
	{
		decl iIndex = strfind(aSel, "Constant");//relevant index of "Constant"

		if (iIndex == -1) //in case "Constant" is not specified
		{
		aSelnew = aSel;
		}
		else
		{
		aSelnew = dropr(aSel, iIndex~iIndex+1~iIndex+2);//drop "Constant" and related zeros
		}
	}
	else
	{
		aSelnew = aSel;
	}

	decl asVars = Database::Select(iGroup, aSelnew);

	return asVars;
break;

default:
	asVars = Database::Select(iGroup, aSel);

	return asVars;
break;
}}

Sfa::SetStart(const vStart)
{switch (m_iMethod){
case 4:
        m_vStart = vStart | zeros(sizer(m_ID),1);
        return TRUE;
break;

default:
        m_vStart = vStart;
        return TRUE;
break;		
}}

Sfa::SetTranslog(const iTl){
        m_iTl = iTl;
        return TRUE;
}

Sfa::SetConfidenceLevel(const alpha){
        m_dAlpha = alpha;
        return TRUE;
}

Sfa::SetLowerBounds(const vBounds){
  m_bUse_maxSQP=TRUE;
  if (m_vUpper != <>){
    if (m_vUpper < vBounds){
      eprint("At least one upper bound is below the lower bound\n");
      eprint("no lower bounds have been applied");
      return FALSE;
    }
  }
  m_vLower = vBounds;
  return FALSE;
}

Sfa::SetUpperBounds(const vBounds){
  m_bUse_maxSQP=TRUE;
  if (m_vLower != <>){
    if (m_vLower > vBounds){
      eprint("At least one lower bound is below the lower bound\n");
      eprint("no upper bounds have been applied");
      return FALSE;
    }
  }
  m_vUpper = vBounds;
  return FALSE;
}

Sfa::fEqRest(const fIn){
  m_bUse_maxSQP=TRUE;
  m_fcon_eq0=fIn;
  return TRUE;
}

Sfa::fIneqRest(const fIn){
  m_bUse_maxSQP=TRUE;
  m_fcon_ge0=fIn;
  return TRUE;
}

Sfa::fKotzetal(const x)
{
	return
			densn(x)
			*prodc(
			probn(
			m_vgreps * (-m_rat) - sqrt(m_lam^2/m_gti)*x));
}

Sfa::fSfa(const vP, const prob, const avScore, const amHessian)
{switch (m_iMethod){
case 4:
	decl epsi  = m_var[][0]-m_var[][1:m_cX]*vP[0:(m_cX-1)];

	decl svsq = exp(vP[m_cX]);//sig-v^2
	decl susq = exp(vP[m_cX+1]);//sig-u^2

	decl lam = sqrt(susq) / sqrt(svsq);//lambda

	decl sigsq = susq + svsq;//sigma^2

	prob[0] = double(meanc(
    -.9189385-0.5.*log(sigsq)- (epsi.^2) ./(2*sigsq) + log(probn(-lam.*epsi./sqrt(sigsq)))
    -log(probn(0))));

	m_mScore = zeros(rows(m_var),(m_cX+1));

	return (!(isnan(prob[0])||isdotinf(prob[0]) ));      // 1 indicates success
break;

case 3:

	epsi  = m_var[][0]-m_var[][1:m_cX]*vP[0:(m_cX-1)];

	svsq = exp(vP[m_cX]);//sig-v^2
	susq = exp(vP[m_cX+1]);//sig-u^2

	lam = sqrt(susq) / sqrt(svsq);//lambda

	sigsq = susq + svsq;//sigma^2

	decl rat = lam / sqrt(sigsq);//  lambda / sigma
	
	m_lam = lam;
	m_rat = rat;

	decl i, veps, epsistar, mS1, lnfx, lnFX, result, abserr;
	lnfx = lnFX = zeros(sizer(m_ID),1);

	oxwarning(FALSE);//warnings off (maximization)

	if (m_fPrintDetails) 
		oxwarning(-1);
	
	for (i = 0; i < sizer(m_ID); i++)//loopL[#firms]
	{
	veps = selectifr(epsi, m_vID .== i+1);	
	m_gti = m_Ti[i][];
	m_vgreps = veps;
	epsistar = veps[:m_Ti[i][]-2][];

	mS1 = 	unit(m_Ti[i][]-1,m_Ti[i][]-1)
			+ constant((-1/m_Ti[i][]),m_Ti[i][]-1,m_Ti[i][]-1);

	lnfx[i][] =
				(-(m_Ti[i][]-1)/2) * log(2*M_PI)

				-0.5*log(determinant(sigsq * mS1))

				-0.5 * epsistar' * invertgen(sigsq * mS1) * epsistar;//pdf

	QPWARN(FALSE);
	QNG(fKotzetal, -10, 10, &result, &abserr);//cdf
	lnFX[i][] = log(result);
	}								   

	decl iC = -meanc(m_Ti) * log(probn(0)) * sizer(m_ID);

	prob[0] = double( iC + sumc(lnfx) + sumc(lnFX) );
	
	
	m_lnLW = iC/sizer(m_ID)  + lnfx + lnFX;

	m_mScore = zeros(rows(m_var),(m_cX+1));

	return (!(isnan(prob[0])||isdotinf(prob[0]) ));      // 1 indicates success
break;

case 1:

	epsi  = m_var[][0]-m_var[][1:m_cX]*vP[0:(m_cX-1)];

decl hit = exp(m_mZ*vP[m_cX+2:m_cX+m_cZ+1]);

	svsq = exp(vP[m_cX]);
	susq = exp(vP[m_cX+1]);

decl mu = 0;
if (m_iDist)			  
	 mu = vP[m_cX+m_cZ+2];

	oxwarning(FALSE);//warnings off (maximization)

	if (m_fPrintDetails) 
		oxwarning(-1);

/*------ w-transf. of hit-equation ------------------------------------------*/
	decl vFm, vFw; 
	decl mMh = zeros(maxc(m_vID), 1);
	decl mMhAll = zeros(sizer(m_vID), 1); 
	decl mWh = zeros(sizer(m_vID), 1); 
//---indiv. mean
	for (i = 0; i < maxc(m_vID); i++)//loopL[#firms] 
	{
	vFm = selectifr(hit, m_vID .== i+1);//firm-specific vector of h-it
	mMh[i][] = (meanc(vFm));//Mean of h-it
	}
//---expansion and w-trans
	decl mMhSel = m_ID ~ mMh;//identifiers + acc means, [#firms][2];
	for (i = 0; i < sizer(m_vID); i++)//loopL[#obs]
	{
	vFw = selectifr(mMhSel, mMhSel[][0] .== m_vID[i][]);//firm data
	mMhAll[i][] = vFw[][1];//expansion
	mWh[i][] = hit[i][] - mMhAll[i][];//w-transf., for use in evaluator => vector (mWh), [#obs][1]
	}			   
/*------ w-t-hit------------------------------------------*/
	decl mFirm, mData = m_vID ~ m_vTi ~ epsi ~ mWh;
/*------ definition of SIGMA(inverse) and aux matrices------------------------------------------*/
	decl s1ti, s2ti, Sig, Siginv;
	decl m1 = zeros(sizer(m_vID),1);
	decl m2 = zeros(sizer(m_vID),1);
	decl m3 = zeros(sizer(m_vID),1);

	for (i = 0; i < sizer(m_vID); i++)//loopL[#obs]
	{								   
	s1ti = s2ti = zeros(m_vTi[i][],m_vTi[i][]);//empty
	//values
	s1ti = constant((-1/m_vTi[i][])*svsq, m_vTi[i][], m_vTi[i][]);
	s2ti[][] = unit(m_vTi[i][])*svsq;
	//Sig = s1ti + s2ti;//the inverse is of interest
	Siginv = invertgen(s1ti + s2ti);
	mFirm = selectifr(mData, mData[i][0] .== m_vID);//indiv. firm data: [2]=eps, [3]=hit
	//M1, M2, M3
	m1[i][] = mFirm[][2]' * Siginv * mFirm[][3];
	m2[i][] = mFirm[][3]' * Siginv * mFirm[][3];
	m3[i][] = mFirm[][2]' * Siginv * mFirm[][2];
	}
//---mustar, sigmastar, LLF---
	//Extract values to get #firms result
	decl im1 = zeros(sizer(m_ID),1);
	decl im2 = zeros(sizer(m_ID),1);
	decl im3 = zeros(sizer(m_ID),1);
	//for results:
	decl mustar = zeros(sizer(m_ID),1);
	decl sigstar = zeros(sizer(m_ID),1);
	decl lnLW = zeros(sizer(m_ID),1);

	for (i = 0; i < sizer(m_ID); i++)//loopL[#firms]
	{
	im1[i][] = selectifr(m1, m_vID .== i+1);
	im2[i][] = selectifr(m2, m_vID .== i+1);
	im3[i][] = selectifr(m3, m_vID .== i+1);

	mustar[i][]= (mu/susq - im1[i][]) / (im2[i][] + 1/susq);//eq(27)

	sigstar[i][]= 1/ (im2[i][] + 1/susq);//eq(28)

	lnLW[i][] =	-(m_Ti[i][]-1)*log(sqrt(svsq)) -0.5*(m_Ti[i][]-1)*log(2*M_PI)
	-0.5*im3[i][] +0.5*(mustar[i][]^2/sigstar[i][])
	-0.5*((mu/sqrt(susq))^2) + 0.5*log(sigstar[i][])
	-0.5*log(sqrt(susq)^2) + log(probn(mustar[i][]/sqrt(sigstar[i][])))
	-log(probn(mu/sqrt(susq)));
	}//LLF, eq(26)

	prob[0] = double(sumc(lnLW));

	m_mScore = zeros(rows(m_var),(m_cX+m_cZ+1));
	m_lnLW = lnLW;//[#firms][1]vector with individual LLF value

	return (!(isnan(prob[0])||isdotinf(prob[0]) ));      // 1 indicates success
break;

default:

	decl my=0;
    decl u  = m_var[][0]-m_var[][1:m_cX]*vP[0:(m_cX-1)];
    if (m_cU>0)
        my = -m_mU[][]*vP[m_cX+m_cZ+1:m_cX+m_cZ+m_cU];
    decl uu = (u-my).*(u-my);
    decl s  = sqrt(exp(2*vP[m_cX])
        +exp(2*m_mZ*vP[m_cX+1:m_cX+m_cZ]));
    decl l = exp(m_mZ*vP[m_cX+1:m_cX+m_cZ]-vP[m_cX]);
    decl ra = my./(s.*l);
    prob[0] = double(meanc(
    -.9189385-0.5.*log(s.^2)-uu./(2*s.^2) + log(probn(-l.*u./s-ra))
    -log(probn(-ra.*sqrt(1+l.^2)))));

if (avScore)
{
    decl grad = zeros(rows(m_var),(m_cX+m_cZ+m_cU+1));
    decl sv2 =exp(2*vP[m_cX]);
    decl su2=exp(2*m_mZ*vP[m_cX+1:m_cX+m_cZ]);
    decl darg = probn(-l.*u./s-ra);
    decl varg = densn(-l.*u./s-ra);
    decl wurz = sqrt(1+(l.^2));
    decl darh = probn(-ra.*wurz);
    decl varh = densn(-ra.*wurz);
    decl tmp  = (u-my)./(s.^2)+l./s.*varg./darg;
    grad[][0:m_cX-1]=tmp.*m_var[][1:m_cX];         // diff nach beta
    tmp = (-su2./s.^2)+uu.*su2./(s.^4)
       +(varg.*((l.*u./s).*(-1+(su2./s.^2)) + ra.*((su2./s.^2)+1))./darg)
       -(varh./darh).*(wurz.*su2./s.^2 + wurz - (l.^2)./wurz).*ra;
    grad[][m_cX+1:m_cX+m_cZ]=tmp.*m_mZ; // diff nach delta
    tmp = -(u-my)./(s.^2) + varg./(s.*l.*darg) - varh.*wurz./(s.*l.*darh);
    if (m_cU>0)
        grad[][m_cX+m_cZ+1:m_cX+m_cZ+m_cU]=tmp.*m_mU; //diff nach gamma
    tmp= (sv2./s.^2).*(-1+uu./s.^2)
      +(varg.*((l.*u./s).*(1+sv2./s.^2)+ra.*((sv2./s.^2)-1))./darg)
      -ra.*varh.*(wurz.*sv2./s.^2-wurz+(l.^2)./wurz)./darh;
    grad[][m_cX] =tmp;
    m_mScore = grad';
    avScore[0] = (meanc(grad))';
}

return (!(isnan(prob[0])||isdotinf(prob[0]) ));      // 1 indicates success
break;
}}

Sfa::TEint(const dAlpha)
{switch (m_iMethod){
case 4:
case 3:
case 2:
case 1:
	println("\n'TEint()' is available for the POOLED model only \n"); return FALSE;
break;

default:
  	decl my=0;
    decl u  = m_var[][0]-m_var[][1:m_cX]*m_par[0:(m_cX-1)];
    if (m_cU>0)
        my = -m_mU[][]*m_par[m_cX+m_cZ+1:m_cX+m_cZ+m_cU];
    decl uu = (u-my).*(u-my);
    decl s  = sqrt(exp(2*m_par[m_cX])
        +exp(2*m_mZ*m_par[m_cX+1:m_cX+m_cZ]));
    decl l = exp(m_mZ*m_par[m_cX+1:m_cX+m_cZ]-m_par[m_cX]);
    decl gam = 1 ./ (1 + l .^ 2); //attention: gam is s_v^2/s, not s_u^2/s!!
    decl su2 = (l .^ 2) .* gam .* (s .^ 2);
    decl sv2 = su2 ./ l .^ 2;
    decl mystar = gam .*  - my + (1 - gam) .* ( - u);
    decl s2star = gam .* su2;
                //here the calculation of confidence bounds starts
    decl zlow=quann(1-(dAlpha/2)*(1-probn(-mystar./sqrt(s2star))));
    decl zupper=quann(1-(1-dAlpha/2)*(1-probn(-mystar./sqrt(s2star))));
//	print(zlow~zupper);exit(0);
    decl lower=exp(-mystar-zlow.*sqrt(s2star));
    decl upper=exp(-mystar-zupper.*sqrt(s2star));
    return( (probn(mystar ./ sqrt(s2star) - sqrt(s2star))
         ./ probn(mystar ./ sqrt(s2star)) .* exp( - mystar + 0.5 .* s2star))
         ~lower~upper);
	//return (!(isnan(eff[0])||isdotinf(eff[0]) ));      // 1 indicates success
break;
}}

Sfa::TE()
{switch (m_iMethod){
case 4:
case 3:
	decl my=0;
	decl u  = GetResiduals();
	decl svsq = exp(m_par[m_cX]);
	decl susq = exp(m_par[m_cX+1]);
	decl s  = sqrt(svsq + susq);
	decl lam =  sqrt(susq) / sqrt(svsq);
	decl gam = 1 ./ (1 + lam .^ 2);

    decl mystar = gam .*  - my + (1 - gam) .* ( - u);
    decl s2star = gam .* susq;

	return( (probn(mystar ./ sqrt(s2star) - sqrt(s2star))
         ./ probn(mystar ./ sqrt(s2star)) .* exp( - mystar + 0.5 .* s2star)));
break;

case 2:
	decl valpha = AiHat();
	decl azero = maxc(valpha);//normalization, S&S1984,p.368
	decl vineff = azero - valpha;//ineffs
//	decl veff = exp(-vineff);//TE est.
	decl mStat = m_ID~exp(-vineff);//data: identifiers~TE scores(cf.K&L2000,p98,eq3.3.5)
	decl i, iSel, vtech = zeros(sizer(m_vID),1);//expand to [#obs]
	for (i=0; i < sizer(m_vID); i++)//LoopL[#obs]
	{
	iSel = selectifr(mStat, mStat[][0] .== m_vID[i][]);//Selection
	vtech[i][] = iSel[][1];//expansion
	}
	return vtech;	
//	return exp(-vineff);//TE est.
break;

case 1:
	decl mu=0;

	if (m_iDist)
	 mu = m_par[m_cX+m_cZ+2];

	decl vhit = exp(m_mZ*m_par[m_cX+2:m_cX+m_cZ+1]);
	decl vepsi  = 0-m_var[][1:m_cX]*m_par[0:(m_cX-1)];//eps-i-tilde-hat, see eq(30)
		//in stata:  gen `res_m'= `y'-`fun2'    /* residuals */
/*---Aux. matrices (similar to 'fWtm')---*/
//	decl i, s1ti, s2ti, Sig, Siginv, mFirm;
	decl s1ti, s2ti, Sig, Siginv, mFirm;
	decl m1 = zeros(sizer(m_vID),1);
	decl m2 = zeros(sizer(m_vID),1);//
	decl mustar = zeros(sizer(m_vID),1);//now [#obs!]
	decl sigstar = zeros(sizer(m_vID),1);

	decl mData = m_vID ~ m_vTi ~ vepsi ~ vhit;

	for (i = 0; i < sizer(m_vID); i++)//loopL[#obs]
	{								   
	s1ti = s2ti = zeros(m_vTi[i][],m_vTi[i][]);//components of SIGMA
	//values
	s1ti = constant((-1/m_vTi[i][])*exp(m_par[m_cX]), m_vTi[i][], m_vTi[i][]);
	s2ti[][] = unit(m_vTi[i][])*exp(m_par[m_cX]);
	//Sig = s1ti + s2ti;
	Siginv = invertgen(s1ti + s2ti);
	mFirm = selectifr(mData, mData[i][0] .== m_vID);//indiv. firm data: [2]=eps, [3]=hit
	//M1, M2
	m1[i][] = mFirm[][2]' * Siginv * mFirm[][3];//epsi*inv*hit
	m2[i][] = mFirm[][3]' * Siginv * mFirm[][3];//hit*inv*hit
	mustar[i][]= (mu/exp(m_par[m_cX+1]) - m1[i][]) / (m2[i][] + 1/exp(m_par[m_cX+1]));//
	sigstar[i][]= 1/ (m2[i][] + 1/exp(m_par[m_cX+1]));
	}
//------------jlms index-------------------
	decl r1 = densn(mustar ./ sqrt(sigstar));
	decl r2 = probn(mustar ./ sqrt(sigstar));
	decl vjlms = vhit .* (mustar + sqrt(sigstar) .* r1 ./ r2);//eq(30)
//------------BC index-------------------
	decl r3 = probn(mustar ./ sqrt(sigstar) - vhit .* sqrt(sigstar));//first PHI of BC-index
	decl bcind = (r3 ./ r2) .* exp(-vhit .* mustar + 0.5 * vhit.^2 .* sigstar); 
//	return vjlms~bcind;
	return bcind;
break;

default:
	my=0;
    u  = m_var[][0]-m_var[][1:m_cX]*m_par[0:(m_cX-1)];
    if (m_cU>0)
        my = -m_mU[][]*m_par[m_cX+m_cZ+1:m_cX+m_cZ+m_cU];
    decl uu = (u-my).*(u-my);
    s  = sqrt(exp(2*m_par[m_cX])
        +exp(2*m_mZ*m_par[m_cX+1:m_cX+m_cZ]));
    decl l = exp(m_mZ*m_par[m_cX+1:m_cX+m_cZ]-m_par[m_cX]);
    gam = 1 ./ (1 + l .^ 2); //attention: gam is s_v^2/s, not s_u^2/s!!
    decl su2 = (l .^ 2) .* gam .* (s .^ 2);
    decl sv2 = su2 ./ l .^ 2;
    mystar = gam .*  - my + (1 - gam) .* ( - u);
    s2star = gam .* su2;
                //here the calculation of confidence bounds starts =>shifted to 'TEint()'
//    decl zlow=quann(1-(dAlpha/2)*(1-probn(-mystar./sqrt(s2star))));
//    decl zupper=quann(1-(1-dAlpha/2)*(1-probn(-mystar./sqrt(s2star))));
////	print(zlow~zupper);exit(0);
//    decl lower=exp(-mystar-zlow.*sqrt(s2star));
//    decl upper=exp(-mystar-zupper.*sqrt(s2star));
//    return( (probn(mystar ./ sqrt(s2star) - sqrt(s2star))
//         ./ probn(mystar ./ sqrt(s2star)) .* exp( - mystar + 0.5 .* s2star))
//         ~lower~upper);
    return( (probn(mystar ./ sqrt(s2star) - sqrt(s2star))
         ./ probn(mystar ./ sqrt(s2star)) .* exp( - mystar + 0.5 .* s2star)));
break;
}}

Sfa::Ineff()
{switch (m_iMethod){
case 4:
case 3:
	decl my=0;
	decl u  = GetResiduals();//eps
	decl svsq = exp(m_par[m_cX]);
	decl susq = exp(m_par[m_cX+1]);
	decl s  = sqrt(svsq + susq);//sigma
	decl lam =  sqrt(susq) / sqrt(svsq);
	decl gam = 1 ./ (1 + lam .^ 2); //attention: gam is s_v^2/s, not s_u^2/s!!
	decl mystar = gam .*  - my + (1 - gam) .* ( - u);//cf. mustar of B&C1988,eq(9)
	decl s2star = gam .* susq;//sigma-star^2 = sigma_u^2*sigma_v^2/sigma^2; cf. B&C1988,eq(10)

	decl vjlms = sqrt(s2star) .* (				  //cf. K&L2000,eq(3.2.50)
								  mystar./sqrt(s2star)+
								  densn(mystar./(sqrt(s2star))) ./ (1-probn(-mystar./(sqrt(s2star))))
								  );
	return vjlms;
break;

case 2:
	decl valpha = AiHat();
	decl azero = maxc(valpha);//normalization, S&S1984,p.368
	decl vui = azero - valpha;//ineffs

	decl i, vineff = zeros(sizer(m_vID),1);
	for (i=0; i < sizer(m_vID); i++)//LoopL[#obs]
	{
	vineff[i][] = selectifr(vui, m_ID .== m_vID[i][]);//Selection
	}
	return vineff;
break;

case 1:
	decl mu=0;

	if (m_iDist)
	 mu = m_par[m_cX+m_cZ+2];

	decl vhit = exp(m_mZ*m_par[m_cX+2:m_cX+m_cZ+1]);
	decl vepsi  = 0-m_var[][1:m_cX]*m_par[0:(m_cX-1)];//eps-i-tilde-hat, see eq(30)
		//in stata:  gen `res_m'= `y'-`fun2'    /* residuals */
/*---Aux. matrices (similar to 'fWtm')---*/
//	decl i, s1ti, s2ti, Sig, Siginv, mFirm;
	decl s1ti, s2ti, Sig, Siginv, mFirm;
	decl m1 = zeros(sizer(m_vID),1);
	decl m2 = zeros(sizer(m_vID),1);//
	decl mustar = zeros(sizer(m_vID),1);//now [#obs!]
	decl sigstar = zeros(sizer(m_vID),1);

	decl mData = m_vID ~ m_vTi ~ vepsi ~ vhit;

	for (i = 0; i < sizer(m_vID); i++)//loopL[#obs]
	{								   
	s1ti = s2ti = zeros(m_vTi[i][],m_vTi[i][]);//components of SIGMA
	//values
	s1ti = constant((-1/m_vTi[i][])*exp(m_par[m_cX]), m_vTi[i][], m_vTi[i][]);
	s2ti[][] = unit(m_vTi[i][])*exp(m_par[m_cX]);
	//Sig = s1ti + s2ti;
	Siginv = invertgen(s1ti + s2ti);
	mFirm = selectifr(mData, mData[i][0] .== m_vID);//indiv. firm data: [2]=eps, [3]=hit
	//M1, M2
	m1[i][] = mFirm[][2]' * Siginv * mFirm[][3];//epsi*inv*hit
	m2[i][] = mFirm[][3]' * Siginv * mFirm[][3];//hit*inv*hit
	mustar[i][]= (mu/exp(m_par[m_cX+1]) - m1[i][]) / (m2[i][] + 1/exp(m_par[m_cX+1]));//
	sigstar[i][]= 1/ (m2[i][] + 1/exp(m_par[m_cX+1]));
	}
//------------jlms index-------------------
	decl r1 = densn(mustar ./ sqrt(sigstar));
	decl r2 = probn(mustar ./ sqrt(sigstar));
	vjlms = vhit .* (mustar + sqrt(sigstar) .* r1 ./ r2);//eq(30)
	return vjlms;//[#obs][1]
break;

default:
	my=0;
    u  = m_var[][0]-m_var[][1:m_cX]*m_par[0:(m_cX-1)];//eps
    if (m_cU>0)
        my = -m_mU[][]*m_par[m_cX+m_cZ+1:m_cX+m_cZ+m_cU];
    decl uu = (u-my).*(u-my);//eps^2
    s  = sqrt(exp(2*m_par[m_cX])
        +exp(2*m_mZ*m_par[m_cX+1:m_cX+m_cZ]));//sigma = sqrt(sigma_v^2+sigma_u^2)
    decl l = exp(m_mZ*m_par[m_cX+1:m_cX+m_cZ]-m_par[m_cX]);//lambda
    gam = 1 ./ (1 + l .^ 2); //attention: gam is s_v^2/s, not s_u^2/s!!
    decl su2 = (l .^ 2) .* gam .* (s .^ 2);//sigma_u^2
    decl sv2 = su2 ./ l .^ 2;//sigma_v^2
    mystar = gam .*  - my + (1 - gam) .* ( - u);//cf. mustar of B&C1988,eq(9)
    s2star = gam .* su2;//sigma-star^2 = sigma_u^2*sigma_v^2/sigma^2; cf. B&C1988,eq(10)
	vjlms = sqrt(s2star) .* 	  (				  //cf. K&L2000,eq(3.2.50)
								  mystar./sqrt(s2star)+
								  densn(mystar./(sqrt(s2star))) ./ (1-probn(-mystar./(sqrt(s2star))))
								  );
	return vjlms;
break;
}}

Sfa::GetResiduals()		
{switch (m_iMethod){
case 3:
case 2:
case 1:
	decl my = GetGroup(Y_VAR);
	decl mx = GetGroup(X_VAR);	
	decl valpha = AiHat();

	decl i,vaexp = zeros(sizer(m_vID),1);

	for (i=0; i < sizer(m_vID); i++)//LoopL[#obs], expand to NTx1
	{
		vaexp[i][] = selectifr(valpha, m_ID .== m_vID[i][]);
	}

	if (m_iTl)
		return my - m_mTldata[][1:] * m_par[0:(m_cX-1)] - vaexp;
	
	return my - mx * m_par[0:(m_cX-1)] - vaexp; //returns eps_it = y_it - y^hat = y_it - X_it'beta - a_i (CFE: a_i^M)
break;

default:			  
	return m_var[][0]-m_var[][1:m_cX]*m_par[0:(m_cX-1)];// returns w= v-u = y - \hat y
break;	
}}

Sfa::InitData()
{switch (m_iMethod){
case 4:
    decl  cp, vp, mh, i;

    m_iT1est = m_iT1sel;  m_iT2est = m_iT2sel;

	decl my = GetGroup(Y_VAR);
	decl mx = GetGroup(X_VAR);	
    //m_cX = columns(mx);
	decl cX = columns(mx);
//-------------------SetTranslog-----------------------------
	if (m_iTl)
    {   print("Constructing Squares and Cross-Products...");
        if (m_iTl>1) cX=m_iTl;
	mx=mx~0.5*(mx[][0:(cX-1)] .* mx[][0:(cX-1)]);//squares

	for (i = 0 ; i < cX-1 ; i++)
    {										
		mx=mx~(mx[][i+1:cX-1] .* mx[][i]);//cross-products
	}
        print("done.\n");

	m_mTldata = my~mx;//constructed data
	}
	
	m_cT = m_iT2est - m_iT1est + 1;
    if (m_cT <= 2)
    {   eprint("Only ", m_cT, " observations. This is not enough.\n");
        return FALSE;
    }
//---------------------dummies-----------------
	decl j, mDum = zeros(sizer(m_vID),sizer(m_ID));//NT x N
	
	for (i = 0; i < sizer(m_vID); i++){//[NT]rows 
		for (j = 0; j < sizer(m_ID); j++)//[N]columns
		{
			if   (m_vID[i][] == m_ID[j][])
			mDum[i][j]=mDum[i][j]+1;						
		}}

	m_var = my~mx~mDum;
	m_cX = columns(mx)+columns(mDum);
	
	m_mY = my;
	m_cY = 1;

	decl mZ = GetGroup(Z_VAR);
    decl cZ = columns(mZ);

    decl mU = GetGroup(U_VAR);
    decl cU = columns(mU);

		if (cZ || cU){
			print("TFE model: neither Z-VAR nor U-VAR \n");
			return FALSE;
		}

	if (!m_cX)
    {   eprint("Need some regressors\n");
        return FALSE;
    }
	
    m_iModelStatus = MS_DATA;

    return TRUE;
break;

case 3:
case 2:
case 1:
    m_iT1est = m_iT1sel;  m_iT2est = m_iT2sel;

	my = GetGroup(Y_VAR);
	mx = GetGroup(X_VAR);	
    m_cX = columns(mx);
//-------------------SetTranslog-----------------------------
	if (m_iTl)
    {   print("Constructing Squares and Cross-Products...");
        if (m_iTl>1) m_cX=m_iTl;
		  
	mx=mx~0.5*(mx[][0:m_cX-1] .* mx[][0:m_cX-1]);//squares
    for (i = 0 ; i < m_cX-1 ; i++)
    {										
		mx=mx~(mx[][i+1:m_cX-1] .* mx[][i]);//cross-products
	}
        m_cX = columns(mx);
        print("done.\n");

	m_mTldata = my~mx;//constructed data
	}
	
	m_cT = m_iT2est - m_iT1est + 1;
    if (m_cT <= 2)
    {   eprint("Only ", m_cT, " observations. This is not enough.\n");
        return FALSE;
    }
//---------------------within-transformation-----------------
	decl mVars = my~mx;
	m_mVars = mVars;
	decl cSel = sizec(mVars);		
	decl mFirm, mStat = zeros(maxc(m_vID), cSel);
	//indiv. mean
	for (i = 0; i < maxc(m_vID); i++)// 
	{
	mFirm = selectifr(mVars, m_vID .== i+1);
	mStat[i][] = (meanc(mFirm[][]));//Mean of all selected vars;
	}
	m_mMeans = mStat;//indivual means [#firms][csel] 
	//Expansion of means
	decl vEmp, mData = m_ID ~ m_Ti ~ m_mMeans;//aux.matrix
	decl mAllMe = zeros(sizer(m_vID),cSel);//for the means
	decl mWtvars = zeros(sizer(m_vID), cSel);//for withins

	for (i=0; i < sizer(m_vID); i++)//LoopL[#obs]
	{
	//expansion
	vEmp = selectifr(mData, mData[][0] .== m_vID[i][]);//Selection
	mAllMe[i][] = vEmp[][2:];//exp. means
	//w-transformation
	mWtvars[i][] = mVars[i][] - mAllMe[i][];
	}
	m_mWithins = mWtvars;//w-trans. vars
	//final values
    m_var = mWtvars;			 
	m_mY = mWtvars[][0];//used internally
	m_cY = 1;

	mZ = GetGroup(Z_VAR);
    cZ = columns(mZ);//for model spec

    mU = GetGroup(U_VAR);
    cU = columns(mU);//for model spec
	
//cond.M2 & M3 
if (m_iMethod==2 || m_iMethod==3){
		if (cZ || cU){
			print("This model: neither Z-VAR nor U-VAR \n");
			return FALSE;
		}}//endM2
else
{//cond.M1

	if (cZ && cU)
	{
		print("WT model: either Z-VAR or U-VAR \n");
		return FALSE;
	}
	if (!cZ && !cU)
	{
		print("WT model: need variables (scaling) \n");
		return FALSE;

	}
	else if (cU)
	{
		m_iDist = TRUE;
		m_mZ = mU;
	}
	else
	{
		m_iDist = FALSE;
		m_mZ = mZ;
	}

	m_cZ = columns(m_mZ);
}//end cond.M1

	if (!m_cX)
    {   eprint("Need some regressors\n");
        return FALSE;
    }
	
    m_iModelStatus = MS_DATA;

    return TRUE;
break;

default:
//    decl  cp, vp, mh, i;
    m_iT1est = m_iT1sel;  m_iT2est = m_iT2sel;

    m_mU = GetGroup(U_VAR);
    m_cU = columns(m_mU);
    m_mY = GetGroup(Y_VAR);
    m_cY = 1;

    if (m_cU < 1 && (m_fPrintDetails))
        print("mu is restricted to zero\n");

    m_var = m_mY~GetGroup(X_VAR);
    m_mZ = GetGroup(Z_VAR);
    m_cZ = columns(m_mZ);
    m_cX = columns(m_var)-1;

    if (!m_cZ)
    {   m_mZ= GetVar("Constant")[:rows(m_var)-1];
        m_cZ=1;
    }

    if (!m_cX)
    {   eprint("Need some regressors\n");
        return FALSE;
    }

    if (m_iTl)
    {   print("Constructing Squares and Cross-Products...");
        if (m_iTl>1) m_cX=m_iTl+1;
        m_var=m_var~0.5*(m_var[][2:m_cX] .* m_var[][2:m_cX]);
        for (i = 2 ; i < m_cX ; i++)
        {   m_var=m_var~(m_var[][i+1:m_cX] .* m_var[][i]);
        }
        m_cX = columns(m_var)-1;
        print("done.\n");

		m_mTldata = m_var[][0]~m_var[][1:];//Y~X, includes 'Constant'
    }
    m_cT = m_iT2est - m_iT1est + 1;
    if (m_cT <= 2)
    {   eprint("Only ", m_cT, " observations. This is not enough.\n");
        return FALSE;
    }

    m_iModelStatus = MS_DATA;

    return TRUE;
break;
}}

Sfa::GridSearch(const beta, const s)
{switch (m_iMethod){
case 4:
case 3:
    decl i, fnow, fbest, fr, par_buffer;
    fnow = fbest = -10000;

    for (i = 0.05 ; i < 0.99 ; i += 0.05)
    {   fr = s/(1 - i * 0.6366198); //2/pi

			par_buffer =				 
			
//			(beta[0]+0.7978846*sqrt(fr*i))|			 
			beta[0:m_cX-1]
			|log((1-i)*fr)
			|log(fr*i);

		if (!fSfa(par_buffer,&fnow,0,0))
        {   print("function evaluation failed at gamma ",i," !\n");
            fnow = -1000;
        }
        else if (fnow > fbest)
        {
            fbest = fnow;
            m_par = par_buffer;
        }
    }
	return TRUE;
break;

case 1:
    fnow = fbest = -10000;

    for (i = 0.05 ; i < 0.99 ; i += 0.01)
    {   fr = s/(1 - i * 0.6366198); //2/pi
//        if (m_cX<2)									  //only for a constant (MC simulations)
//            par_buffer = (beta[0]+0.7978846*sqrt(fr*i))
//            |log(sqrt((1-i)*fr))|log(sqrt(fr*i))
//            // |zeros(m_cZ+m_cU-1,1);
//            |zeros(m_cZ-1,1);
//        else	//this is relevant:
            par_buffer =				 //Here: adjustment of vP
//			(beta[0]+0.7978846*sqrt(fr*i))|			 
			beta[0:m_cX-1]
			|log((1-i)*fr)
			|log(fr*i)
            // |zeros(m_cZ+m_cU-1,1);//original sfamb
            |zeros(m_cZ,1);

	if (m_iDist)		  
            par_buffer =
//			(beta[0]+0.7978846*sqrt(fr*i))|			 
			beta[0:m_cX-1]
			|log((1-i)*fr)
			|log(fr*i)
            |zeros(m_cZ+1,1);//adjustment for 'mu'
			
		if (!fSfa(par_buffer,&fnow,0,0))
        //if (!fWtm(par_buffer,&fnow,0,0))//merging
        {   print("function evaluation failed at gamma ",i," !\n");
            fnow = -1000;
        }
        else if (fnow > fbest)
        {
            fbest = fnow;
            m_par = par_buffer;
        }
    }
	return TRUE;
break;

default:
//    decl i, fnow, fbest, fr, par_buffer;
    fnow = fbest = -1000;

    for (i = 0.05 ; i < 0.99 ; i += 0.01)
    {   fr = s/(1 - i * 0.6366198); //2/pi
        if (m_cX<2)
            par_buffer = (beta[0]+0.7978846*sqrt(fr*i))
            |log(sqrt((1-i)*fr))|log(sqrt(fr*i))
            |zeros(m_cZ+m_cU-1,1);
        else
            par_buffer = (beta[0]+0.7978846*sqrt(fr*i))
            |beta[1:m_cX-1]|log(sqrt((1-i)*fr))|log(sqrt(fr*i))
            |zeros(m_cZ+m_cU-1,1);
        if (!fSfa(par_buffer,&fnow,0,0))
        {   print("function evaluation failed at gamma ",i," !\n");
            fnow = -1000;
        }
        else if (fnow > fbest)
        {
            fbest = fnow;
            m_par = par_buffer;
        }
    }
    return TRUE;
break;
}}

// fcon(const avF, const vP){
//       avF[0]=matrix(1-vP[1]-vP[2]-vP[3]-vP[4]);
//       return 1;
//     }

Sfa::DoEstimation(vPar)
{switch (m_iMethod){
case 4:
  decl  temp,cp,sd,gam,time;

  if (m_vStart==0)
    ols2c(m_var[][0], m_var[][1:m_cX], &temp);

  else 
    temp=m_vStart;

  sd=
    ((m_var[][0]-m_var[][1:m_cX]*temp[:m_cX-1])'(m_var[][0]-m_var[][1:m_cX]*temp[:m_cX-1]))/m_cT;//'
  m_dLoglikOLS = -m_cT/2*(log(sd)+2.837877066409); /* log-likelihood OLS*/

  if (sizeof(temp) == m_cX) {  
// check for length, if TRUE, only function pars are given, hence  perform gridsearch
	if (!GridSearch(temp,sd))
	  print("Error in GridSearch");
      }
  else m_par = temp;
  
  if (m_fPrintDetails) 
    print("starting values: ",
			m_par[0:(m_cX-sizer(m_ID)-1)] | m_par[m_cX:(m_cX+1)]);

  cp = sizer(m_par);
  SetParCount(cp);
  m_mCovP = unit(cp);
    
  time=timer();
  SetResult(5);

  SetResult(MaxBFGS(fSfa, &m_par, &m_dLogLik, &m_mCovP, TRUE));
	  m_dLogLik *= m_cT;

	  if (m_fPrintDetails)
	    print("Elapsed time: ",timespan(time),"\n");
	  vPar=m_par;

		//SetResult(MAX_CONV); //uncomment to get SE's even after non-convergence
	  return vPar;
break;

case 2:
    decl mRes, iSig;
	cp = columns(m_var[][1:m_cX]);			  
    SetParCount(cp);//must be called for the other functions

	olsc(m_var[][0], m_var[][1:m_cX], &vPar, &m_mCovar);

	mRes = m_var[][0] - m_var[][1:m_cX] * vPar;
	decl iTim = sumc(m_Ti)/maxc(m_ID);//avg. T-i
    iSig = mRes'mRes / (maxc(m_ID)*(iTim-1)-m_cX);//sigma-e^2 (corrected for panel)
    m_mCovar *= iSig;														   

	decl iN = maxc(m_ID);
	decl iF = iN*iTim/2;
	decl iSigNC = mRes'mRes / (maxc(m_ID)*iTim);//without correction

	m_dLogLik = double(-iF*log(2*M_PI)-iF*log(iSigNC)-(1/(2*iSigNC))*mRes'mRes);//LLF-OLS
	m_par = vPar;

    SetResult(MAX_CONV);
	
    return vPar;
break;

case 3:
case 1:

  if (m_vStart==0)//standard
    ols2c(m_var[][0], m_var[][1:m_cX], &temp);

  else 
    temp=m_vStart;

  if (m_fPrintDetails) 
    print("starting values(OLS on x~): ", temp);

  sd=
    ((m_var[][0]-m_var[][1:m_cX]*temp[:m_cX-1])'(m_var[][0]-m_var[][1:m_cX]*temp[:m_cX-1]))/m_cT;//'
  m_dLoglikOLS = -m_cT/2*(log(sd)+2.837877066409); /* log-likelihood OLS*/

  if (sizeof(temp) == m_cX) {  
// check for length, if TRUE, only function pars are given, hence  perform gridsearch
	if (!GridSearch(temp,sd))
	  print("Error in GridSearch");
      }
  else m_par = temp;
  
//  if (m_fPrintDetails) 
//    print("starting values(after GridSearch): ", m_par);

  cp = sizer(m_par);
  SetParCount(cp);
  m_mCovP = unit(cp);
    
  time=timer();
  SetResult(5);//Modelbase-function

//  SetResult(MaxBFGS(fWtm, &m_par, &m_dLogLik, &m_mCovP, TRUE));//TRUE=numerical gradients
  SetResult(MaxBFGS(fSfa, &m_par, &m_dLogLik, &m_mCovP, TRUE));//merging
//	  m_dLogLik *= m_cT;//original sfamb

	  if (m_fPrintDetails)
	    print("Elapsed time: ",timespan(time),"\n");
	  vPar=m_par;

		//SetResult(MAX_CONV); //uncomment to get SE's even after non-convergence
	  return vPar;
break;

default:
//  decl  temp,cp,sd,gam,time;
    
  /*    if (!InitPar())
	{   eprint("Failed to load data\n");
        m_par = zeros(m_cX+m_cZ+m_cU+1,1);
        return;
	}*/

  if (m_vStart==0)
    ols2c(m_var[][0], m_var[][1:m_cX], &temp);
  else 
    temp=m_vStart;

  sd=
    ((m_var[][0]-m_var[][1:m_cX]*temp[:m_cX-1])'(m_var[][0]-m_var[][1:m_cX]*temp[:m_cX-1]))/m_cT;//'
  m_dLoglikOLS = -m_cT/2*(log(sd)+2.837877066409); /* log-likelihood OLS*/
//  println(sd~sizeof(temp)~m_cX);
  //sd=sqrt(sd);
  
  if (sizeof(temp) == m_cX) {
// check for length, if TRUE, only function pars are given, hence  perform gridsearch
	if (!GridSearch(temp,sd))
	  print("Error in GridSearch");
      }
  else m_par = temp;
  
  if (m_fPrintDetails) 
    print("starting values: ", m_par);

  cp = sizer(m_par);
  SetParCount(cp);
  m_mCovP = unit(cp);
    
  time=timer();
  SetResult(5);
  if (!m_bUse_maxSQP){
    if (m_bUse_maxSQPF){
      print("Using SQPF");
      SetResult(MaxSQPF(fSfa, &m_par, &m_dLogLik, &m_mCovP, FALSE,
			m_fcon_ge0,m_fcon_eq0,m_vLower,m_vUpper));
    }
    
    else 
      SetResult(MaxBFGS(fSfa, &m_par, &m_dLogLik, &m_mCovP, FALSE));
  }
  else SetResult(MaxSQP(fSfa, &m_par, &m_dLogLik, &m_mCovP, FALSE,
			m_fcon_ge0,m_fcon_eq0,m_vLower,m_vUpper));
  m_dLogLik *= m_cT;
  
  if (m_fPrintDetails)
    print("Elapsed time: ",timespan(time),"\n");
  vPar=m_par;
    //    SetResult(MAX_CONV); //uncomment to get SE's even after non-convergence
  return vPar;
break;
}}

Sfa::Covar()				 
{switch (m_iMethod){
case 3:
case 1:
  m_mCovar = diag(ones(m_cPar,1));
    if (GetResult() == MAX_CONV || GetResult() == MAX_WEAK_CONV)
    {															  
        decl mXprod, dfunc;
        //fSfa(m_par, &dfunc, &mXprod, 0);						  
        //mXprod = m_mScore * m_mScore' /m_cT; //'						  
		//if (Num2Derivative(fSfa, m_par, &m_mCovP))				  
//		if (Num2Derivative(fWtm, m_par, &m_mCovP))//relevant here				  
		if (Num2Derivative(fSfa, m_par, &m_mCovP))//merging
        {															  
            m_mCovP = invertgen(-m_mCovP,30) ;
			println("");
                if (m_fPrintDetails)									  
                    println("standard errors from Hessian.");
        }																  
        else															  
        {   m_mCovP = m_mCovP ;										  
                if (m_fPrintDetails)									  
                    println("standard errors from last update:");		  
        }																  
		
//		m_mCovar = m_mCovP / m_cT;//original sfamb
		m_mCovar = m_mCovP;//adjustment

//		if (m_bUseRobustStdErr)
//            m_mCovarRobust = m_cT * m_mCovar * mXprod * m_mCovar ;
    }																	  
break;

default:
  m_mCovar = diag(ones(m_cPar,1));
    if (GetResult() == MAX_CONV || GetResult() == MAX_WEAK_CONV)
    {
        decl mXprod, dfunc;
        //fSfa(m_par, &dfunc, &mXprod, 0);
        //println(rows(mXprod)," ",columns(mXprod));
        mXprod = m_mScore * m_mScore' /m_cT; //'
        //print(rows(mXprod)," ",columns(mXprod));exit(0);
        if (Num2Derivative(fSfa, m_par, &m_mCovP))
        {
            m_mCovP = invertgen(-m_mCovP,30) ; //println(m_fPrintDetails);
			println("");
                if (m_fPrintDetails)
                    println("standard errors from Hessian.");
        }
        else
        {   m_mCovP = m_mCovP ;
                if (m_fPrintDetails)
                    println("standard errors from last update:");
        }

		m_mCovar = m_mCovP / m_cT;
		
	if (m_iMethod==0)//orig.sfamb:
		{
		if (m_bUseRobustStdErr)
            m_mCovarRobust = m_cT * m_mCovar * mXprod * m_mCovar ;
		}			
    }
break;
}}

Sfa::GetParNames()
{switch (m_iMethod){
case 4:
case 3:
    decl i, j, iSize, asy = {}, asx = {};
	decl asvu = {"ln{\\sigma_v^2}","ln{\\sigma_u^2}"};

	GetGroupLagNames(Y_VAR, 1, 100, &asy);
    if (m_cX) // regressors
        GetGroupNames(X_VAR, &asx);
    if (m_iTl){  // print squares and crossproducts, if necessary; adjusted for WT model
        iSize = (m_iTl > 1 ? m_iTl : sizeof(asx));
        for (i = 0 ; i < iSize; i++)    //print squares
                asx ~= sprint(".5*",asx[i],"^2");
        for (i = 0 ; i < iSize-1 ; i++)    //print cross terms
        {
				for (j = i; j < iSize-1; j++)
						asx ~= sprint(asx[i],"*",asx[j+1]);
        }
	}

	asx ~= asvu;
	
	for (i = sizeof(asx); i < m_cPar; ++i)
        asx ~= sprint("Par ", "-%2d", i + 1);
    return asy ~ asx;
break;

case 2:
    decl asu = {}, asz= {};
	
	GetGroupLagNames(Y_VAR, 1, 100, &asy);
    if (m_cX) // regressors
        GetGroupNames(X_VAR, &asx);
    if (m_iTl){  // print squares and crossproducts, if necessary; adjusted for WT model
        iSize = (m_iTl > 1 ? m_iTl : sizeof(asx));
        for (i = 0 ; i < iSize; i++)    //print squares
                asx ~= sprint(".5*",asx[i],"^2");
        for (i = 0 ; i < iSize-1 ; i++)    //print cross terms
        {
				for (j = i; j < iSize-1; j++)
						asx ~= sprint(asx[i],"*",asx[j+1]);
        }
	}

	for (i = sizeof(asx); i < m_cPar; ++i)
        asx ~= sprint("Par ", "-%2d", i + 1);
    return asy ~ asx;
break;

case 1:
//  decl i, j, iSize, asy = {}, asx = {}, asu = {}, asz= {};
	asvu = {"ln{\\sigma_v^2}","ln{\\sigma_u^2}"};
	decl smu = {};
	
    GetGroupLagNames(Y_VAR, 1, 100, &asy);
    if (m_cX) // regressors
        GetGroupNames(X_VAR, &asx);

	if (m_iTl){  // print squares and crossproducts, if necessary; adjusted for WT model
        iSize = (m_iTl > 1 ? m_iTl : sizeof(asx));
        for (i = 0 ; i < iSize; i++)    //print squares
                asx ~= sprint(".5*",asx[i],"^2");
        for (i = 0 ; i < iSize-1 ; i++)    //print cross terms
        {
				for (j = i; j < iSize-1; j++)
						asx ~= sprint(asx[i],"*",asx[j+1]);
        }
	}
//	GetGroupNames(Z_VAR, &asz);//original WT

	if (m_iDist)
	{GetGroupNames(U_VAR, &asz);}
	else
	{GetGroupNames(Z_VAR, &asz);}

	if (m_iDist)
		smu = "mu";

	asx ~= asvu ~= asz ~= smu;

	for (i = sizeof(asx); i < m_cPar; ++i)
        asx ~= sprint("Par ", "-%2d", i + 1);
    return asy ~ asx;
break;

default:
//    decl i, j, iSize, asy = {}, asx = {}, asu = {}, asz= {};
    GetGroupLagNames(Y_VAR, 1, 100, &asy);
    if (m_cX) // regressors
        GetGroupNames(X_VAR, &asx);
    if (m_iTl){  // print squares and crossproducts, if necessary
        iSize = (m_iTl > 1 ? m_iTl : sizeof(asx)-1);
        for (i = 0 ; i < iSize; i++)    //print squares
                asx ~= sprint(".5*",asx[i+1],"^2");
        for (i = 1 ; i < iSize ; i++)    //print cross terms
        {
                for (j = i; j < iSize; j++)
                        asx ~= sprint(asx[i],"*",asx[j+1]);
        }

    }
    asx ~= "ln{\\sigma_v}";
    if (m_cZ > 1)
        GetGroupNames(Z_VAR, &asz);
    else asz = "ln{\\sigma_u}";
    asx ~= asz;
    if (m_cU){
        GetGroupNames(U_VAR, &asu);
        asx ~= asu;
    }
    for (i = sizeof(asx); i < m_cPar; ++i)
        asx ~= sprint("Par ", "-%2d", i + 1);
    return asy ~ asx;
break;
}}

Sfa::Output()
{switch (m_iMethod){
case 4:
	println("-TFE model-");
	/*Modelbase::Output(); =>include Modelbase source code: */
	if (!OutputHeader(GetPackageName()))     // returns FALSE if no estimation
        return;
    //OutputPar();{...									
    decl i, tval, abstval, tprob;
    decl aspar = GetParNames();
    decl aspartypes = GetParTypes();
	decl cpartypes = isarray(aspartypes) ? sizeof(aspartypes) : 0;
    decl vstderr = GetStdErr(), vrobstderr = GetStdErrRobust(), bcovarrobust = FALSE;
    decl ct = GetcT();
    decl mpar, vp = GetPar();
	decl cdfloss = GetcDfLoss();
    //mpar = vp ~ vstderr ~ vp ./ vstderr; //original modelbase
    decl mparaux = vp ~ vstderr ~ vp ./ vstderr;//K x 3
	
	mpar = mparaux[0:(m_cX-sizer(m_ID)-1)][] | mparaux[m_cX:m_cX+1][];//(betas+sv2+su2) x 3	
	decl cp = rows(mpar);

	if (rows(vrobstderr) > 1)
    {   mpar ~= vrobstderr ~ vp ./ vrobstderr;
        bcovarrobust = TRUE;
    }
    print("%29s", "Coefficient", "%11s", "Std.Error");
    if (bcovarrobust)
        print("%11s", "robust-SE");
    println("%9s", "t-value", "%8s", "t-prob");
    
    for (i = 0; i < cp; ++i)
    {
        if (!m_vIsFreePar[i] && mpar[i][0] == 0)
            continue;
        tval = mpar[i][2];
		if (i < cpartypes)
	        print("%-13s", aspar[i], " ", "%-2s", aspartypes[i], "%#13.6g", mpar[i][0]);
		else
	        print("%-16s", aspar[i], "%#13.6g", mpar[i][0]);
        if (!m_vIsFreePar[i])
            println("    (fixed)");
        else if (mpar[i][1] > 0)
        {
            print("%#11.4g", mpar[i][1]);
            if (bcovarrobust)
            {
                tval = mpar[i][4];
                print("%#11.4g", mpar[i][3]);
            }
			abstval = fabs(tval);
            tprob = 2 * tailt(abstval, ct - cdfloss);
            println(abstval <= 1e-4 ? "%9.2f" : abstval >= 1e3 ? "%#9.4g"
				: "%#9.3g", tval, "%8.3f", tprob);
        }
        else
            print("\n");
    }//end of 'OutputPar()'

    OutputLogLik();//orig. ModelBase
	
	decl sv = sqrt(exp(m_par[m_cX]));
	decl su = sqrt(exp(m_par[m_cX+1]));
	decl lambda = double(su/sv);
	println("lambda  ","%21.4g", lambda);
break;

case 3:
	println("-CFE model-");

	Modelbase::Output();

	sv = sqrt(exp(m_par[m_cX]));
	su = sqrt(exp(m_par[m_cX+1]));
	lambda = double(su/sv);
	println("lambda  ","%21.4g", lambda);	
break;

case 2:
	println("\n-LSDV model-");

	if (!OutputHeader(GetPackageName()))//orig. Modelbase
        return;

	println("");
    Modelbase::OutputPar();
/*adjustment of 'OutputLogLik()'*/
    ct = GetcT();
    cdfloss = GetcDfLoss();
/*LSDV part*/
	decl mRes = m_var[][0] - m_var[][1:m_cX] * m_par;
	decl iSumEps = double(mRes'mRes);
	decl iTim = sumc(m_Ti)/maxc(m_ID);//avg. T-i
	decl iSigma = double(sqrt(mRes'mRes / (maxc(m_ID)*(iTim-1)-m_cX)));//sigma-e^2 (corrected)
	
    println("\nlog-likelihood", "%15.9g", m_dLogLik);
    println("no. of observations", "%10d", ct,
          "  no. of parameters  ", "%10d", cdfloss);
//    println("AIC.T         ", "%15.9g", -2 * m_dLogLik + 2 * cdfloss,
//          "  AIC           ", "%15.9g", (-2 * m_dLogLik + 2 * cdfloss) / ct);
    println("AIC1 (all obs)", "%15.9g", double(-2 * m_dLogLik + 2 * (cdfloss+1)),//acc to stata
          "  AIC2          ", "%15.9g", double(log(iSumEps/ct) + (2*(m_cX+maxc(m_ID)))/ct));//acc to limdep

	if (m_cY == 1)
    {
        decl asy;
        GetGroupLagNames(Y_VAR, 0, 0, &asy);
        println("%-16s", sprint("mean(", asy[0],")"), "%13.6g", double(meanc(m_mY)),
                "%-18s", sprint("  var(", asy[0],")"),  "%13.6g", double(varc(m_mY)));
    }

	println("sigma_e       ","%15.6g", iSigma,
		  "  SSR           ","%15.6g", iSumEps);
break;

case 1:
	println("-WT model-");

    Modelbase::Output();

	sv = sqrt(exp(m_par[m_cX]));
	decl su2 = exp(m_par[m_cX+1]);
	
	if (m_iDist==FALSE)
		su2 = meanc(exp(m_par[m_cX+1] + m_mZ*m_par[m_cX+2:m_cX+m_cZ+1]));
	su = sqrt(su2);
	lambda = double(su/sv);
	println("lambda  ","%21.4g", lambda);
break;

default:
	println("-Pooled model-");

    Modelbase::Output();
    decl v_v=exp(2*m_par[m_cX]);
    decl s_u2= double(meanc(exp(2*m_mZ*m_par[m_cX+1:m_cX+m_cZ])));
    println("\\gamma:       ", "%15.4g", s_u2/(v_v+s_u2),
    "  VAR(u)/VAR(total)", "%12.4g", ((M_PI-2)/M_PI*s_u2)/(v_v+(M_PI-2)/M_PI*s_u2));
    println("Test of one-sided err","%8.5g",
        double(2*(m_dLogLik-m_dLoglikOLS)),"  mixed Chi^2 !!");
break;
}}

Sfa::TestGraphicAnalysis()
{switch (m_iMethod){
case 4:
case 3:
case 2:
case 1:
    decl mEff = TE();
    // Histogram and boxplot
    DrawDensity(0, (mEff[][0])', "Technical efficiency", 0, 1, 0);
	DrawBoxPlot(1, (mEff[][0])', "Technical efficiency");
	ShowDrawWindow();
break;

default:
//    decl mEffs = TE(m_dAlpha), i, mFan, vQuan, mTmp;
    decl mEffs = TEint(m_dAlpha), i, mFan, vQuan, mTmp;
	// setup fan matrix 
	decl dStep = 0.05;
	mFan=zeros(m_cT,2/dStep-1);
	vQuan = range(dStep, 1-dStep , dStep);
	for (i=0; i < 1/dStep-1; ++i){
//		mTmp = TE(vQuan[i]);
		mTmp = TEint(vQuan[i]);
		mFan[][i]= mTmp[][1] - mTmp[][0];
		mFan[][2/dStep-2-i]= mTmp[][2] - mTmp[][0];
		}
//	print(range(m_dAlpha/2, 1-m_dAlpha/2 , .025)|mFan[0][]);exit(0);
        // Histogram
    DrawDensity(0, (mEffs[][0])', "Technical efficiency", 0, 1, 0);//'
//    ShowDrawWindow();
    DrawBoxPlot(1, (mEffs[][0])', "Technical efficiency");//'
    DrawMatrix(2, sortbyc(mEffs,0)', {"TE", "Lower bound", "Upper bound"}, 0, 1);//'
    DrawLegend(2, 1, 1, FALSE);
	mFan = sortbyc((mEffs[][0]~mFan),0);
    ShowDrawWindow();
    Draw(0, mFan[][0]');
    DrawZ(mFan[][1:]', "FAN2", ZMODE_FAN, range(dStep/2, 1-dStep/2 , dStep/2));
 	ShowDrawWindow();
//    DrawAdjust(ADJ_MINMAX, 0, 1);
break;
}}

Sfa::GetResults(const par, const eff, const fct, const v)
{
    decl temp;
//    eff[0]=TE(m_dAlpha);        // return estimated efficiencies
    eff[0]=TEint(m_dAlpha);        // return estimated efficiencies
/*    par[0]=m_par~sqrt(diagonal(m_mCovP)')      // return estimated coefficients
        ~2*tailt(fabs(m_par)./ sqrt(diagonal(m_mCovP)'),m_cT-sizer(m_par));
*/
    par[0]=m_par;
    if (m_bUseRobustStdErr) temp = GetStdErrRobust(); else temp = GetStdErr();
    par[0] ~= temp;
    par[0]~= temp .== 0 .? 0 .: 2*tailt(fabs(m_par)./temp,m_cT-sizer(m_par));
    if (m_fPrintDetails)
        print(par[0]);
    v[0]=m_mCovP;       //return variance-covariance matrix
        // return some statistics (3x2):
        // first col likelihood: llfOLS | llfSFA | LR-Test of one-sided error
        // second col variance: VAR(u) |Var(v) | correct decomposition
    decl v_v=exp(2*m_par[m_cX]);
    decl v_u=(M_PI-2)/M_PI*meanc(exp(2*m_mZ*m_par[m_cX+1:m_cX+m_cZ]));
    fct[0]=(m_dLoglikOLS|m_dLogLik|2*(m_dLogLik-m_dLoglikOLS))
        ~(v_u|v_v|(v_u/(v_v+v_u)));
}

Sfa::SetPrintDetails(const bool)
{
  m_fPrintDetails=bool;
  return 1;
}

Sfa::SetUse_maxSQPF(const bool)
{
  if (bool)
    m_bUse_maxSQP=FALSE;
  m_bUse_maxSQPF=bool;
  return 1;
}

Sfa::SetRobustStdErr(const bool)
{
  m_bUseRobustStdErr = bool;
  return 1;
}

Sfa::AiHat()
{switch (m_iMethod){
case 4:
	return m_par[(m_cX-sizer(m_ID)):m_cX-1];
break;

case 3:
	return m_mMeans[][0]
	- m_mMeans[][1:]*m_par[0:(m_cX-1)]
	+ sqrt( 2/M_PI) * sqrt( exp(m_par[m_cX+1]));//eq(21) mean-adj. alpha
break;

case 2:
	return m_mMeans[][0]-m_mMeans[][1:]*m_par[0:(m_cX-1)];
break;

case 1:
	decl vjlms = Ineff();
/*inactive code related to eq(31) of W&H(2010), not verified!*/

//	decl mu=0;	  
//	if (m_iDist)
//	 mu = m_par[m_cX+m_cZ+2];
//	decl vhit = exp(m_mZ*m_par[m_cX+2:m_cX+m_cZ+1]);//h_it
//	decl veit = m_mVars[][0]-m_mVars[][1:m_cX]*m_par[0:(m_cX-1)];//eps_it from "orignal" data; (m_mVars from InitData)
//	decl mData = m_vjlms~vhit~(vhit.^2)~(veit.*vhit);
	decl vaihat, i, mFirm, vmine = zeros(maxc(m_vID),1);//for mean of jlms
//	decl vmhit = zeros(maxc(m_vID),1);//meanc(h_it)
//	decl vhsq = zeros(maxc(m_vID),1);//sumc(h_it^2)
//	decl vpeh = zeros(maxc(m_vID),1);//sumc(eps*hit)
	for (i = 0; i < maxc(m_vID); i++)//loopL[#obs]
	{
//	mFirm = selectifr(mData, m_vID .== i+1);
	mFirm = selectifr(vjlms, m_vID .== i+1);
	vmine[i][] = (meanc(mFirm[][0]));//mean of u_it
//	vmhit[i][] = (meanc(mFirm[][1]));
//	vhsq[i][] = (sumc(mFirm[][2]));
//	vpeh[i][] = (sumc(mFirm[][3]));
	}
//	decl vaihat = m_mMeans[][0]-m_mMeans[][1:]*m_par[0:(m_cX-1)]+vmine;//robust formula, paper eq(2.5)

	vaihat = m_mMeans[][0]-m_mMeans[][1:]*m_par[0:(m_cX-1)]+vmine;//robust formula, paper eq(2.5)

	return vaihat;
//	//for results: //related to eq(31)
//	decl im1 = zeros(sizer(m_ID),1);//aux1
//	decl im2 = zeros(sizer(m_ID),1);//aux2
//	decl mustar = zeros(sizer(m_ID),1);//mu_***
//	decl sigstar = zeros(sizer(m_ID),1);//sig_***
//	decl ahat = zeros(sizer(m_ID),1);//a_i_hat
//	//actual calc:
//	for (i = 0; i < sizer(m_ID); i++)//loopL[#firms]  <=> eq.s (31-33)
//	{
//	im1[i][] = (sqrt(exp(m_par[m_cX]))^-2*m_Ti[i][])*vpeh[i][];//part of eq(32)numerator
//	im2[i][] = (sqrt(exp(m_par[m_cX]))^-2*m_Ti[i][])*vhsq[i][];//part of eq(32)denom.
//	//mustar, sigstar
//	mustar[i][]= (mu*sqrt(exp(m_par[m_cX+1]))^-2 - im1[i][])/(im2[i][] + sqrt(exp(m_par[m_cX+1]))^-2);//eq(32)
//	sigstar[i][]= sqrt(exp(m_par[m_cX]))^2*m_Ti[i][] /
//				(vhsq[i][] + sqrt(exp(m_par[m_cX]))^2*m_Ti[i][]*sqrt(exp(m_par[m_cX+1]))^-2);//eq(33) 
//	}
//	//aux. matrices
//	decl im3 = densn(mustar ./ sqrt(sigstar));
//	decl im4 = probn(mustar ./ sqrt(sigstar));
//
//	decl vaihat2 = m_mMeans[][0]-m_mMeans[][1:]*m_par[0:(m_cX-1)]
//				+ mustar.*vmhit + sqrt(sigstar).*vmhit.*im3./im4;//eq(31), alpha-i-hat
//	return vaihat~vaihat2;
break;

default:
	return FALSE;
break;
}}

Sfa::Elast(const sXname)
{switch (m_iMethod){
case 2:
  decl i, maxParIdx, sTmp, vIdx, vData; //asX = {"Constant",sXname}, 
  decl asX = {sXname};
  decl asNames;
  asNames=GetParNames();
  maxParIdx = sizec(asNames);
  vIdx = strfind(asNames,sXname);// we need the index of the first order term
  //find indices of second order terms
  for (i=vIdx+1;i< maxParIdx; ++i){
  	vIdx ~= strfind(asNames[i],sXname) == -1 ? <> : i;
  }
  for (i=2;i<sizerc(vIdx);++i){
    sTmp=asNames[vIdx[i]];//println(sTmp);
  	asX ~= strfind(sTmp,sXname) == 0
		? sTmp[sizeof(sXname)+1:]
		: sTmp[:strfind(sTmp,"*")-1];
  }
  vData = ones(GetcT(),1)~GetVar(asX);//vector of ones(panel models: "Constant" was set off)
  decl vEstimates = GetPar()[vIdx];
//  decl sigma2, vCovar = GetCovarRobust()[vIdx][vIdx];
  decl sigma2, vCovar = GetCovar()[vIdx][vIdx];
  decl vTvalue=zeros(GetcT(),1);
	for (i=0;i<GetcT();i++){
	  sigma2 = (vData[i][] * vCovar * (vData[i][])');
	  vTvalue[i]=((vData * vEstimates)[i]/sqrt(sigma2));
	}
  
  return (vData * vEstimates)~(vTvalue);
break;

case 4:
case 3:
case 1:
  asX = {sXname};
  asNames=GetParNames();
// get index of sigma_v^2 to exclude re-occuring X in Z or U
// when integrating this into the package, use a better var
  maxParIdx = strfind(asNames,"ln{\sigma_v^2}");
  if (maxParIdx == -1){
    eprint("Error with finding sigma_v in the list of vars, exiting...");
	exit(1);
  }
// we need the index of the first order term
  vIdx = strfind(asNames,sXname);
// find indices of second order terms
  for (i=vIdx+1;i< maxParIdx; ++i){
  	vIdx ~= strfind(asNames[i],sXname) == -1 ? <> : i;
  }
  for (i=2;i<sizerc(vIdx);++i){
    sTmp=asNames[vIdx[i]];//println(sTmp);
  	asX ~= strfind(sTmp,sXname) == 0
		? sTmp[sizeof(sXname)+1:]
		: sTmp[:strfind(sTmp,"*")-1];
  }
  vData = ones(GetcT(),1)~GetVar(asX);//vector of ones added
  vEstimates = GetPar()[vIdx];
//  decl sigma2, vCovar = GetCovarRobust()[vIdx][vIdx];//not for WT
  vCovar = GetCovar()[vIdx][vIdx];
  vTvalue=zeros(GetcT(),1);
	for (i=0;i<GetcT();i++){
	  sigma2 = (vData[i][] * vCovar * (vData[i][])');
	  vTvalue[i]=((vData * vEstimates)[i]/sqrt(sigma2));
	}
  
  return (vData * vEstimates)~(vTvalue);
break;

default:
//  decl i, maxParIdx, sTmp, vIdx, asX = {"Constant",sXname}, vData; 
//  decl asNames;
  asX = {"Constant",sXname};
  asNames=GetParNames();
// get index of sigma_v to exclude re-occuring X in Z or U
// when integrating this into the package, use a better var
  maxParIdx = strfind(asNames,"ln{\sigma_v}");
  if (maxParIdx == -1){						  
    eprint("Error with finding sigma_v in the list of vars, exiting...");
	exit(1);
  }
// we need the index of the first order term (i.e. of the considered variable)
  vIdx = strfind(asNames,sXname);
// find indices of second order terms
  for (i=vIdx+1;i< maxParIdx; ++i){
  	vIdx ~= strfind(asNames[i],sXname) == -1 ? <> : i;
  }
  for (i=2;i<sizerc(vIdx);++i){
    sTmp=asNames[vIdx[i]];
  	asX ~= strfind(sTmp,sXname) == 0
		? sTmp[sizeof(sXname)+1:]
		: sTmp[:strfind(sTmp,"*")-1];
  }
  if (sizerc(asX) == 2) asX = {"Constant"}; // accommodate CD ...
  //variables and coeffs  
  vData = GetVar(asX);//constant(coeff(.*1)
//  decl vEstimates = GetPar()[vIdx];
//  decl sigma2, vCovar = GetCovarRobust()[vIdx][vIdx];
//  decl vTvalue=zeros(GetcT(),1);
  vEstimates = GetPar()[vIdx];
  vCovar = GetCovarRobust()[vIdx][vIdx];
  vTvalue=zeros(GetcT(),1);
	for (i=0;i<GetcT();i++){ 
	  sigma2 = (vData[i][] * vCovar * (vData[i][])');
	  vTvalue[i]=((vData * vEstimates)[i]/sqrt(sigma2));
	}

  return (vData * vEstimates)~(vTvalue);
break;
}}

Sfa::GetTldata()
{
  return m_mTldata;
}

Sfa::GetLLFi()
{switch (m_iMethod){
case 4:
	decl i, vLLFi = zeros(maxc(m_vID),1);

	for (i = 0; i < maxc(m_vID); i++){ 
		vLLFi[i][] = m_dLogLik/m_cT * m_Ti[i][];}
	return vLLFi;
break;

case 2:
	return ones(GetcT(),1) * GetLogLik()/GetcT();
break;

case 3:
case 1:
	return m_lnLW;
break;

default:
	decl my=0;
	decl u  = m_var[][0]-m_var[][1:m_cX]*m_par[0:(m_cX-1)];

    if (m_cU>0)
        my = -m_mU[][]*m_par[m_cX+m_cZ+1:m_cX+m_cZ+m_cU];

	decl uu = (u-my).*(u-my);
	decl s  = sqrt(exp(2*m_par[m_cX])
        		+exp(2*m_mZ*m_par[m_cX+1:m_cX+m_cZ]));
				
	decl l = exp(m_mZ*m_par[m_cX+1:m_cX+m_cZ]-m_par[m_cX]);
	decl ra = my./(s.*l);
    
	return (-.9189385-0.5.*log(s.^2)-uu./(2*s.^2) + log(probn(-l.*u./s-ra))
    			-log(probn(-ra.*sqrt(1+l.^2))));
break;
}}

Sfa::GetMeans()
{
  return m_mMeans;
}

Sfa::GetWithins()
{
  return m_mWithins;
}