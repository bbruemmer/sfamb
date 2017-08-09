#include <oxstd.h>
#import  <packages/sfamb/sfamb>

/*
 Example of function DropGroupIf(const mifr) using an arbitrary condition.
 Use after function Ident(const vID, const vPer).
*/


main(){
/*new object of class 'Sfa', Load data*/
    decl fob = new Sfa();fob.Load("Sample2.xls");
//	fob.Info(); exit(1);

/*Model specification*/  
	fob.SetMethod(LSDV);
	fob.SetConstant();
/*Identification of panel structure*/
	fob.Ident(fob.GetVar("id"), fob.GetVar("time"));

/*Exemplary use of DropGroupIf(const mifr)*/
	decl iMeanY = meanc(fob.GetVar("y"));
	decl iStdvY = sqrt(varc(fob.GetVar("y")));

	fob.DropGroupIf(fob.GetVar("y") .> (iMeanY + iStdvY));

/*Data preparation*/	
	fob.Renew(fob.PrepData(fob.GetVar("y"), 0), "lny");
	fob.Renew(fob.PrepData(fob.GetVar("x1"), 0), "lnx1");
	fob.Renew(fob.PrepData(fob.GetVar("x2"), 0), "lnx2");
	fob.Renew(fob.PrepData(fob.GetVar("x3"), 0), "lnx3");
	
/*Set up model*/
    fob.Select(Y_VAR, {"lny",0,0});			// Select dependent variable

    fob.Select(X_VAR, {						// Select regressors
						"Constant",0,0, 
						"lnx1",0,0, 
						"lnx2",0,0,
						"lnx3",0,0,
						"time",0,0
						});

	// Select estimation sample
	fob.SetSelSample(-1, 1, -1, 1);  // full sample
	fob.SetPrintSfa(TRUE);
	MaxControl(1000,10,TRUE);
	fob.SetTranslog(0);

	fob.Estimate();

	delete fob;
}