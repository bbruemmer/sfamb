#include <oxstd.h>
#include <packages/gnudraw/gnudraw.h>
#import <packages/sfamb/sfamb>	

main(){								   
/*new object of class 'Sfa', Load data*/
    decl fob = new Sfa(); fob.Load("USDAafrica.xls");
//	fob.Info(); exit(1);	
/*Model specification*/  
	fob.SetMethod(POOLED); fob.SetConstant();
/*Data preparation*/	
	decl inorm = 1;//choose logs; or normalization and logs
	fob.Renew(fob.PrepData(fob.GetVar("output"), inorm), "lny");
	fob.Renew(fob.PrepData(fob.GetVar("labour"), inorm), "lnlab");
	fob.Renew(fob.PrepData(fob.GetVar("land"), inorm), "lnland");
	fob.Renew(fob.PrepData(fob.GetVar("machinery"), inorm), "lnmac");
	fob.Renew(fob.PrepData(fob.GetVar("fertilizer"), inorm), "lnfert");
	fob.Renew(fob.GetVar("time") - meanc(fob.GetVar("time")), "trend");
//	fob.Info(); exit(1);	
	
/*Set up model*/
    fob.Select(Y_VAR, {"lny", 0, 0});			// Select dependent variable

    fob.Select(X_VAR, {						// Select regressors
						"Constant", 0, 0, 
						"lnlab", 0, 0, 
						"lnland", 0, 0, 
						"lnmac", 0, 0, 
						"lnfert", 0, 0,
//						"time", 0, 0						
						"trend",0 ,0						
						});

	fob.Select(U_VAR, {						//Shifting mu	
						"Constant", 0, 0         
						});
	
	fob.Select(Z_VAR, {					 	//Scaling sigma_u		
						"Constant", 0, 0,
						"lnlab", 0, 0, 
						"lnland", 0, 0, 
						"lnmac", 0, 0, 
						"lnfert", 0, 0
						});
						
	fob.SetSelSample(-1, 1, -1, 1);  // full sample
    fob.SetPrintSfa(TRUE);
    MaxControl(1000, 10, TRUE);
	fob.SetTranslog(1);
	fob.Estimate();

	fob.Renew(fob.TEint(0.05), {"TE", "lower", "upper"});
	fob.Renew(fob.Ineff(), {"jlms"});
	fob.Save("out.xls");

	fob.SetConfidenceLevel(0.05);	
	fob.TestGraphicAnalysis();

    delete fob;
}