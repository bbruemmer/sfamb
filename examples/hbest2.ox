#include <oxstd.h>
#include <packages/gnudraw/gnudraw.h>
//#import <packages/sfamb/sfamb>	
#include "../sfamb/sfamb.ox"

main(){								   
/*new object of class 'Sfa', Load data*/
    decl fob = new Sfa();fob.Load("USDAafrica.xls");
//	fob.Info(); exit(1);	
/*Model specification*/  
	fob.SetMethod(CFE); fob.SetConstant();
/*Identification of panel structure*/
	fob.Ident(fob.GetVar("ID"), fob.GetVar("time"));
/*Data preparation*/	
	decl inorm = 1;//choose logs; or normalization and logs
	fob.Renew(fob.PrepData(fob.GetVar("output"), inorm), "lny");
	fob.Renew(fob.PrepData(fob.GetVar("labour"), inorm), "lnlab");
	fob.Renew(fob.PrepData(fob.GetVar("land"), inorm), "lnland");
	fob.Renew(fob.PrepData(fob.GetVar("machinery"), inorm), "lnmac");
	fob.Renew(fob.PrepData(fob.GetVar("fertilizer"), inorm), "lnfert");
	fob.Renew(fob.GetVar("time") - meanc(fob.GetVar("time")), "trend");//trend norm.
//	fob.Info(); exit(1);
	
/*Set up model*/
    fob.Select("Y", {"lny", 0, 0});			// Select dependent variable

    fob.Select("X", {						// Select regressors
						"Constant", 0, 0, 
						"lnlab", 0, 0, 
						"lnland", 0, 0, 
						"lnmac", 0, 0, 
						"lnfert", 0, 0,
//						"time", 0, 0						
						"trend",0 ,0						
						});

	// Select estimation sample
	fob.SetSelSample(-1, 1, -1, 1);  // full sample
    fob.SetPrintSfa(TRUE);
    MaxControl(1000, 10, TRUE);
	fob.SetTranslog(1);
	fob.Estimate();

	fob.Renew(fob.TE(),{"TE"});
	fob.Renew(fob.Ineff(),{"jlms"});
	fob.TestGraphicAnalysis();	

	delete fob;
}