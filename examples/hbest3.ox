#include <oxstd.h>
#include <packages/gnudraw/gnudraw.h>
//#import <packages/sfamb/sfamb>	
#include  "../sfamb/sfamb.ox"	

main(){								   
/*new object of class 'Sfa', Load data*/
    decl fob = new Sfa(); fob.Load("Sample2.xls");
/*Model specification*/  
	fob.SetMethod(WT); fob.SetConstant();
/*Identification of panel structure*/
	fob.Ident(fob.GetVar("id"), fob.GetVar("time"));
/*Data preparation*/	
	fob.Renew(fob.PrepData(fob.GetVar("y"), 0), "lny");
	fob.Renew(fob.PrepData(fob.GetVar("x1"), 0), "lnx1");
	fob.Renew(fob.PrepData(fob.GetVar("x2"), 0), "lnx2");
	fob.Renew(fob.PrepData(fob.GetVar("x3"), 0), "lnx3");
//	fob.Renew(fob.IDandPer(), {"Firmno", "groupti"});	
//	fob.Info();	exit(1);

/*Set up model*/
    fob.Select("Y", {"lny", 0, 0});			// Select dependent variable

    fob.Select("X", {						// Select regressors
						"Constant", 0, 0, 
						"lnx1", 0, 0, 
						"lnx2", 0, 0,
						"lnx3", 0, 0,
						"time", 0, 0
						});	     
	
	fob.Select("Z", {							
						"Constant", 0, 0,
						"z1", 0, 0 
						});

	// Select estimation sample
	fob.SetSelSample(-1, 1, -1, 1);  // full sample
    fob.SetPrintSfa(TRUE);
    MaxControl(1000, 10, TRUE);
	fob.SetTranslog(0);

	//Estimate the model, get TE scores
	fob.Estimate();
	fob.Renew(fob.TE(), {"TE"});
	fob.Renew(fob.Ineff(), {"jlms"});

//	fob.Save("out2.xls");
//	savemat("alpha2.xls", fob.AiHat());
	fob.TestGraphicAnalysis();

    delete fob;
}