#include <oxstd.h>
#include <packages/gnudraw/gnudraw.h>
//#include  <packages/sfamb/sfamb.ox>	
#import <packages/sfamb/sfamb>

main(){								   
/*new object of class 'Sfa', Load data*/
    decl fob = new Sfa();fob.Load("Sample2.xls");
//	fob.Info();	exit(1);

/*Model specification*/  
	fob.SetMethod(POOLED);
//	fob.SetMethod(LSDV);
	fob.SetConstant();
/*Identification of panel structure*/
	fob.Ident(fob.GetVar("id"), fob.GetVar("time"));
/*Data preparation*/	
	decl inorm = 0;//choose logs; or normalization and logs
	fob.Renew(fob.PrepData(fob.GetVar("y"), inorm), "lny");
	fob.Renew(fob.PrepData(fob.GetVar("x1"), inorm), "lnx1");
	fob.Renew(fob.PrepData(fob.GetVar("x2"), inorm), "lnx2");
	fob.Renew(fob.PrepData(fob.GetVar("x3"), inorm), "lnx3");
	fob.Renew(fob.GetVar("time")-meanc(fob.GetVar("time")), "trend");//trend norm.
//	fob.Info();	exit(1);

/*Set up model*/
    fob.Select(Y_VAR, {"lny",0,0});			// Select dependent variable

    fob.Select(X_VAR, {						// Select regressors
						"Constant",0,0, 
						"lnx1",0,0, 
						"lnx2",0,0,
						"lnx3",0,0,
						"time",0,0
//						"trend",0,0
						});	     
	
//	fob.Select(Z_VAR, {							
//						"Constant",0,0,
//						"z1",0,0 
//						});

//	fob.Select(U_VAR, {						
//						"Constant",0,0,          
//						"lnx1",0,0, 
//						"lnx2",0,0,
//						"lnx3",0,0
//						});
						
	// Select estimation sample
	fob.SetSelSample(-1, 1, -1, 1);  // full sample
    fob.SetPrintSfa(TRUE);
    MaxControl(1000,10,TRUE);
//	fob.SetPrintDetails(1);//display starting values and elapsed time
	fob.SetTranslog(inorm);

	//Estimate the model, get TE scores
	fob.Estimate();
	fob.Renew(fob.TE(),{"TE"});
	
/*Additional options*/
//	fob.Save("out1.xls");

//	fob.TestGraphicAnalysis();
//	savemat("jlms.xls",fob.Ineff());
//	savemat("Tldata.xls",fob.GetTldata());
//	savemat("IndivLLF.xls", fob.GetLLFi());
//	savemat("Elast.xls",fob.Elast("lnx2"));
//	savemat("covar.xls", fob.GetCovar());
//	fob.Renew(fob.Ineff(),{"jlms"});
//	fob.Renew(fob.GetResiduals(),{"eps"});

/*Pooled model*/
//	fob.Renew(fob.TEint(0.05), {"TE", "lower", "upper"});
//	fob.SetConfidenceLevel(0.05);
//	fob.TestGraphicAnalysis();
//	savemat("covarRob.xls", fob.GetCovarRobust());

/*Panel models*/
//	fob.Renew(fob.IDandPer(), {"Firmno","groupti"}); //Append firm number and individual T_i
//	savemat("alpha.xls", fob.AiHat());//individual alpha_i

//	savemat("Means.xls",fob.GetMeans());
//	savemat("Withins.xls",fob.GetWithins());

//	fob.Save("out2.xls");
	
    delete fob;
}
