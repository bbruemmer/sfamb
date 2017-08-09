#include <oxstd.h>
#include <packages/gnudraw/gnudraw.h>
#import <packages/sfamb/sfamb>

main(){								   
/*new object of class 'Sfa', Load data*/
    decl fob = new Sfa(); fob.Load("Sample2.xls");
//	fob.Info();	exit(1);

/*Model specification*/  
	fob.SetMethod(POOLED); fob.SetConstant();
/*Data preparation*/	
	decl inorm = 1;//choose logs; or normalization and logs
	fob.Renew(fob.PrepData(fob.GetVar("y"), inorm), "lny");
	fob.Renew(fob.PrepData(fob.GetVar("x1"), inorm), "lnx1");
	fob.Renew(fob.PrepData(fob.GetVar("x2"), inorm), "lnx2");
	fob.Renew(fob.PrepData(fob.GetVar("x3"), inorm), "lnx3");
	fob.Renew(fob.GetVar("time") - meanc(fob.GetVar("time")), "trend");//trend norm.
//	fob.Info();	exit(1);

/*Set up model*/
    fob.Select(Y_VAR, {"lny", 0, 0});			

    fob.Select(X_VAR, {					
						"Constant", 0, 0, 
						"lnx1", 0, 0, 
						"lnx2", 0, 0,
						"lnx3", 0, 0,
						"trend", 0, 0
						});	     
							
	fob.SetSelSample(-1, 1, -1, 1);
    fob.SetPrintSfa(TRUE);
    MaxControl(1000,10,TRUE);

	fob.SetTranslog(inorm);
	fob.Estimate();

	decl vEps1 = fob.Elast("lnx1");
	decl vEps2 = fob.Elast("lnx2");
	decl vEps3 = fob.Elast("lnx3");
	decl vEpst = fob.Elast("trend");

	DrawDensity(0, vEps1[][0]', {"eps1"}, 1, 1, 0, 0, 0, 0, 1, 0, 1);
	DrawDensity(1, vEps2[][0]', {"eps2"}, 1, 1, 0, 0, 0, 0, 1, 0, 1);
	DrawDensity(2, vEps3[][0]', {"eps3"}, 1, 1, 0, 0, 0, 0, 1, 0, 1);
	DrawDensity(3, vEpst[][0]', {"epst"}, 1, 1, 0, 0, 0, 0, 1, 0, 1);

	ShowDrawWindow();
	
    delete fob;
}