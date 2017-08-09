#include <oxstd.h>
#import <packages/sfamb/sfamb>

main(){
    decl fob;
      // create an object of class Frontier
    fob = new Sfa();
    fob.LoadIn7("data2.in7");
      // create constant
    fob.Deterministic(FALSE);

      // Formulate the model
    fob.Select(Y_VAR, {"lny",0,0});
      // regressors
    fob.Select(X_VAR, {"Constant",0,0,
                        "lnx1",0,0,
                        "lnx2",0,0,
                        "lnx3",0,0,
                        "lnx4",0,0});
      // heteroscedastic part
   fob.Select(Z_VAR, {"Constant",0,0,
                        "lnx1",0,0});
      // technical effects, or \mu
    fob.Select(U_VAR, {"Constant",0,0,
                        "age",0,0});
      // Select estimation sample
    fob.SetSelSample(-1, 1, -1, 1);  // full sample
    fob.Estimate();
    decl mR, vR;
    mR = zeros(2,10);
    mR[0][1:4]=1;
    mR[1][1]=1;mR[1][4]=-1;
    vR=<1;0>;
    fob.TestRestrictions(mR,vR);
    delete fob;
}
