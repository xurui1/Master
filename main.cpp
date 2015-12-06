
#include "global.h"
#include "parameters.h"
#include "filename.h"
#include "omega.h"
#include "TDMA.h"
#include "solvediffeq.h"
#include "vol.h"
#include "phi.h"
#include "Q_partition.h"
#include "polymers.h"
#include "loop.h"
#include "conc.h"
#include "Incomp.h"
#include "output.h"
#include "fE.h"
#include "homogfE.h"
#include "secant.h"
#include "FreeEnergy.h"
#include "curvefitting.h"


int main( ){
    
    double **w;
    double *eta;
    double **phi;
    double *chi;
    double *f;
    double *mu;
    double ds;
    int *Ns;
    double dr;
    double volume;
    double **chiMatrix;
    double *dFE;
    double *Rad;
    double fE_hom;
    int radius,imax;
    double OP;
    
    //Allocate memory
    w=create_2d_double_array(ChainType,Nr,"w");          //Auxiliary potential fields
    eta=create_1d_double_array(Nr,"eta");                //Incompressibility field
    phi=create_2d_double_array(ChainType,Nr,"phi");      //Concentration fields
    chi=create_1d_double_array(ChainType,"chi");            //Interaction parameters
    f=create_1d_double_array(ChainType,"f");                //Chain fractions
    Ns=create_1d_integer_array(ChainType, "Ns");            //Chain lengths
    mu=create_1d_double_array(3, "mu");                     //Chemical potentials
    dFE=create_1d_double_array(20, "dFE");                  //Bending free energy
    Rad=create_1d_double_array(20, "Rad");
    chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");
    
    
    //Initial time for random number generator
    long iseed;
    time_t t;
    iseed=time(&t);
    srand48(iseed);
    
    //Set parameters
    parameters(chi,f,&ds,Ns,&dr,mu);
    //Set interaction matrix
    Xmatrix(chiMatrix,chi);
    
    fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
    omega(w);                                       //Initiate omega field
    //secant(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,f);  //Find tensionless mmb
    volume=vol(dr);                                 //calculate volume
    OP = calcOP(phi,dr,volume);                     //calculate order parameter
    
    //open main output file
    ofstream outFile2;
    string filename2;
    filename2="./results/fe(r)_OP_" + DoubleToStr(OP)+ ".dat";
    outFile2.open(filename2.c_str());
    
    
    for (radius=0;radius<20;radius++){
        volume=vol(dr);
        omega(w);
        
        dFE[radius]=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f);
        OP = calcOP(phi,dr,volume);                    //calculate order parameter
        imax=mmbcentre(phi);
        Rad[radius]=r_0+imax*dr;

        
        outFile2 <<OP<<" "<< r_0 << " "<<r_0+(double)imax*dr<<" "<<dFE[radius]<<std::endl;
        outputphi(phi,dr);
        
        r_0*=1.2;
        
    }
    
    outFile2.close();
    
    
    curvefit(dFE,Rad,20);

    

    
    //Destroy memory allocations------------
    destroy_2d_double_array(w);
    destroy_1d_double_array(eta);
    destroy_2d_double_array(phi);
    destroy_1d_double_array(chi);
    destroy_1d_integer_array(Ns);
    destroy_1d_double_array(f);
    destroy_2d_double_array(chiMatrix);
    //-------------------------------------
    
    return 0;
}
