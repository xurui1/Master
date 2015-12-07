
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
    double *Curvsq;
    double *kappaM;
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
    Rad=create_1d_double_array(20, "Rad");                  //Radius for fitting
    Curvsq=create_1d_double_array(20, "Curvsq");            //Radius for fitting
    kappaM=create_1d_double_array(20, "kappaM");
    chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");
    
    
    //Initial time for random number generator
    long iseed;
    time_t t;
    iseed=time(&t);
    srand48(iseed);
    
    //Set interaction matrix
    parameters(chi,f,&ds,Ns,&dr,mu);
    Xmatrix(chiMatrix,chi);
    
    //open main output file
    ofstream outFile2;
    string filename2;
    filename2="./results/fA_test.dat";
    outFile2.open(filename2.c_str());
    int counter=0;
    
    for (double df=0.0;df<=0.4;df+=0.02){
        counter+=1;
        //Set parameters
        updateparameters(f,Ns,df);
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field
        secant(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,f);  //Find tensionless mmb
        volume=vol(dr);                                 //calculate volume
        OP = calcOP(phi,dr,volume);                     //calculate order parameter
    
        for (radius=0;radius<20;radius++){
            volume=vol(dr);
            omega(w);
        
            dFE[radius]=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f);
            OP = calcOP(phi,dr,volume);                    //calculate order parameter
            imax=mmbcentre(phi);
            Rad[radius]=r_0+imax*dr;
            Curvsq[radius] = (4.3/Rad[radius])*(4.3/Rad[radius]);
        
            outFile2 <<f[0]<<" "<< r_0 << " "<<r_0+(double)imax*dr<<" "<<dFE[radius]<<std::endl;
            outputphi(phi,dr);
        
            r_0*=1.5;
        
        }
    
        
    
        kappaM[counter-1] = curvefit(Curvsq,dFE,20);
        
    }
    outFile2.close();
    
    outputkappa(kappaM);

    
    //Destroy memory allocations------------
    destroy_2d_double_array(w);
    destroy_1d_double_array(eta);
    destroy_2d_double_array(phi);
    destroy_1d_double_array(chi);
    destroy_1d_integer_array(Ns);
    destroy_1d_double_array(f);
    destroy_2d_double_array(chiMatrix);
    destroy_1d_double_array(Rad);
    destroy_1d_double_array(Curvsq);
    //-------------------------------------
    
    return 0;
}
