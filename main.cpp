
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
#include "radius.h"
#include "calcexcess.h"
#include "FreeEnergy.h"
#include "curvefitting.h"
#include "mod_width.h"
#include "mod_radius.h"
#include "mod_main.h"


int main( ){
    
    double **w=create_2d_double_array(ChainType,Nr,"w");          //Auxiliary potential fields
    double *eta=create_1d_double_array(Nr,"eta");                //Incompressibility field;
    double **phi=create_2d_double_array(ChainType,Nr,"phi");      //Concentration fields
    double *chi=create_1d_double_array(ChainType,"chi");            //Interaction parameters
    double *f=create_1d_double_array(ChainType,"f");                //Chain fractions
    double *mu=create_1d_double_array(3, "mu");                     //Chemical potentials
    int *Ns=create_1d_integer_array(ChainType, "Ns");            //Chain lengths
    double **chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");
    int nradii=15,nfa=21;                                       //number of radius & fa measurements
    double *mu_vec=create_1d_double_array(nfa,"mu_vec");
    double *dFE=create_1d_double_array(nradii, "dFE");                  //Bending free energy
    double *A=create_1d_double_array(1,"A");
    double *B=create_1d_double_array(1,"B");
    double *C=create_1d_double_array(1,"C");
    double ds,dr;
    
    
    //Initial time for random number generator
    long iseed;
    time_t t;
    iseed=time(&t);
    srand48(iseed);
    
    //Set parameters & interaction matrix
    parameters(chi,f,&ds,Ns,&dr,mu);
    Xmatrix(chiMatrix,chi);
    
    //mod_width(f,mu,chiMatrix,w,phi,eta,Ns,ds,chi,dr,nfa);
    
    mod_radius(f,mu,chiMatrix,w,phi,eta,Ns,ds,chi,dr,A,B,C,nfa,mu_vec);
    
    ofstream outputrad;
    outputrad.open("./results/radius_fit.dat");
    outputrad<<A[0]<<" "<<B[0]<<" "<<C[0]<<endl;
    outputrad.close();
    
    //reset parameters
    parameters(chi,f,&ds,Ns,&dr,mu);

    //main function for finding bending moduli
    mod_main(f,mu,chiMatrix,w,phi,eta,Ns,ds,chi,dr,nfa,A,B,nradii,dFE,mu_vec);
   

    
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
