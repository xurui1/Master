
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
    
    double **w=create_2d_double_array(ChainType,Nr,"w");          //Auxiliary potential fields
    double *eta=create_1d_double_array(Nr,"eta");                //Incompressibility field;
    double **phi=create_2d_double_array(ChainType,Nr,"phi");      //Concentration fields
    double *chi=create_1d_double_array(ChainType,"chi");            //Interaction parameters
    double *f=create_1d_double_array(ChainType,"f");                //Chain fractions
    double *mu=create_1d_double_array(3, "mu");                     //Chemical potentials
    int *Ns=create_1d_integer_array(ChainType, "Ns");            //Chain lengths
    double **chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");
    int nradii=20,nfa=21;                                       //number of radius & fa measurements
    double *dFE=create_1d_double_array(nradii, "dFE");                  //Bending free energy
    double *Rad=create_1d_double_array(nradii, "Rad");                  //Radius for fitting
    double *Curvsq=create_1d_double_array(nradii, "Curvsq");            //Radius for fitting
    double *a1=create_1d_double_array(nfa, "a1");            //bending moduli
    double *a2=create_1d_double_array(nfa, "a2");            //bending moduli
    double *a3=create_1d_double_array(nfa, "a3");            //bending moduli
    double fE_hom;
    double avgradius;
    int imax;
    double OP,ds,dr,volume;
    
    
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
        r_0=0.1;                                        //reset radius
        avgradius=0.0;                                  //reset avgradius
    
        for (int radius=0;radius<nradii;radius++){
            volume=vol(dr);
            omega(w);
        
            dFE[radius]=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f);
            OP = calcOP(phi,dr,volume);                    //calculate order parameter
            imax=mmbcentre(phi);
            avgradius+=imax*dr;
            Rad[radius]=r_0;
            outFile2 <<f[0]<<" "<< r_0 << " "<<r_0+(double)imax*dr<<" "<<dFE[radius]<<std::endl;
            outputphi(phi,dr);
        
            r_0*=1.5;
        
        }
        avgradius/=nradii;
        for (int radius=0;radius<nradii;radius++){
            Rad[radius]+=avgradius;
            Curvsq[radius] = (4.3/Rad[radius])*(4.3/Rad[radius]);
        }
    
        
    
        curvefit(Curvsq,dFE,nradii,counter,a1,a2,a3);
        
    }
    outFile2.close();
    
    outputkappa(a1,a2,a3,nfa);

    
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
