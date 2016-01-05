/*************************This is my main function for calculating free energies*************************/

double FreeEnergy(double **w, double **phi, double *eta, int *Ns, double ds, double *chi, double dr, double **chiMatrix, double *mu, double volume, double *f, int pin_location, int out_loop){
    
    
    double  currentfE, oldfE, deltafE;
    int     maxIter=100000;
    double precision=1e-5;          //convergence condition
    int     mmb;                    //Turns on pinning condition
    double  Q;                      //Chain partition functions
    double  fE_int, fES;            //interaction free energy and chain partition function fE
    double  epsilon, gamma;         //update step size
    double  *delphi;                //change in phi
    double  *sigma;
    double  *loop;
    double  **delW;
    double  **newW;
    double  deltaW;
    double fE_hom;
    fE_hom=homogfE(mu,chiMatrix,f);

    
    //Arrays for updating the omega fields
    delW=create_2d_double_array(ChainType,Nr,"delW");
    delphi=create_1d_double_array(Nr,"delphi");
    sigma=create_1d_double_array(Nr,"sigma");
    newW=create_2d_double_array(ChainType,Nr,"newW");
    loop=create_1d_double_array(2,"loop");
    
    //set energies to zero
    currentfE=0.0;
    deltafE=0.0;
    
    //update parameters
    epsilon=0.05;
    gamma=0.05;
    
    //Turn pinning condition on for membrane
    mmb=1;

        
    for (int iter=0;iter<maxIter;iter++){
        
        fE_int=0.0;
        fES=0.0;
        deltaW=0.0;

        
        Q=Conc(phi,w,Ns,ds,dr,mu,volume,loop);      //Calculate Chain partition functions
        
        
        Incomp(eta,phi,delphi);           //Enforce incompressibility condition
        
        if (iter%100==0){
            output(dr,phi);                   //Output concentration data to file
        }
        
        if (mmb==1){
            Pin(sigma, phi, pin_location);
        }
        
        //Calculate components for new field and interaction free energies
        for(int i=0;i<Nr;i++){
            for(int ii=0;ii<ChainType;ii++){
                newW[ii][i]=0.0;            //set field update to zero
                for(int jj=0;jj<ChainType;jj++){
                    newW[ii][i]+=(chiMatrix[ii][jj]*phi[jj][i]);
                }
                newW[ii][i]+=eta[i];
                
                if (mmb==1){
                    if (ii==0 || ii == 2 || ii==4){
                        newW[ii][i]-=sigma[i];
                    }
                    else if (ii==1 || ii==3){
                        newW[ii][i]+=sigma[i];
                    }
                }
                delW[ii][i]=newW[ii][i]-w[ii][i];                    //change in omega field
                w[ii][i]+=(gamma*delW[ii][i]-epsilon*delphi[i]);     //update omega field
                deltaW+=fabs(delW[ii][i])*dV(i,dr);                  //total change
                }
        }
        fE_int=fE(newW,phi,chiMatrix,dr,volume);
        
        //Normalize by box size
        deltaW/=volume;
 
        //Update free energy
        fES=Q;
        oldfE=currentfE;
        currentfE=-fES+fE_int;
        deltafE=fabs(currentfE-oldfE);
        
        //Print free energy, difference in free energy, change in omega field to screen
        if (iter%100==0){std::cout<<iter<<" fE:"<<currentfE<< " dfE:"<<currentfE-fE_hom<<" " << deltaW<<" "<<fE_hom<<" "<<loop[0]<<" "<<loop[1]<<std::endl;}
        

        if (deltafE<precision && deltaW<precision){break;} //Convergence condition
        
    }

    //output loop data if certain conditions are met
    if ((f[0] == 0.3 || f[0] == 0.4 || f[0] == 0.5 || f[0] == 0.6 || f[0] == 0.7) && out_loop==1){
        string filename="./results/loop/loopr"+IntToStr((int)10*f[0])+".dat";
        std::ofstream outputloop;
        int imax=mmbcentre(phi);
    
        outputloop.open(filename.c_str(), std::ofstream::app);
        outputloop <<4.3/(r_0+imax*dr)<<" "<<loop[0]<<" "<<loop[1]<<" "<<0.5*(loop[0]+loop[1])<<" "<<1.0-0.5*(loop[0]+loop[1])<<endl;
        outputloop.close();
    }
    
    //deallocate arrays
    destroy_1d_double_array(delphi);
    destroy_1d_double_array(loop);
    destroy_1d_double_array(sigma);
    destroy_2d_double_array(delW);
    destroy_2d_double_array(newW);
    
    return currentfE-fE_hom;
    
}
