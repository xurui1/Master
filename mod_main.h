void mod_main(double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns,double ds,double *chi,double dr,int nfa,double *A, double *B, int nradii, double *dFE){
    
    double *r_0vector=create_1d_double_array(nradii+1, "r_0vector");
    double *Rad=create_1d_double_array(nradii, "Rad");                  //Radius for fitting
    double *Curv=create_1d_double_array(nradii, "Curvsq");            //Curv sq. for fitting
    double *Curvsq=create_1d_double_array(nradii, "Curvsq");            //Curv sq. for fitting
    double *a1=create_1d_double_array(nfa, "a1");            //bending moduli const term
    double *a2=create_1d_double_array(nfa, "a2");            //bending moduli linear term
    double *a3=create_1d_double_array(nfa, "a3");            //bending moduli quad term
    double *a4=create_1d_double_array(nfa, "a4");            //bending moduli const term
    double *a5=create_1d_double_array(nfa, "a5");            //bending moduli quad term
    double *a6=create_1d_double_array(nfa, "a6");            //bending moduli quart term
    
    //open main output file
    ofstream outFile2;
    string filename2;
    filename2="./results/fA_test.dat";
    outFile2.open(filename2.c_str());
    
    //radius ouput
    ofstream radiout;
    radiout.open("./results/main_radius.dat");
    
    int counter=0;
    double volume;
    double OP;
    double fE_hom;
    
    for (int dds=0 ;dds<=80;dds+=4){
        counter+=1;
        //Set parameters
        updateparameters(f,Ns,dds);
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field
        
        double pin_location=10.8-A[0]-B[0]*f[0];
        int pin = pin_location/dr;
        
        secant(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,f,pin);  //Find tensionless mmb
        volume=vol(dr);                                 //calculate volume
        OP = calcOP(phi,dr,volume);                     //calculate order parameter
        //Set radius vector
        set_radius(r_0vector,nradii,f[0],A,B);
        
        r_0=r_0vector[0];                                        //reset radius
        int avgradius=0.0;                                  //reset avgradius
        
        for (int radius=0;radius<nradii;radius++){
            volume=vol(dr);
            omega(w);
            

            
            dFE[radius]=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f,pin,1);
            OP = calcOP(phi,dr,volume);                    //calculate order parameter
            int imax=mmbcentre(phi);
            avgradius+=(double)imax*dr;
            Rad[radius]=r_0;
            outFile2 <<f[0]<<" "<< r_0 << " "<<r_0+(double)imax*dr<<" "<<dFE[radius]<<std::endl;
            outputphi(phi,dr);
            
            r_0=r_0vector[radius+1];
            
        }
        avgradius/=nradii;
        radiout<<f[0]<<" "<<avgradius<<endl;
        
        
        for (int radius=0;radius<nradii;radius++){
            Rad[radius]+=6.0;   //membrane should be centered
            
            Curv[radius] =(4.3/Rad[radius]);
            Curvsq[radius] = (4.3/Rad[radius])*(4.3/Rad[radius]);
        }
        
        
        
        curvefit(Curv,dFE,nradii,counter,a1,a2,a3);
        curvefit(Curvsq,dFE,nradii,counter,a4,a5,a6);
        
        
    }
    outFile2.close();
    radiout.close();
    
    outputkappa(a1,a2,a3,a4,a5,a6,nfa);
    
    
}