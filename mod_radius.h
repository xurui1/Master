void mod_radius(double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns,double ds,double *chi,double dr, double *A, double *B, double *C, int nfa){
    
    
    double fE_hom;
    double volume;
    double displacer;
    
    double *fA=create_1d_double_array(nfa,"fA");
    double *test_rad=create_1d_double_array(nfa,"test_rad");
    
    ofstream outputrad_fa;
    outputrad_fa.open("./results/outputrad_fa.dat");
    
    int counter=0;
    for (int dds=0; dds<=80;dds+=4){
        //Set parameters
        updateparameters(f,Ns,dds);
        fA[counter]=f[0];
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field
        secant(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,f,2*Nr/5);  //Find tensionless mmb
        volume=vol(dr);                                 //calculate volume
        r_0=1.0;
        double avgradius=0.0;
        
        for (int radius=0;radius<4;radius++){
            volume=vol(dr);
            omega(w);
            
            displacer=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f,2*Nr/5,0);
            int imax=mmbcentre(phi);
            avgradius+=(double)imax*dr;
            r_0*=3.0;
        }
        avgradius/=4.0;
        test_rad[counter]=avgradius;
        outputrad_fa<<f[0]<<" "<<avgradius<<endl;
        
        counter++;
        
    }
    curvefit(fA,test_rad,nfa,0,A,B,C);
    
    outputrad_fa.close();
    
    
}
