void mod_radius(double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns,double ds,double *chi,double dr, double *A, double *B, double *C, double *D, double *E, double *F, int nfa, double *mu_vec){
    
    
    double fE_hom;
    double volume;
    double displacer;
    
    double *fA=create_1d_double_array(nfa,"fA");
    double *test_rad1=create_1d_double_array(nfa,"test_rad1");
    double *test_rad2=create_1d_double_array(nfa,"test_rad2");

    
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
        mu_vec[counter] = mu[2];                        //Save chemical potential
        
        volume=vol(dr);                                 //calculate volume
        r_0=1.0;
        double avgradius=0.0;
        double avgmiddle=0.0;
        
        for (int radius=0;radius<4;radius++){
            volume=vol(dr);
            omega(w);
            
            displacer=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f,2*Nr/5,0);
            int imax=mmbcentre(phi);
            int ihalf=(phi,2*Nr/5);

            avgradius+=(double)imax*dr;
            avgmiddle+=(double)ihalf*dr;
            r_0*=3.0;
        }
        avgradius/=4.0;
        avgmiddle/=4.0;
        
        test_rad1[counter]=avgradius;
        test_rad2[counter]=avgmiddle;

        outputrad_fa<<f[0]<<" "<<avgradius<<endl;
        cout<<"fA: "<<f[0]<<"radius: "<<avgradius<<endl;
        
        counter++;
        
    }
    curvefit(fA,test_rad1,nfa,1,A,B,C);
    curvefit(fA,test_rad2,nfa,1,D,E,F);

    
    outputrad_fa.close();
    
    
}
