/**************My function for calculating the concentration profile for an FA = 0.5 membrane***********************************************/

void mod_phif50(double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns,double ds,double *chi,double dr){
    
    f[0] = 0.5;
    f[1] = 0.5;
    Ns[0] = 50;
    Ns[1] = 50;
    Ns[2] = 100;
    
    double volume;
    double fE_hom;
    
    //Initiate omega field
    omega(w);

    
    //Find tensionless mmb
    secant(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,f,2*Nr/5);

    
    //calculate homogeneous free energy
    fE_hom=homogfE(mu,chiMatrix,f);
    
    
    //define pinning location
    int pin = 2*Nr/5;
    
    //calculate volume
    volume=vol(dr);
    
    //reset radius
    r_0=0.5;
        
    //calculate free energy minus homogeneneous free energy
    double dFE=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f,pin,1);
    
    
        
    //output concentration profile
    ofstream outfile;
    outfile.open("./results/phi_FA50.dat");
    for (int i=0;i<Nr;i++){
        outfile<< (double)i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<endl;
    }
    
        
        
    
    
}