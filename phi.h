void phi_total(double **phi, double dr, double volume){
    
    //Here I am calculating the total concentration of each species using a trapezoidal (?) rule.
    //This is an ugly function, and I'd like to rewrite it
    
    double *phi_tot=create_1d_double_array(ChainType,"phi_tot");
    
    //integrate concentrations and normalize by volume
    for (int j=0;j<ChainType;j++){
        phi_tot[j] = integratedV(phi[j],0,Nr,dr);
        phi_tot[j] /= volume;
    }
    
    
    
    double total=0.0;
    
    for (int j=0;j<ChainType;j++){
        total += phi_tot[j];
    }

    //output total average concentration, if ya want!
    //cout<<total<<endl;
    
    destroy_1d_double_array(phi_tot);

}

/***************Here I calculate the various concentration profiles from the propagators******************/
void phi_calc(double **phi,double **qA1,double **qdagA1,double **qB1,double **qdagB1,double **qA2,double **qB2,double **qA3,double **qC,int *Ns,double*mu,double ds){

    
    for(int i=0;i<Nr;i++){
        
        //Empty array elements
        phi[0][i]=0.0;
        phi[1][i]=0.0;
        phi[2][i]=0.0;
        phi[3][i]=0.0;
        phi[4][i]=0.0;
        phi[5][i]=0.0;
        
        
            //phiA1 integration
            for(int s=0;s<(int)Ns[0]+1;s++){
                if(s==0 || s==(int)Ns[0]){
                    phi[0][i]+=0.5*qA1[i][s]*qdagA1[i][Ns[0]-s]*ds;
                }
                else{
                    phi[0][i]+=qA1[i][s]*qdagA1[i][Ns[0]-s]*ds;
                }
            }
            
            //phiB1 integration
            for(int s=0;s<(int)Ns[1]+1;s++){
                if(s==0 || s==(int)Ns[1]){
                    phi[1][i]+=0.5*qB1[i][s]*qdagB1[i][Ns[1]-s]*ds;
                }
                else{
                    phi[1][i]+=qB1[i][s]*qdagB1[i][Ns[1]-s]*ds;
                }
            }
        
        //phiA2 integration
        for(int s=0;s<(int)Ns[0]+1;s++){
            if(s==0 || s==(int)Ns[0]){
                phi[2][i]+=0.5*qA2[i][s]*qA3[i][Ns[0]-s]*ds;
            }
            else{
                phi[2][i]+=qA2[i][s]*qA3[i][Ns[0]-s]*ds;
            }
        }
        
        //phiB2 integration
        for(int s=0;s<2*(int)Ns[1]+1;s++){
            if(s==0 || s==2*(int)Ns[1]){
                phi[3][i]+=0.5*qB2[i][s]*qB2[i][2*Ns[1]-s]*ds;
            }
            else{
                phi[3][i]+=qB2[i][s]*qB2[i][2*Ns[1]-s]*ds;
            }
        }
        
        //phiA3 integration
        for(int s=0;s<(int)Ns[0]+1;s++){
            if(s==0 || s==(int)Ns[0]){
                phi[4][i]+=0.5*qA3[i][s]*qA2[i][Ns[0]-s]*ds;
            }
            else{
                phi[4][i]+=qA3[i][s]*qA2[i][Ns[0]-s]*ds;
            }
        }
        
        
            //phiC integration
            for(int s=0;s<(int)Ns[2]+1;s++){
                if(s==0 || s==(int)Ns[2]){
                    phi[5][i]+=0.5*qC[i][s]*qC[i][Ns[2]-s]*ds;
                }
                else{
                    phi[5][i]+=qC[i][s]*qC[i][Ns[2]-s]*ds;
                }
            }
            
            //Grand canonical relation
            phi[0][i]=exp(mu[0])*phi[0][i];
            phi[1][i]=exp(mu[0])*phi[1][i];
            phi[2][i]=exp(2.0*mu[1])*phi[2][i]/2.0;
            phi[3][i]=exp(2.0*mu[1])*phi[3][i]/2.0;
            phi[4][i]=exp(2.0*mu[1])*phi[4][i]/2.0;
            phi[5][i]=exp((mu[2])*kappa)*phi[5][i]*(1.0/kappa);
    }
    
    
}

/************Here I calculate the diblock/triblock order parameter************/
double calcOP(double **phi, double dr, double volume){
    
    //define average concentrations
    double phi_ABA;
    double phi_AB;
    double OP;
    
    double *phi_tot=create_1d_double_array(ChainType,"phi_tot");

    //integrate concentrations and normalize by volume
    for (int j=0;j<ChainType;j++){
        phi_tot[j] = integratedV(phi[j],0,Nr,dr);
        phi_tot[j] /= volume;
    }
    
    phi_AB = phi_tot[0]+phi_tot[1];
    phi_ABA = phi_tot[0]+phi_tot[1]+phi_tot[2];
    OP = (phi_AB-phi_ABA)/(phi_AB+phi_ABA);
    
    destroy_1d_double_array(phi_tot);
    
    return OP;
    
}

//calculate centre hydrophobic maximum
int mmbcentre(double **phi){
    int imax;
    double phiB1B2,phiB1B2new;
    imax=0;
    phiB1B2=phi[1][0]+phi[3][0];
    
    for (int i=0;i<Nr;i++){
        phiB1B2new=phi[1][i]+phi[3][i];
        
        if (phiB1B2new>phiB1B2){
            imax=i;
            phiB1B2=phiB1B2new;
        }
    }
    
    return imax;
    
}

//calculate right hydrophilic maximum
int mmbright(double **phi,int imax){
    int iright=imax;
    double phiA1A2A3,phiA1A2A3new;
    
    phiA1A2A3 = phi[0][imax]+phi[2][imax]+phi[4][imax];
    
    for (int i=imax;i<Nr;i++){
        phiA1A2A3new = phi[0][i]+phi[2][i]+phi[4][i];
        
        if (phiA1A2A3new>phiA1A2A3){
            iright = i;
            phiA1A2A3=phiA1A2A3new;
        }
    }
    
    return iright;
    
}

//calculate left hydrophobic maximum
int mmbleft(double **phi,int imax){
    int ileft=imax;
    double phiA1A2A3,phiA1A2A3new;
    
    phiA1A2A3 = phi[0][imax]+phi[2][imax]+phi[4][imax];
    
    for (int i=imax;i>=0;i--){
        phiA1A2A3new = phi[0][i]+phi[2][i]+phi[4][i];
        
        if (phiA1A2A3new>phiA1A2A3){
            ileft = i;
            phiA1A2A3=phiA1A2A3new;
        }
    }
    
    return ileft;
    
}

//calculate the midpoint between the two locations where phiA = phi B
int mmb_half(double **phi, int imax, int pin){
    
    //int pin is the pinning location, which is predetermined
    
    int outer_intersection = 0;
    double del_phi=1.0;
    double del_phi_new = 1.0;
    
    for (int i=imax;i<Nr;i++){
        
        del_phi_new = phi[0][i] + phi[2][i]+ phi[4][i] - phi[1][i] - phi[3][i];
        
        if (del_phi_new<del_phi){
            del_phi = del_phi_new;
            outer_intersection = i;
        }
            
    }
    
    double result = ((double)outer_intersection+(double)pin)/2.0;
    
    int output = result;
    
    return output;
    
}
