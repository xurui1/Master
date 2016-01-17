void phi_total(double **phi, double dr, double volume){
    
    //Here I am calculating the total concentration of each species using a trapezoidal (?) rule.
    //This is an ugly function, and I'd like to rewrite it
    
    double phiA1_tot,phiB1_tot,phiC_tot;
    double phiA2_tot,phiB2_tot,phiA3_tot;
    phiA1_tot=0.0;
    phiB1_tot=0.0;
    phiA2_tot=0.0;
    phiB2_tot=0.0;
    phiA3_tot=0.0;
    phiC_tot=0.0;
    double phi_tot=0.0;
    int i;
    
    
    //ends
    phiA1_tot+=0.5*phi[0][0]*dV(0,dr);
    phiB1_tot+=0.5*phi[1][0]*dV(0,dr);
    phiA2_tot+=0.5*phi[2][0]*dV(0,dr);
    phiB2_tot+=0.5*phi[3][0]*dV(0,dr);
    phiA3_tot+=0.5*phi[4][0]*dV(0,dr);
    phiC_tot+=0.5*phi[5][0]*dV(0,dr);
    
    phiA1_tot+=0.5*phi[0][Nr-1]*dV(Nr-1,dr);
    phiB1_tot+=0.5*phi[1][Nr-1]*dV(Nr-1,dr);
    phiA2_tot+=0.5*phi[2][Nr-1]*dV(Nr-1,dr);
    phiB2_tot+=0.5*phi[3][Nr-1]*dV(Nr-1,dr);
    phiA3_tot+=0.5*phi[4][Nr-1]*dV(Nr-1,dr);
    phiC_tot+=0.5*phi[5][Nr-1]*dV(Nr-1,dr);
    
 
    for (i=1;i<(int)Nr-1;i++){
        phiA1_tot+=phi[0][i]*dV(i,dr);
        phiB1_tot+=phi[1][i]*dV(i,dr);
        phiA2_tot+=phi[2][i]*dV(i,dr);
        phiB2_tot+=phi[3][i]*dV(i,dr);
        phiA3_tot+=phi[4][i]*dV(i,dr);
        phiC_tot+=phi[5][i]*dV(i,dr);
    }
    
    //normalize by volume
    phiA1_tot/=volume;
    phiB1_tot/=volume;
    phiA2_tot/=volume;
    phiB2_tot/=volume;
    phiA3_tot/=volume;
    phiC_tot/=volume;
    phi_tot=phiA1_tot+phiB1_tot+phiA2_tot+phiB2_tot+phiA3_tot+phiC_tot;

    
    
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
    
    double phi_ABA;
    double phi_AB;
    double OP;
    
    double phiA1_tot,phiB1_tot;
    double phiA2_tot,phiB2_tot,phiA3_tot;
    phiA1_tot=0.0;
    phiB1_tot=0.0;
    phiA2_tot=0.0;
    phiB2_tot=0.0;
    phiA3_tot=0.0;

    
    //ends
    phiA1_tot+=0.5*phi[0][0]*dV(0,dr);
    phiB1_tot+=0.5*phi[1][0]*dV(0,dr);
    phiA2_tot+=0.5*phi[2][0]*dV(0,dr);
    phiB2_tot+=0.5*phi[3][0]*dV(0,dr);
    phiA3_tot+=0.5*phi[4][0]*dV(0,dr);
    
    phiA1_tot+=0.5*phi[0][Nr-1]*dV(Nr-1,dr);
    phiB1_tot+=0.5*phi[1][Nr-1]*dV(Nr-1,dr);
    phiA2_tot+=0.5*phi[2][Nr-1]*dV(Nr-1,dr);
    phiB2_tot+=0.5*phi[3][Nr-1]*dV(Nr-1,dr);
    phiA3_tot+=0.5*phi[4][Nr-1]*dV(Nr-1,dr);
    
    //middle
    for (int i=1;i<(int)Nr-1;i++){
        phiA1_tot+=phi[0][i]*dV(i,dr);
        phiB1_tot+=phi[1][i]*dV(i,dr);
        phiA2_tot+=phi[2][i]*dV(i,dr);
        phiB2_tot+=phi[3][i]*dV(i,dr);
        phiA3_tot+=phi[4][i]*dV(i,dr);
    }
    
    //normalize by volume
    phiA1_tot/=volume;
    phiB1_tot/=volume;
    phiA2_tot/=volume;
    phiB2_tot/=volume;
    phiA3_tot/=volume;
    
    phi_AB = phiA1_tot+phiB1_tot;
    phi_ABA = phiA2_tot+phiB2_tot+phiA3_tot;
    OP = (phi_AB-phi_ABA)/(phi_AB+phi_ABA);
    
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
    
    return (outer_intersection+pin)/2;
    
}
