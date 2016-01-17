/*************This is my function for calulating the looping fraction for triblock************************/

void calcloop(double **qA2, double **qA2LL, double **qB2LL,double **qA2LR, double **qB2LR,double **qA3,int *Ns,double dr,double ds, double **w, double *mu,int imax, double *loop){
    
    
    double Q_ABA;
    double *ploopL,*ploopR;
    ploopL = create_1d_double_array(Nr, "ploopL");
    ploopR = create_1d_double_array(Nr, "ploopR");
    loop[0]=0.0;
    loop[1]=0.0;

    
    //Generate constrained case
    for (int i=0;i<Nr;i++){
        if (i<imax){
            qA2LL[i][Ns[0]] = qA2[i][Ns[0]];
            qA2LR[i][Ns[0]] = 0.0;
        }
        else{
            qA2LL[i][Ns[0]]=0.0;
            qA2LR[i][Ns[0]] = qA2[i][Ns[0]];
        }
    }
    
    
    //solve diffusion equation
    for (int i=0;i<Nr;i++){
        qB2LL[i][0]=qA2LL[i][Ns[0]];
        qB2LR[i][0]=qA2LR[i][Ns[0]];
    }
    solvediffyQ(qB2LL,w[3],ds,2*Ns[1],dr);
    solvediffyQ(qB2LR,w[3],ds,2*Ns[1],dr);
    
   //Calculate ABA chain partition function
    Q_ABA+=0.5*qA3[0][Ns[0]]/**dV(0,dr)*/;
    Q_ABA+=0.5*qA3[(int)Nr-1][Ns[0]]/**dV(Nr-1,dr)*/;
    for(int i=1;i<(int)Nr-1;i++){
        Q_ABA+=qA3[i][Ns[0]]/**dV(i,dr)*/;
    }
    Q_ABA=exp(mu[1]*2.0)*Q_ABA/2.0;
    
    
    //Calculate probability of looping
    for (int i=0;i<imax;i++){
        ploopL[i] = (qB2LL[i][2*Ns[1]])*(qA2[i][Ns[0]])/Q_ABA;
        loop[0]+=ploopL[i];
    }
    for (int i=imax;i<Nr;i++){
        ploopR[i] = (qB2LR[i][2*Ns[1]])*(qA2[i][Ns[0]])/Q_ABA;
        loop[1]+=ploopR[i];
    }
   
    
    destroy_1d_double_array(ploopL);
    destroy_1d_double_array(ploopR);
}

/*************This is my function for calulating the bridging fraction for triblock************************/

void calcbridge(double **qA2, double **qA2BL, double **qB2BL,double **qA2BR, double **qB2BR,double **qA3,int *Ns,double dr,double ds, double **w, double *mu,int imax, double *bridge){
    
    
    double Q_ABA;
    double *pbridgeL,*pbridgeR;
    pbridgeL = create_1d_double_array(Nr, "pbridgeL");
    pbridgeR = create_1d_double_array(Nr, "pbridgeR");
    bridge[0]=0.0;
    bridge[1]=0.0;
    
    
    //Generate constrained case
    for (int i=0;i<Nr;i++){
        if (i<imax){
            qA2BL[i][Ns[0]] = qA2[i][Ns[0]];
            qA2BR[i][Ns[0]] = 0.0;
        }
        else{
            qA2BL[i][Ns[0]]=0.0;
            qA2BR[i][Ns[0]] = qA2[i][Ns[0]];
        }
    }
    
    
    //solve diffusion equation
    for (int i=0;i<Nr;i++){
        qB2BL[i][0]=qA2BL[i][Ns[0]];
        qB2BR[i][0]=qA2BR[i][Ns[0]];
    }
    solvediffyQ(qB2BL,w[3],ds,2*Ns[1],dr);
    solvediffyQ(qB2BR,w[3],ds,2*Ns[1],dr);
    
    //Calculate ABA chain partition function
    Q_ABA+=0.5*qA3[0][Ns[0]]/**dV(0,dr)*/;
    Q_ABA+=0.5*qA3[(int)Nr-1][Ns[0]]/**dV(Nr-1,dr)*/;
    for(int i=1;i<(int)Nr-1;i++){
        Q_ABA+=qA3[i][Ns[0]]/**dV(i,dr)*/;
    }
    Q_ABA=exp(mu[1]*2.0)*Q_ABA/2.0;
    
    
    //Calculate probability of bridging
    for (int i=0;i<imax;i++){
        pbridgeR[i] = (qB2BR[i][2*Ns[1]])*(qA2[i][Ns[0]])/Q_ABA;
        bridge[1]+=pbridgeL[i];
    }
    for (int i=imax;i<Nr;i++){
        pbridgeR[i] = (qB2BL[i][2*Ns[1]])*(qA2[i][Ns[0]])/Q_ABA;
        bridge[0]+=pbridgeR[i];
    }
    
    
    destroy_1d_double_array(pbridgeL);
    destroy_1d_double_array(pbridgeR);
}