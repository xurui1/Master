

double calcloop(double **qA2, double **qA2L, double **qB2L,double **qA3,int *Ns,double dr,double ds, double **w, double *mu){
    
    double Q_ABA;
    double loop;
    double *p_loop;
    p_loop = create_1d_double_array(Nr, "p_loop");
    
    //Generate constrained case
    for (int i=0;i<Nr;i++){
        if (i<(Nr/2)){
            qA2L[i][Ns[0]] = qA2[i][Ns[0]];
        }
        else{
            qA2L[i][Ns[0]]=0.0;
        }
    }
    
    //solve diffusion equation
    for (int i=0;i<Nr;i++){
        qB2L[i][0]=qA2L[i][Ns[0]];
    }
    solvediffyQ(qB2L,w[3],ds,2*Ns[1],dr);
    
   //Calculate ABA chain partition function
    Q_ABA+=0.5*qA3[0][Ns[0]]/**dV(0,dr)*/;
    Q_ABA+=0.5*qA3[(int)Nr-1][Ns[0]]/**dV(Nr-1,dr)*/;
    for(int i=1;i<(int)Nr-1;i++){
        Q_ABA+=qA3[i][Ns[0]]/**dV(i,dr)*/;
    }
    Q_ABA=exp(mu[1]*2.0)*Q_ABA/2.0;
    //Q_ABA/=vol(dr);

    
    
    //Calculate probability of looping
    for (int i=0;i<Nr/2;i++){
        p_loop[i] = (qB2L[i][2*Ns[1]])*(qA2[i][Ns[0]])/Q_ABA;
        loop+=p_loop[i];
    }
   
    //loop/= Ucellvol(dr);
    
    
    destroy_1d_double_array(p_loop);
    
    return loop;

}