/**********Here I calculate three chain partition functions by integrating over space****************/

double q_partition(double **qB1,double **qA3,double **qC, double dr, int *Ns, double *mu, double volume){
    
    double Q,Q_AB,Q_C,Q_ABA;
    int i;
 
    Q=0.0;
    Q_AB=integrate2d_dV(qB1,0,Nr,Ns[1],dr);
    Q_C=integrate2d_dV(qA3,0,Nr,Ns[0],dr);
    Q_ABA=integrate2d_dV(qC,0,Nr,Ns[2],dr);
    
    Q_AB=exp(mu[0])*Q_AB;
    Q_ABA=exp(mu[1]*2.0)*Q_ABA/2.0;
    Q_C=(exp(mu[2]*kappa)*Q_C)/kappa;
    
    //I'm adding the three single chain partition functions together for the return function
    Q=Q_AB+Q_C+Q_ABA;
    
    // Normalizing with respect to box volume
    Q/=volume;
    
    
    return Q;
}