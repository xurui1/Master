double integrate (double *function, int N, int M){
    
    double result=0.0;
    
    result+=0.5*function[N];
    result+=0.5*function[M-1];
    
    
    
    for (int i=M+1;i<N-1;i++){
        result+=function[i];
    }
    
    
    return result;
    
    
}

double integratedV (double *function, int N, int M, double dr){
    
    double result=0.0;
    
    result+=0.5*function[N]*dV(N,dr);
    result+=0.5*function[M-1]*dV(M,dr);
    
    
    
    for (int i=M+1;i<N-1;i++){
        result+=function[i]*dV(i,dr);
    }
    
    
    return result;
    
    
}


double integrate2d (double **function, int N, int M, int k){
    
    double result=0.0;
    
    result+=0.5*function[N][k];
    result+=0.5*function[M-1][k];
    
    for (int i=M+1;i<N-1;i++){
        result+=function[i][k];
    }
    
    
    return result;
    
    
}

double integrate2d_dV (double **function, int N, int M, int k, double dr){
    
    double result=0.0;
    
    result+=0.5*function[N][k]*dV(N,dr);
    result+=0.5*function[M-1][k]*dV(M,dr);

    for (int i=M+1;i<N-1;i++){
        result+=function[i][k]*dV(i,dr);
    }
    
    
    return result;
    
    
}