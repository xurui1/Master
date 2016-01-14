/****************************Here I set some parameters *****************************/

void parameters(double *chi,double *f,double *ds,int *Ns,double *dr,double *mu){
    
    //number of monomers
    int Ds=200;
    
    double chi_input;
    
    //inside radius
    r_0=1.0;
    
    //space step
    double delr;
    
    //initial settings
    initial=3;
    Coord=2; //if 1->Cartesian, if 2->Cylindrical, if 3->Spherical coordinate system
    
    //Length ratio of c homopolymer to diblock copolymer
    kappa=1.0;
    
    ifstream inputchi;
    inputchi.open("chi.dat");
    inputchi >> chi_input;
    inputchi.close();
        
    //Interaction parameters
    chi[0]=chi_input;        //Chi_AB
    chi[1]=chi_input;        //Chi_BC
    chi[2]=0.0;         //Chi_AC
    
    //Type of system
    //poly=0;
    
    //Chemical potential array
    if (poly==0){
        mu[0]=-20.0;      //AB
        mu[1]=0.0;      //ABA
        mu[2]=-4.48757;    //C
    }
    else if (poly==1){
        mu[0]=0.0;      //AB
        mu[1]=-20.0;      //ABA
        mu[2]=-4.48757;    //C
    }
    else{
        mu[0]=0.0;      //Other
        mu[1]=0.0;
        mu[2]=0.0;
    }
    
    //Chain length array
    Ns[0]=60;            //A blocks
    Ns[1]=140;            //B blocks
    Ns[2]=200;            //C blocks
    
    //Chain fraction array
    f[0]=Ns[0]/Ds;      //A
    f[1]=Ns[1]/Ds;      //B
    f[2]=kappa*1.0;     //C
    
    //cout<<Ns[0]<<" "<<Ns[1]<<" "<<Ns[2]<<endl;
    
    //Step size in r,z direction
    *dr=0.06;
    delr=*dr;
    
    
    //Step length along polymer
    *ds=1.0/(double)Ds;
    
}


/****************Here I build the interaction Matrix************/

void Xmatrix(double **chiMatrix, double *chi){
    //Interaction Matrix
    chiMatrix[0][0]=0.0;    //ChiA1,A1
    chiMatrix[0][1]=chi[0]; //ChiA1,B1
    chiMatrix[0][2]=0.0;    //ChiA1,A2
    chiMatrix[0][3]=chi[0]; //ChiA1,B2
    chiMatrix[0][4]=0.0;    //ChiA1,A3
    chiMatrix[0][5]=chi[2]; //ChiA1,C
    
    chiMatrix[1][0]=chi[0]; //ChiB1,A1
    chiMatrix[1][1]=0.0;    //ChiB1,B1
    chiMatrix[1][2]=chi[0]; //ChiB1,A2
    chiMatrix[1][3]=0.0;    //ChiB1,B2
    chiMatrix[1][4]=chi[0]; //ChiB1,A3
    chiMatrix[1][5]=chi[1]; //ChiB1,C
    
    chiMatrix[2][0]=0.0;    //ChiA2,A1
    chiMatrix[2][1]=chi[0]; //ChiA2,B1
    chiMatrix[2][2]=0.0;    //ChiA2,A2
    chiMatrix[2][3]=chi[0]; //ChiA2,B2
    chiMatrix[2][4]=0.0;    //ChiA2,A3
    chiMatrix[2][5]=chi[2]; //ChiA2,C
    
    chiMatrix[3][0]=chi[0]; //ChiB2,A1
    chiMatrix[3][1]=0.0;    //ChiB2,B1
    chiMatrix[3][2]=chi[0]; //ChiB2,A2
    chiMatrix[3][3]=0.0;    //ChiB2,B2
    chiMatrix[3][4]=chi[0]; //ChiB2,A3
    chiMatrix[3][5]=chi[1]; //ChiB2,C
    
    chiMatrix[4][0]=0.0;    //ChiA3,A1
    chiMatrix[4][1]=chi[0]; //ChiA3,B1
    chiMatrix[4][2]=0.0;    //ChiA3,A2
    chiMatrix[4][3]=chi[0]; //ChiA3,B2
    chiMatrix[4][4]=0.0;    //ChiA3,A3
    chiMatrix[4][5]=chi[2]; //ChiA3,C
    
    
    chiMatrix[5][0]=chi[2]; //ChiC,A1
    chiMatrix[5][1]=chi[1]; //ChiC,B1
    chiMatrix[5][2]=chi[2]; //ChiC,A2
    chiMatrix[5][3]=chi[1]; //ChiC,B2
    chiMatrix[5][4]=chi[2]; //ChiC,A3
    chiMatrix[5][2]=0.0;    //ChiC,C
    
}

/******Here I update parameters for looping through fA*******/

void updateparameters(double *f, int *Ns, int dds){
    
    double Ds=200.0;
    
    //Chain length array
    Ns[0]=60+dds;
    Ns[1]=140-dds;
    
    //Chain fraction array
    f[0]=(double)Ns[0]/Ds;    //A
    f[1]=(double)Ns[1]/Ds;  //B
    f[2]=kappa*1.0; //C

    
}
