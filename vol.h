//Here I calculate the volume of an element in certain coordinates
double dV(int i,double dr){
    //Definition of volume elements
    double d_V;
    
    if (Coord==1){
        d_V=pow(dr,3.0);                                //Cartesian
    }
    else if (Coord==2){
        d_V=((double)i*dr+r_0)*pow(dr,3.0);             //cylindrical
    }
    else if (Coord==3){
        d_V=(pow((double)i*dr+r_0,2.0))*pow(dr,3.0);    //spherical
    }
    
    return d_V;
}

double vol(double dr){
    //Volume
    if (Coord==1){
        return ((double)Nr-1.0)*pow(dr,3.0);            //Cartesian
    }
    else if (Coord==2){
        return (1.0/2.0)*(pow(r_0+((double)Nr-1.0)*dr,2.0)-pow(r_0,2.0))*pow(dr,2.0); //cylindrical
    }
    else if (Coord==3){
        return (1.0/3.0)*(pow(r_0+((double)Nr-1.0)*dr,3.0)-pow(r_0,3.0))*pow(dr,2.0); //spherical
    }
    else{
        cout<<"Choose another integer for Coordinate"<<endl;
        return 0;
    }
    
}

double Ucellvol(double dr){
    //Volume
    if (Coord==1){
        return ((double)Nr/2.0-1.0)*pow(dr,3.0);        //Cartesian
    }
    else if (Coord==2){
        return (1.0/2.0)*(pow(r_0+((double)Nr/2.0-1.0)*dr,2.0)-pow(r_0,2.0))*pow(dr,2.0); //cylindrical
    }
    else if (Coord==3){
        return (1.0/3.0)*(pow(r_0+((double)Nr/2.0-1.0)*dr,3.0)-pow(r_0,3.0))*pow(dr,2.0); //spherical
    }
    else{
        cout<<"Choose another integer for Coordinate"<<endl;
        return 0;
    }
    
}