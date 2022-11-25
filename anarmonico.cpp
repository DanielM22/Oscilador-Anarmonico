#include<iostream>
#include<math.h>
#include<cmath>

using namespace std;

double h=1;
double m=1;
double k_1=1;
double k_2=0.1;

double V(double x){
  
  double V_x=(0.5*k_1*x*x)+(0.5*k_2*pow(x,4));
  
  return V_x;
}

double psiED(double psi, double x, double E){

  double ED=((2*m)/(h*h))*(V(x)-E)*psi;

  return ED;
}

int main(){
  
  double N=1000; //Pasos
  double xm=1; //x match
  double x=-6; 
  double xl=x; //DEF x left
  double xr=-x; //Def x right
  double hl=(-x+xm)/N; //Paso parte izquierda
  double hr=(x+xm)/N;  //Paso parte derecha
  double E=3.0;  //Energía inicial
  double deltaE=0.001;  //Delta de energía
  double tol=1e-8;  //Tolerancia
  double psi0L=0;   //Condición inical izquierda
  double dpsi0L=0.001;  //Condición inicial, derivada, izquierda
  double psi0R=0;  //Condición inicial derecha
  double dpsi0R=0.001;  //Condición inicial, derivada, derecha
  double Delta_b=1;  //Factor multiplicación para definición del primer delta de energía
  double k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,q1,q2,q3,q4,Delta,factor;

  for (int s=0; s<=N; s++){

    for (int i=0; i<=N; i++){

      //Método de RK-4

      //Zona izquierda
      k1=hl*dpsi0L;
      l1=hl*psiED(psi0L,xl,E);
      k2=hl*(dpsi0L+0.5*l1);
      l2=hl*psiED((psi0L+0.5*k1),(xl+0.5*hl),E);
      k3=hl*(dpsi0L+0.5*l2);
      l3=hl*psiED((psi0L+0.5*k2),(xl+0.5*hl),E);
      k4=hl*(dpsi0L+l3);
      l4=hl*psiED((psi0L+k3),(xl+hl),E);

     

      //Zona derecha
      q1=-hr*dpsi0R;
      m1=-hr*psiED(psi0R,xr,E);
      q2=-hr*(dpsi0R+0.5*m1);
      m2=-hr*psiED((psi0R+0.5*q1),(xr+0.5*hr),E);
      q3=-hr*(dpsi0R+0.5*m2);
      m3=-hr*psiED((psi0R+0.5*q2),(xr+0.5*hr),E);
      q4=-hr*(dpsi0R+m3);
      m4=-hr*psiED((psi0R+q3),(xr+hr),E);

      //cout<<xl<<"\t"<<psi0L<<"\t"<<xr<<"\t"<<psi0R<<endl;

      psi0L+=(k1+2*k2+2*k3+k4)/6;
      dpsi0L+=(l1+2*l2+2*l3+l4)/6;

      xl=xl+hl;

      psi0R+=(q1+2*q2+2*q3+q4)/6;
      dpsi0R+=(m1+2*m2+2*m3+m4)/6;

      xr=xr+hr;
      
    }//For RK

    //Nueva realización de RK para hallar psi(x+deltax)
    
    //Zona izquierda
    k1=hl*dpsi0L;
    l1=hl*psiED(psi0L,xl,E);
    k2=hl*(dpsi0L+0.5*l1);
    l2=hl*psiED((psi0L+0.5*k1),(xl+0.5*hl),E);
    k3=hl*(dpsi0L+0.5*l2);
    l3=hl*psiED((psi0L+0.5*k2),(xl+0.5*hl),E);
    k4=hl*(dpsi0L+l3);
    l4=hl*psiED((psi0L+k3),(xl+hl),E);
    
    
    
    //Zona derecha
    q1=-hr*dpsi0R;
    m1=-hr*psiED(psi0R,xr,E);
    q2=-hr*(dpsi0R+0.5*m1);
    m2=-hr*psiED((psi0R+0.5*q1),(xr+0.5*hr),E);
    q3=-hr*(dpsi0R+0.5*m2);
    m3=-hr*psiED((psi0R+0.5*q2),(xr+0.5*hr),E);
    q4=-hr*(dpsi0R+m3);
    m4=-hr*psiED((psi0R+q3),(xr+hr),E);
    
    
    double psiL=psi0L+(k1+2*k2+2*k3+k4)/6;

    double dpsiL=(psiL-psi0L)/(hl);
    
    double psiR=psi0R+(q1+2*q2+2*q3+q4)/6;

    double dpsiR=(psiR-psi0R)/(hr);

    //Función Delta (error)
    
    Delta=((dpsiL/psi0L)-(dpsiR/psi0R))/((dpsiL/psi0L)+(dpsiR/psi0R));
    factor=Delta_b*Delta;

    //Realización del método de RK-4 paral aenergía aceptada, impresión de datos
    if(abs(Delta)<=tol){
      s=N;
      psi0L=0;
      dpsi0L=0.001;
      xl=x;
      psi0R=0;
      dpsi0R=0.001;
      xr=-x;

      for (int n=0; n<=N; n++){

	//Zona izquierda
      k1=hl*dpsi0L;
      l1=hl*psiED(psi0L,xl,E);
      k2=hl*(dpsi0L+0.5*l1);
      l2=hl*psiED((psi0L+0.5*k1),(xl+0.5*hl),E);
      k3=hl*(dpsi0L+0.5*l2);
      l3=hl*psiED((psi0L+0.5*k2),(xl+0.5*hl),E);
      k4=hl*(dpsi0L+l3);
      l4=hl*psiED((psi0L+k3),(xl+hl),E);

     

      //Zona derecha
      q1=-hr*dpsi0R;
      m1=-hr*psiED(psi0R,xr,E);
      q2=-hr*(dpsi0R+0.5*m1);
      m2=-hr*psiED((psi0R+0.5*q1),(xr+0.5*hr),E);
      q3=-hr*(dpsi0R+0.5*m2);
      m3=-hr*psiED((psi0R+0.5*q2),(xr+0.5*hr),E);
      q4=-hr*(dpsi0R+m3);
      m4=-hr*psiED((psi0R+q3),(xr+hr),E);

      cout<<xl<<"\t"<<psi0L<<"\t"<<xr<<"\t"<<psi0R<<endl;

      psi0L+=(k1+2*k2+2*k3+k4)/6;
      dpsi0L+=(l1+2*l2+2*l3+l4)/6;

      xl=xl+hl;

      psi0R+=(q1+2*q2+2*q3+q4)/6;
      dpsi0R+=(m1+2*m2+2*m3+m4)/6;

      xr=xr+hr;

      }//For Datos
      
      //cout<<"Energía"<<" "<<E<<endl;
      
    }//If toletancia
    else{
      //cout<<"Datos"<<" "<<Delta<<" "<<psi0L<<" "<<psi0R<<" "<<E<<endl;

      //Redefinición de las condiciones iniciales para ejecutar el método con la nueva energía
      
      psi0L=0;
      dpsi0L=0.001;
      xl=x;
      psi0R=0;
      dpsi0R=0.001;
      xr=-x;
      
      //Actualización deltaE y de la energía
      if(factor>0){
	deltaE=deltaE;
      }
      else{
	deltaE=-deltaE/2;
      }
      Delta_b=Delta;
      E+=deltaE;
    }
    
  }//For Energía
  
  return 0;
}
