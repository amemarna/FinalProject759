#include<math.h>
#include<algorithm>
#include<stdio.h>
#include<iostream>


using namespace std;


void TensorProduct(double* A,double* B,double* C,int n);
void InitializeArray(double* A, int row, int col);
void AddArray(double* A, double* B, double* C, int row, int col, int op, int id);
void matvecmult(double* a, double* b, double* c, int row, int col);
void TensorProduct(double* A,double* B,double* C,int n);





void MyVonMises(double* sig_el,double* kap_el,double* eps_el,double* Etng_el,double* oldsig_el,double* oldkap_el,double* oldeps_el,double* Deps,double E,double\
 nu,double H,double sy)
{

  double G = E/(2*(1+nu));             // Shear modulus                                                                                                                            

  double L = E*nu/((1+nu)*(1-2*nu));   // Lame constant                                                                                                                            

  double K = L + 2*G/3;                // Bulk modulus                                                                                                                             


  double* I   = (double*) malloc(9*9*sizeof(double));
  double* I2  = (double*) malloc(9*sizeof(double));
  double* II  = (double*) malloc(9*9*sizeof(double));
  double* IIm  = (double*) malloc(9*9*sizeof(double));
  double* Idev  = (double*) malloc(9*9*sizeof(double));
  double* Ee  = (double*) malloc(9*9*sizeof(double));

                                                                
 for(int i=0;i<9;i++){
   for(int j=0;j<9;j++){
     if(i==j){
       I[i*9+j]=1.000000;}
     else{
       I[i*9+j]=0.000000;}
   }
 }

 for (int i=0;i<9;i++){
   if(i<3){
     I2[i]=1.000000;}
   else{
     I2[i]=0.000000;}
 }


                                                               
 TensorProduct(I2,I2,II,9);
                                                               


 for(int i=0;i<9;i++){
   for(int j=0;j<9;j++){
     IIm[i*9+j] = II[i*9+j]/3.0;
   }
 }


          AddArray(I,IIm,Idev,9,9,1,1);



          for(int i=0;i<9;i++){
     for(int j=0;j<9;j++){
       Ee[i*9+j] = Idev[i*9+j]*2*G + II[i*9+j]*K;
     }
   }




   for(int i=0;i<9;i++){
   eps_el[i] = oldeps_el[i] + Deps[i];
   }



   double* EeDeps  = (double*) malloc(9*sizeof(double));
   InitializeArray(EeDeps,1,9);
   matvecmult(Ee, Deps, EeDeps, 9, 9);


   double* sigtr  = (double*) malloc(9*sizeof(double));
   for (int i=0;i<9;i++){
   sigtr[i] = oldsig_el[i] + EeDeps[i];
   }
sigtr[2]=0;



double* sigtr_dev  = (double*) malloc(9*sizeof(double));
   matvecmult(Idev, sigtr, sigtr_dev, 9, 9);
   sigtr_dev[2]=0;


   double sigtr_e = 0.0;
   sigtr_e = sqrt(3.0/2.0) * sqrt(sigtr_dev[0]*sigtr_dev[0] + sigtr_dev[1]*sigtr_dev[1] + sigtr_dev[3]*sigtr_dev[3]*2.0);


   double phi_tr = sigtr_e - sy - oldkap_el[0];
   double mu,c1,sig_m;
   double* sig_dev  = (double*) malloc(9*sizeof(double));
   double* oppr  = (double*) malloc(9*9*sizeof(double));
   double* Q  = (double*) malloc(9*9*sizeof(double));


	  if (phi_tr<=0){

     for(int i=0;i<9;i++){
       for(int j=0;j<9;j++){
         Etng_el[i*9+j]=Ee[i*9+j];
       }
     }
                                                               

     kap_el[0] = oldkap_el[0];

     for (int i=0;i<9;i++){
       sig_el[i] = sigtr[i];
     }

   }

   else{
     mu = phi_tr/(3.0*G+H);
     c1 = 1.0-(3.0*mu*G/sigtr_e);

     for (int i=0;i<9;i++){
       sig_dev[i] = sigtr_dev[i]*c1;
     }

     sig_m = (sigtr[0]+sigtr[1]+sigtr[2])*(1.0/3.0);//mean trial stress                                                                                                               
     kap_el[0] = oldkap_el[0] + mu*H;                                                                                                                                
                                                               

     TensorProduct(sigtr_dev,sigtr_dev,oppr,9);
	  

     for(int i=0;i<9;i++){
       for(int j=0;j<9;j++){
         Q[i*9+j] = Idev[i*9+j] - oppr[i*9+j]*(3.0/(2.0*sigtr_e*sigtr_e));
       }
     }

     
     double bb = 3.0*G*mu/(sigtr_e+(3.0*G*mu));

     for(int i=0;i<9;i++){
       for(int j=0;j<9;j++){
         Etng_el[i*9+j] = Ee[i*9+j] - Q[i*9+j]*2.0*G*bb - oppr[i*9+j]*(9.0*G*G/((3.0*G+H)*sigtr_e*sigtr_e));
       }
     }

    for (int i=0;i<9;i++){
      sig_el[i] = sig_dev[i] + sig_m*I2[i];
    }


   }




  free(I); free(EeDeps);
  free(I2); free(sigtr);
  free(II); free(sigtr_dev);
  free(IIm);  free(sig_dev);
  free(Idev); free(oppr);
  free(Ee); free(Q);

}





void TensorProduct(double* A,double* B,double* C,int n)
{

  for(int i=0;i<n;i++){
   for(int j=0;j<n;j++){
     C[i*n+j] = A[i] * B[j];
   }
  }

}
