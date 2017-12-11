//void TensorProduct(double* A,double* B,double* C,int n);

#include<math.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>

#include "MyVonMises.cpp"


using namespace std;



double invmatrix(double *A,double *Ainv,int n);
void invmatrix2(double *A,double *Ainv,int n);
void LUdecomp(double *A,double *L,double *U,double *P,double *I,int n);
void SolveSysOEqs(double* L,double* U,double* b,double* x, int n);
void checkmat(double* A, int n);
void matvecmult(double* a, double* b, double* c, int row, int col);
void matrixmult(double* a, double* b, double* c, const int row, const int inner, const int col);
void Transpose(double *a,double *at,int row, int col);
void InitializeArray(double* A, int row, int col);
void Emodsix2three(double* De, double* Ce3);
void MyVonMises(double* sig_el,double* kap_el,double* eps_el,double* Etng_el,double* oldsig_el,double* oldkap_el,double* oldeps_el,double* Deps,double E, const double nu,double H,double sy);
void AddArray(double* A, double* B, double* C, int row, int col, int op, int id);
void ApplyBCs(double* A,double* Af,int ndofs,int nfr,int* frdofs,int id);
double CompResNorm(double* A,int n);
void GlbAsmbly(double* Kel,double* fel,double* Kg,double* fg,int* dofs_el,int ndofs);



void readmatrix(FILE* fp, double* in, const int nrows, const int ncols)
{
  for (int i=0;i<nrows;i++){
    for (int j=0;j<ncols;j++){
      double r=fscanf(fp,"%lf",&in[i*ncols+j]);
      //fp << in[i*nrows+j];
    }
  }
  fclose(fp);
}


void Write2File(double* A, int r, int c, const char* fp)
{
  std::ofstream ofile(fp);
	for(unsigned int i = 0; i < r; i++){
	  for(unsigned int j = 0; j < c; j++){
	    ofile<<A[i*c+j]<<" ";
	    //double wb=fwrite(&A[i*c+j],"%f \n",fp);
	  }
	  cout<<'\n';
	}
	ofile.close();
	//fclose(fp);
}


//void hooke(double* De, const double E, const double nu)
void hooke( const double E, const double nu)
{

  // double* De=(double*)malloc(6*6*sizeof(double));
  double fc = E/((1+nu)*(1-2*nu));
  double D_e[6][6] = {{fc*(1-nu),   fc*nu,       fc*nu,        0,            0,                0},
		     {fc*nu,       fc*(1-nu),   fc*nu,        0,            0,                0},
		     {fc*nu,       fc*nu,       fc*(1-nu),    0,            0,                0},
		     {   0,        0,           0,       fc*(1-2*nu)/2,     0,                0},
		     {   0,        0,           0,            0,      fc*(1-2*nu)/2,          0},
		     {   0,        0,           0,            0,            0,         fc*(1-2*nu)/2} };
  

  //  return De;
}


void coordxtr(double* edof,double* coords, double* ex, double* ey, int nels)
{
  //int xi, yi;
  int nod;
  
  for (int i=0;i<nels;i++){
    for (int j=0;j<3;j++){
      //xi = edof[i*7+(2*j+1)]+1;
      nod = edof[i*7+(2*j+2)] / 2;
      ex[i*3+j]=coords[(nod-1)*2+0];
      ey[i*3+j]=coords[(nod-1)*2+1];
    }
  }

}


double prescribedu(double t, double maxdisp)
{
  return maxdisp*t;
}


void checkmat(double* A, int r, int c)
{
  for(int i=0;i<r;i++){
    for(int j=0;j<c;j++){

      cout<< A[i*c+j] << " ";

       }
    cout<<'\n';
  }

}






int main()
{  


  int bcload[11][2] ={{20, 1},
		      {22, 1},
		      {34, 1},
		      {36, 1},
		      {38, 1},
		      {40, 1},
		      {42, 1},
		      {44, 1},
		      {46, 1},
		      {48, 1},
		      {50, 1}};

  int bc0[23][2]={{17, 0}, 
	     {19, 0},
	     {25, 0},
	     {27, 0},
	     {29, 0},
	     {31, 0},
	     {16, 0},
	     {24, 0},
	     {70, 0},
	     {72, 0},
	     {74, 0},
	     {76, 0},
	     {20, 1},
	     {22, 1},
	     {34, 1},
	     {36, 1},
	     {38, 1},
	     {40, 1},
	     {42, 1},
	     {44, 1},
	     {46, 1},
	     {48, 1},
	     {50, 1}};



  
  //% Rectangelar plate with a centered hole given in meters
  //% width  = 0.2        
  //% height = 0.2
  //% radius = width/4
  //load metal_sheet
  FILE *fpe = fopen("edofdata.inp","r");
  FILE *fpc = fopen("coordata0.inp","r");
  FILE *fph = fopen("hooke0.inp","r");
  FILE *fpg = fopen("gamma.inp","r");
  //FILE *fpew = fopen("edof.out","w");
  const char* fpew = "edof.out";
  
  int nels=176;
  int neldofs=6;
  int nods=108;
  double* edofs  = (double*) malloc(nels*(neldofs+1)*sizeof(double));
  double* coords = (double*) malloc(2*nods*sizeof(double));
  double* De     = (double*) malloc(6*6*sizeof(double));
  double* gamma  = (double*) malloc(9*6*sizeof(double));
  readmatrix(fpe,edofs,nels,neldofs+1);
  //checkmat(edofs,176,7);
  //Write2File(edofs, nels, neldofs+1, fpew);  
  readmatrix(fpc,coords,nods,2);
  checkmat(coords,nods,2);
  readmatrix(fph,De,6,6);
  //checkmat(De,6,6);
  readmatrix(fpg,gamma,9,6);
  //checkmat(gamma,9,6);
  
  double th    = 0.001;   // % plate thickness
  //  int ndofs  = max_element(edofs,edofs+nels*(neldofs+1));//numel(Dof);
// int nrelem = size(edof,1);
  int ndofs = 216;
  
//% Analysis parameters
double NrOfTsteps = 100.000000000000000; // % total number of timesteps
double reltol     = 1e-6; //% tolerance for equilibrium iterations
int ptype      = 1;   // % 1 = plane stress 
int magn       = 10;  // % magnification factor for deformed mesh
int maxiter    = 200;  // % max number of iterations
int iter;
 
//% Material parameters for elastic plastic material
double E          = 200*pow(10,9); 
double nu         = 0.3000000000000001;
//double nu         = 3.0/10; 
double H          = 0.2*E;
double sy         = E/(1e2);//%!!remember u changed it
//matmod     = 'plast'; % Change this to use "plast" instead of elast.
//double De[6][6];
//hooke(De,E,nu);
//hooke(E,nu);
//% Converts from gamma notation to pure epsilon 9x1
//gam = [eye(3),zeros(3);zeros(3),eye(3)/2;zeros(3),eye(3)/2];

//double gam[9][6]={ {1,0,0,0,0,0},
//                  {0,1,0,0,0,0},
//                  {0,0,1,0,0,0},
//                  {0,0,0,0.5,0,0},
//                  {0,0,0,0,0.5,0},
//                  {0,0,0,0,0,0.5},
//                  {0,0,0,0.5,0,0},
//                  {0,0,0,0,0.5,0},
//                  {0,0,0,0,0,0.5} };



//% sc = controlled stresses 
//sc = logical([0,0,1,0,1,1,0,1,1]); 
//ec = ~sc; %strain-controlled

 double* ex = (double*) malloc(3*nels*sizeof(double));
 double* ey = (double*) malloc(3*nels*sizeof(double));
//[ex,ey] = coordxtr(Edof,Coord,Dof,nen);
coordxtr(edofs,coords,ex,ey,nels);
checkmat(ex,nels,3);
checkmat(ey,nels,3); 
free(coords);
//% Solves for a linear material for reference.
//% As mentions in the project support, you may use this, or check during the main iterations.
//%Klin = spalloc(nrdof,nrdof,nen*2*nrelem);
//%for i=1:nrelem
//%    Ke = plante(ex(i,:),ey(i,:),[ptype t],De);
//%    Klin = assem(Edof(i,:),Klin,Ke);
//%end
//%u0 = solveq(Klin,zeros(nrdof,1),bc0); % Solves isotropic, elastic, problem for a unit displacement bc0.
//%u_el = ??? % Part of assignment 1.c to find this.
double u_el = 0.001; //% You can still run with just any value. Might need to adjust in order to get decent time steps.
double maxdisp = 5.0*u_el;// % Maximum displacement 5 times the elastic limit.

  //modpar.E = E;
  //modpar.H = H;
  //modpar.nu = nu;
  //modpar.sy = sy;

double* u   =(double*) malloc(ndofs*sizeof(double));//zeros(nrdof,1);
//double* sig =(double*) malloc(9*nels*sizeof(double));// zeros(9,nrelem);
//double* eps =(double*) malloc(9*nels*sizeof(double));// zeros(9,nrelem);
double* oldsig =(double*) malloc(9*nels*sizeof(double));// zeros(9,nrelem);
double* oldeps =(double*) malloc(9*nels*sizeof(double));// zeros(9,nrelem);
//double* oldeps_p =(double*) malloc(1*nels*sizeof(double));// zeros(1,nrelem);
//double* eps_p = (double*) malloc(1*nels*sizeof(double));//=oldeps_p;
double* oldkap = (double*) malloc(nels*sizeof(double));//zeros(1,nrelem);
double* kap = (double*) malloc(nels*sizeof(double));//=oldkap; 
double* u_hist    =(double*) malloc(NrOfTsteps*sizeof(double));// zeros(NrOfTsteps,1); 
double* FReaction =(double*) malloc(NrOfTsteps*sizeof(double));// zeros(NrOfTsteps,1);

double* eldisps=(double*)malloc(nels*(neldofs+1)*sizeof(double));

 InitializeArray(oldkap,1,nels);
 
 //for(int i=0;i<nels;i++){
 //	       eldisps[i*7+0]=edofs[i*7+0];
 //}
 
//prescribedu = @(t) maxdisp*t;%sin(3*2*pi*t);

//% Global timestepping
double olddisp = 0.0;
double* du = (double*) malloc(ndofs*sizeof(double));// zeros(nrdof,1);

 double dispcontrolled = 0.0;
 int nbcld=11;
 int nbc = 23;


      int* cstd =(int*) malloc(nbc*sizeof(int));
      int nfr0 = ndofs-nbc;
      int* frdofs =(int*) malloc(nfr0*sizeof(int));
      int ii = 0;
      int jj = 0;
      int nfr = 0;      
      for(int i=0;i<ndofs;i++){
	int cp = 0;
	for(int j=0;j<nbc;j++){
		if(i==bc0[j][0]-1){
		cstd[jj]=bc0[j][0]-1;
		jj++;
		cp=1;}
		//else{
		//frdofs[ii]=i;
		//ii++;
		//nfr++;}
	  }
	  if(cp==0){
	    frdofs[ii]=i;
		  ii++;
		  nfr++;}
      }


      cout<<"cstd"<<'\n';
	  for(int i=0;i<23;i++){
	    cout<<cstd[i]<<'\n';}

	  cout<<"frdofs"<<'\n';
	  for(int i=0;i<nfr;i++){
	    cout<<frdofs[i]<<'\n';}

	  
      	int ngp = 1;
	double gp1[1] = {0.333333333333333333};
	double gp2[1] = {0.333333333333333333};
	double gw[1]  = {0.50};
	// double gp1[3]={0.16666666, 0.66666666, 0.16666666};
	//double gp2[3]={0.16666666, 0.16666666, 0.66666666};
	//double gw[3]={0.16666666, 0.16666666, 0.16666666};


  
	InitializeArray(FReaction,1,NrOfTsteps);

	double err;
	double Tstepc=-1.000000000000000;
 for (int Tstep=0;Tstep<NrOfTsteps;Tstep++){
   Tstepc++;
   double prp = (Tstepc)/NrOfTsteps + 1e-17;
    dispcontrolled = prescribedu(prp,maxdisp); //% Load scalar value of prescribed di\

    //double* du =(double*) malloc(ndofs*sizeof(double));// zeros(nrdof,1);
    //du(bcload(:,1)) = dispcontrolled - olddisp;
    for (int m=0;m<ndofs;m++){
      for (int n=0;n<nbcld;n++){
	if(m==bcload[n][0]-1){
	  du[m]=dispcontrolled - olddisp;}
      }
    }


    cout<<"du"<<'\n';
	  for(int i=0;i<ndofs;i++){
	    cout<<du[i]<<'\n';}

	  
	  //double* eldisps=(double*)malloc(nels*(neldofs+1)*sizeof(double));
    
    for(int i=0;i<nels;i++){
	       eldisps[i*7+0]=edofs[i*7+0];
    }
    
    //int pos;
    //for (int k=0;k<nbcld;k++){
    //for(int i=0;i<nels;i++){
    //for(int j=1;j<7;j++){
    //if( edofs[i*7+j]==bcload[k][0]){
    //pos = edofs[i*7+j];
    //eldisps[i*7+j]=du[pos];}
    //}
    //}
    //}
    
      //printf("Timestep %f Displacement %f \n",Tstep,dispcontrolled);                                        

    double* sig =(double*) malloc(9*nels*sizeof(double));// zeros(9,nrelem);
    double* eps =(double*) malloc(9*nels*sizeof(double));// zeros(9,nrelem);
    double* kap = (double*) malloc(1*nels*sizeof(double));//=oldkap;
    
    
    for(iter=0; iter<maxiter; iter++){

      
      //if(iter>0){
      //int pos;
      for (int k=0;k<ndofs;k++){
	  for(int i=0;i<nels;i++){
	    for(int j=0;j<7;j++){
	      if( edofs[i*7+j]-1==k && j!=0){
		//pos = edofs[i*7+j];
		      eldisps[i*7+j]=du[k];}
	    }
	  }
      }
      //}
      

      cout<<"eldisps"<<'\n';
	  for(int i=0;i<nels;i++){
	    for(int j=0;j<7;j++){
	      cout<<eldisps[i*7+j]<<" ";
	    }
	    cout<<'\n';
	  }

  double* fg = (double*) malloc(ndofs*sizeof(double));
  double* Kg = (double*) malloc(ndofs*ndofs*sizeof(double));          
      InitializeArray(Kg,ndofs,ndofs);
      InitializeArray(fg,1,ndofs);
      
      for(int iel=0; iel<nels; iel++){
      // [sig(:,i),kap(:,i),eps(:,i),E_tang]=mymises(oldsig(:,i),oldkap(:,i),oldeps(:,i),Deps,modpar);%,oldeps_p(:,i)


	double* oldsig_el = (double*) malloc(9*sizeof(double));
	double* oldeps_el = (double*) malloc(9*sizeof(double));
	for(int j=0; j<9; j++){
	oldsig_el[j] = oldsig[iel*9+j];
	oldeps_el[j] = oldeps[iel*9+j];
	}

	//double oldkap_el=oldkap[iel];
	double* oldkap_el = (double*) malloc(1*sizeof(double));
	oldkap_el[0]=oldkap[iel];
	
	
	//Se = inv(De);
	//Ce = inv(Se[1 2 4],[1 2 4]);
	//trsig([1 2 4]) = Ce*;

	cout<<"eldisps"<<'\n';
	  for(int i=0;i<nels;i++){
	    for(int j=0;j<7;j++){
	      cout<<eldisps[i*7+j]<<" ";
	    }
	    cout<<'\n';
	  }
	  
	double* disp_el = (double*) malloc(neldofs*sizeof(double));
	for(int j=1;j<7;j++){
	  disp_el[j-1]=eldisps[iel*7+j];}

	cout<<"disp_el"<<'\n';
	  for(int i=0;i<6;i++){
	    cout<<disp_el[i]<<'\n';}
	  

	//double elcrds[3][2];
	double* elcrds =(double*) malloc(3*2*sizeof(double));
	for(int i=0;i<3;i++){
	  // elcrds[i][0]=ex[iel*3+i];
	  // elcrds[i][1]=ey[iel*3+i];
	    elcrds[i*2+0]=ex[iel*3+i];
	    elcrds[i*2+1]=ey[iel*3+i];
	}

	cout<<"elcrds"<<'\n';
	  for(int i=0;i<3;i++){
	    for(int j=0;j<2;j++){
	      cout<<elcrds[i*2+j]<<" ";
	    }
	    cout<<'\n';
	  }

	double ksi,etha;
	double N1,N2,N3;
	double dN1dksi, dN1detha, dN2dksi, dN2detha, dN3dksi, dN3detha;
	double* dNdksi_v =(double*) malloc(6*sizeof(double));
	double* dNdksi =(double*) malloc(2*3*sizeof(double));
	double* Jmat  =(double*) malloc(2*2*sizeof(double));
	double* Jinv  =(double*) malloc(2*2*sizeof(double));
	double* dNdx  =(double*) malloc(2*3*sizeof(double));
	double* Bmat  =(double*) malloc(3*6*sizeof(double));
	double* Bmat_v  =(double*) malloc(18*sizeof(double));

	double* Kel = (double*) malloc(6*6*sizeof(double));
        double* fel = (double*) malloc(6*sizeof(double));

	InitializeArray(Kel,6,6);
	InitializeArray(fel,1,6);
	
	double* sig_el   = (double*) malloc(9*sizeof(double));
	double* eps_el   = (double*) malloc(9*sizeof(double));
	//double* Etang_el = (double*) malloc(9*9*sizeof(double));
	//double kap_el = oldkap_el;
	double* kap_el = (double*) malloc(1*sizeof(double));;

	for (int igp=0;igp<ngp;igp++){
	  //for (int jgp=0;jgp<ngp;jgp++){

	  
	  ksi=gp1[igp];
	  etha=gp2[igp];

	  N1=ksi;
	  N2=etha;
	  N3=1-ksi-etha;

	  dN1dksi=1;
	  dN1detha=0;
	  dN2dksi=0;
	  dN2detha=1;
	  dN3dksi=-1;
	  dN3detha=-1;


	  // dNdksi={ {dN1dksi, dN2dksi, dN3dksi},
	  //	   {dN1detha, dN2detha, dN3detha} };

	  dNdksi_v[0]=dN1dksi;dNdksi_v[1]=dN2dksi;dNdksi_v[2]=dN3dksi;
	  dNdksi_v[3]=dN1detha;dNdksi_v[4]=dN2detha;dNdksi_v[5]=dN3detha;
	  //dNdksi_v={dN1dksi, dN2dksi, dN3dksi,dN1detha,dN2detha,dN3detha};

	  int cnt=0;
	  for(int i=0;i<2;i++){
	    for(int j=0;j<3;j++){
	      dNdksi[i*3+j]=dNdksi_v[cnt];
	      cnt+=1;
	    }
	  }

	  cout<<"dNdksi="<<'\n';
	  for(int i=0;i<2;i++){
	    for(int j=0;j<3;j++){
	      cout<<dNdksi[i*3+j]<<" ";
	    }
	    cout<<'\n';
	  }
	  
	  // Jmat = dNdksi*[ex(i,:)^T, ey(i,:)^T];
	  matrixmult(dNdksi,elcrds,Jmat,2,3,2);

	  cout<<"Jmat="<<'\n';
	  for(int i=0;i<2;i++){
	    for(int j=0;j<2;j++){
	      cout<<Jmat[i*2+j]<<" ";
	    }
	    cout<<'\n';
	  }

	  invmatrix2(Jmat,Jinv,2);

	  cout<<"Jinv="<<'\n';
	  for(int i=0;i<2;i++){
	    for(int j=0;j<2;j++){
	      cout<<Jinv[i*2+j]<<" ";
	    }
	    cout<<'\n';
	  }
	  
	  //dNdx=Jinv*dNdksi;
	  matrixmult(Jinv,dNdksi,dNdx,2,2,3);

	  cout<<"dNdx="<<'\n';
	  for(int i=0;i<2;i++){
	    for(int j=0;j<3;j++){
	      cout<<dNdx[i*3+j]<<" ";
	    }
	    cout<<'\n';
	  }
	  
	  // Bmat = { {dN1dx, 0, dN2dx, 0, dN3dx, 0},
	  //         {0, dN1dy, 0, dN2dy, 0, dN3dy},
	  // 	   {dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx} };

	  // Bmat_v = {dN1dx, 0, dN2dx, 0, dN3dx, 0 , 0, dN1dy, 0, dN2dy, 0, dN3dy, dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx};


 	for(int i=0;i<3;i++){
        for(int j=0;j<6;j++){
	Bmat[i*6+j]=0.0;
	}
	}


	Bmat[0]=dNdx[0];Bmat[2]=dNdx[1];Bmat[4]=dNdx[2];
	Bmat[7]=dNdx[3];Bmat[9]=dNdx[4];Bmat[11]=dNdx[5];
	Bmat[12]=dNdx[3];Bmat[14]=dNdx[4];Bmat[16]=dNdx[5];
	Bmat[13]=dNdx[0];Bmat[15]=dNdx[1];Bmat[17]=dNdx[2];
	  
	//int cnte=0;
	//int cnto=3;
	//for(int i=0;i<3;i++){
	//for(int j=0;j<6;j++){
	//if(i==0 && j%2==0){
	//Bmat[i*6+j]=dNdx[cnt];
	//cnt+=1;}
	//else if (i==1 && j%2>0){
	//Bmat[i*6+j]=dNdx[cnt];
	//cnt+=1;}
	//else if (i==2 && j%2==0){
	//Bmat[i*6+j]=dNdx[cnto];
	//cnto+=1;}
	//else if (i==2 && j%2>0){
	//Bmat[i*6+j]=dNdx[cnte];
	//cnte+=1;}
	//}
	//}


	cout<<"Bmat="<<'\n';
	  for(int i=0;i<3;i++){
	    for(int j=0;j<6;j++){
	      cout<<Bmat[i*6+j]<<" ";
	    }
	    cout<<'\n';
	  }

	  
	  //cnt=0;
	  //for(int i=0;i<3;i++){
	  //for(int j=0;j<6;j++){
	  //Bmat[i*6+j]=Bmat_v[cnt];
	  //cnt+=1;
	  //}
	  //}


	  //int row=0;
	  //int col=0;

	  double* Se = (double*) malloc(6*6*sizeof(double));
	  //double* Se3= (double*) malloc(3*3*sizeof(double));

	  invmatrix(De,Se,6);
	  
	  //Emodsix2three(Se,Se3);
	  
	  //for(int i=0;i<6;i++){
	  //  for(int j=0;j<6;j++){
	  //    if(i==0 || i==1 || i==3){
	  //	row+=1;
	  //	 if(j==0 || j==1 || j==3){
	  //	   col+=1;
	  //	   Se3[row*3+col]=Se[i*6+j];
	  //	 }
	  //	 col=0;
	  //    }
	  //  }
	  //}
	  
	  
	  double* Ce3 = (double*) malloc(3*3*sizeof(double));
	  //invmatrix(Se3,Ce3,3);
	  Emodsix2three(De,Ce3);
	  
	  
	  //Deps3 = Bmat*disp_el;
	  double* Deps3= (double*) malloc(3*sizeof(double));
	  matvecmult(Bmat,disp_el,Deps3,3,6);

	  cout<<"Deps3="<<'\n';
	  for (int i=0;i<3;i++){
	    cout<<Deps3[i]<<'\n';
	  }
	  
	  double* t_sig= (double*) malloc(6*sizeof(double));
          for(int j=0;j<6;j++){
	  t_sig[j]=0;
	  }

	  double* t_sig3= (double*) malloc(3*sizeof(double));
	  
	  //t_sig3 = Ce*Deps;
	  matvecmult(Ce3,Deps3,t_sig3,3,3);

	  cout<<"t_sig3="<<'\n';
	  for (int i=0;i<3;i++){
	    cout<<t_sig3[i]<<'\n';
	  }
	  
	  cnt=-1;
	  for(int j=0;j<6;j++){
		 if(j==0 || j==1 || j==3){
		   cnt+=1;
		   t_sig[j]=t_sig3[cnt];
		 }
	  }

	  cout<<"t_sig="<<'\n';
	  for (int i=0;i<6;i++){
	    cout<<t_sig[i]<<'\n';
	  }
	
	  double* Depsv= (double*) malloc(6*sizeof(double));	  
	  matvecmult(Se,t_sig,Depsv,6,6);

	  cout<<"Depsv="<<'\n';
	  for (int i=0;i<6;i++){
	    cout<<Depsv[i]<<'\n';
	  }
	  
	  double* Deps= (double*) malloc(9*sizeof(double));	  
	  matvecmult(gamma,Depsv,Deps,9,6);

	  cout<<"Deps="<<'\n';
	  for (int i=0;i<9;i++){
	    cout<<Deps[i]<<'\n';
	  }
	  
	  double* Etng_el = (double*) malloc(9*9*sizeof(double));
		  
	  MyVonMises(sig_el,kap_el,eps_el,Etng_el,oldsig_el,oldkap_el,oldeps_el,Deps,E,nu,H,sy);//%,oldeps_p(:,i);

	  cout<<"Etng_el="<<'\n';
          for(int i=0;i<9;i++){
            for(int j=0;j<9;j++){
              cout<<Etng_el[i*9+j]<<" ";
            }
            cout<<'\n';
          }


	  cout<<"sig_el="<<'\n';
          for (int i=0;i<9;i++){
            cout<<sig_el[i]<<'\n';
          }
	  
	  double* Etngr = (double*) malloc(6*6*sizeof(double));
	  
	  for(int i=0;i<6;i++){
	    for(int j=0;j<6;j++){
	      if(i>2 && j>2){
	  	Etngr[i*6+j] = Etng_el[i*9+j]/2.0;}
	    else {
	      Etngr[i*6+j] = Etng_el[i*9+j];}
	    }
	  }

	  

	  cout<<"Etngr  ="<<'\n';
          for(int i=0;i<6;i++){
            for(int j=0;j<6;j++){
              cout<<Etngr[i*6+j]<<" ";
            }
            cout<<'\n';
          }

	  
	  double* BmatT = (double*) malloc(6*3*sizeof(double));
	  Transpose(Bmat,BmatT,3,6);

	  cout<<"BmatT="<<'\n';
          for(int i=0;i<6;i++){
            for(int j=0;j<3;j++){
              cout<<BmatT[i*3+j]<<" ";
            }
            cout<<'\n';
          }

	  
	  double* CE3 = (double*) malloc(3*3*sizeof(double));
	  Emodsix2three(Etngr,CE3);

	  cout<<"CE3="<<'\n';
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
              cout<<CE3[i*3+j]<<" ";
            }
            cout<<'\n';
          }

	  
	  double detJ=Jmat[0]*Jmat[3]-Jmat[1]*Jmat[2];

	  
	  //Kel+ = BmatT*CE3*Bmat*detJ*gw[igp]*gw[jgp]*th;
	  double* BTC = (double*) malloc(6*3*sizeof(double));
	  matrixmult(BmatT,CE3,BTC,6,3,3);

	  cout<<"BTC="<<'\n';
          for(int i=0;i<6;i++){
            for(int j=0;j<3;j++){
              cout<<BTC[i*3+j]<<" ";
            }
            cout<<'\n';
          }
	  
	  double* BTCB = (double*) malloc(6*6*sizeof(double));
	  matrixmult(BTC,Bmat,BTCB,6,3,6);

	  cout<<"BTCB="<<'\n';
          for(int i=0;i<6;i++){
            for(int j=0;j<6;j++){
              cout<<BTCB[i*6+j]<<" ";
            }
            cout<<'\n';
          }
	  
	  for(int i=0;i<6;i++){
	    for(int j=0;j<6;j++){
	      Kel[i*6+j] = Kel[i*6+j] + BTCB[i*6+j]*detJ*gw[igp]*th;
	    }
	  }
	    //Kel+ = BTCB*detJ*gw[igp]*gw[jgp]*th;

	  cout<<"Kel="<<'\n';
          for(int i=0;i<6;i++){
            for(int j=0;j<6;j++){
              cout<<Kel[i*6+j]<<" ";
            }
            cout<<'\n';
          }

	  const char* fpkw = "Kel.out";
	  Write2File(Kel,6,6,fpkw);
	  
     	  double* sig_el3 = (double*) malloc(3*sizeof(double));
	  InitializeArray(sig_el3,1,3);
	  cnt=0;
	  for(int j=0;j<9;j++){
		 if(j==0 || j==1 || j==3){
		   sig_el3[cnt]=sig_el[j];
		   cnt+=1;
		 }
	  }

	  cout<<"sig_el3="<<'\n';
	  for (int i=0;i<3;i++){
	    cout<<sig_el3[i]<<'\n';
	  }

	  double* BTsig = (double*) malloc(6*sizeof(double));
	  InitializeArray(BTsig,1,6);
	  matvecmult(BmatT,sig_el3,BTsig,6,3);

	  cout<<"fel="<<'\n';
	  for (int i=0;i<6;i++){
	    cout<<fel[i]<<'\n';
	  }
	  
	  for(int i=0;i<6;i++){
	      fel[i] = fel[i] + BTsig[i]*detJ*gw[igp]*th;
	  }
	  //fel+ = BTsig*detJ*gw[igp]*th;
	  //fel+ = BmatT*sig_el3*detJ*gw[igp]*th;

	  cout<<"fel="<<'\n';
	  for (int i=0;i<6;i++){
	    cout<<fel[i]<<'\n';
	  }

 	free(BTsig);free(sig_el3);
	free(BTCB);free(BTC);free(CE3);free(BmatT);
	free(Etngr);free(Etng_el);free(Deps);free(Depsv);
	free(t_sig3);free(t_sig);free(Deps3);free(Ce3);free(Se);//free(Se3);
	free(oldsig_el);free(oldeps_el);free(disp_el);free(elcrds);
	free(dNdksi_v);free(dNdksi);free(Jmat);free(Jinv);free(dNdx);free(Bmat);free(Bmat_v); 
	}
	//}   //end of gp integration

	int* dofs_el = (int*) malloc(6*sizeof(int));
	int cnt=0;
	for (int j=1;j<7;j++){
	  dofs_el[cnt]=edofs[iel*7+j];
	  cnt++;
	}

	cout<<"dofs_el="<<'\n';
	  for (int i=0;i<neldofs;i++){
	    cout<<dofs_el[i]<<'\n';
	  }
	
	GlbAsmbly(Kel,fel,Kg,fg,dofs_el,ndofs);

	cout<<"fg="<<'\n';
	  for (int i=0;i<ndofs;i++){
	    cout<<fg[i]<<'\n';
	  }

	  const char* kgw = "Kg.dat";
	  Write2File(Kg,ndofs,ndofs,kgw);

	  
	
	//for(int i=0;i<6;i++){
	//for(int j=0;j<6;j++){
	//    Kg[dofs_el[i]*ndofs+dofs_el[j]] = Kg[dofs_el[i]*ndofs+dofs_el[j]] + Kel[i*6+j];
	//}
	//}

	for(int j=0; j<9; j++){
	sig[iel*9+j] = sig_el[j];
	eps[iel*9+j] = eps_el[j];
	}

	kap[iel]=kap_el[0];

	free(dofs_el);
	free(sig_el);free(eps_el);free(kap_el);
	free(Kel);free(fel);
      }  //end of element loop

      double* fext = (double*) malloc(ndofs*sizeof(double));
      InitializeArray(fext,1,ndofs);
      double* res  = (double*) malloc(ndofs*sizeof(double));
      InitializeArray(res,1,ndofs);
      double* resf = (double*) malloc(nfr*sizeof(double));
      InitializeArray(resf,1,nfr);
      
      cout<<"fext="<<'\n';
	  for (int i=0;i<ndofs;i++){
	    cout<<fext[i]<<'\n';
	  }
	  
      cout<<"fg="<<'\n';
	  for (int i=0;i<ndofs;i++){
	    cout<<fg[i]<<'\n';
	  }
      
	  AddArray(fext,fg,res,ndofs,1,1,0);
	//res = fext - fint;
      
      
      ApplyBCs(res,resf,ndofs,nfr,frdofs,0);

      cout<<"resf="<<'\n';
	  for (int i=0;i<nfr;i++){
	    cout<<resf[i]<<'\n';
	  }
      
      err = CompResNorm(resf,nfr);
      double tol = 1e-6;
      double tv1 = *std::max_element(fg,fg+ndofs);
      double tv2 = *std::min_element(fg,fg+ndofs);
      //cout<< tv <<"***********" <<fabs(tv) << endl;
      if (err<=tol*max(fabs(tv1),fabs(tv2))){
	break;}
      printf("%d %f \n",iter,err);

      double* Kgf = (double*) malloc(nfr*nfr*sizeof(double));
      ApplyBCs(Kg,Kgf,ndofs,nfr,frdofs,1);

      const char* kgfw = "Kgf.dat";
      Write2File(Kgf,nfr,nfr,kgfw);
	  
      double* L = (double*) malloc(nfr*nfr*sizeof(double));
      double* U = (double*) malloc(nfr*nfr*sizeof(double));
      double* I = (double*) malloc(nfr*nfr*sizeof(double));
      double* P = (double*) malloc(nfr*nfr*sizeof(double));

      LUdecomp(Kgf,L,U,P,I,nfr);

      double* ddu = (double*) malloc(nfr*sizeof(double));
      SolveSysOEqs(L,U,resf,ddu,nfr);      

      cout<<"ddu="<<'\n';
	  for (int i=0;i<nfr;i++){
	    cout<<ddu[i]<<'\n';
	  }

	  int pos;
      for(int i=0;i<nfr;i++){
	pos  = frdofs[i];
	du[pos] = du[pos] + ddu[i];
      }
      //AddArray(du,ddu,du,ndofs,1,1);

      cout<<"du="<<'\n';
	  for (int i=0;i<ndofs;i++){
	    cout<<du[i]<<'\n';
	  }
	  
      free(L);free(U);free(I);free(P);free(Kgf);free(Kg);free(fg);
      free(ddu); free(res);free(resf);free(fext);
    }//end of iteration loop

    AddArray(u,du,u,ndofs,1,0,0);

    for(int i=0;i<nels;i++){
      for(int j=0; j<9; j++){
	oldsig[i*9+j] = sig[i*9+j];
	oldeps[i*9+j] = eps[i*9+j];
      }
    }

    for(int i=0;i<nels;i++){
      oldkap[i]=kap[i];
    }
    //free(du);
    //free(eldisps);
    free(sig);
    free(eps);free(kap);
    olddisp = dispcontrolled;
    u_hist[Tstep]=dispcontrolled;
    int pos=0;
    for(int i=0;i<nbcld;i++){
      pos = bcload[i][0];
      //FReaction[Tstep] = FReaction[Tstep] + 2*fg[pos];
    }
    
 }//end of loadsteps

 return 0;

}

//}


    
//return 0;

 
//}


void GlbAsmbly(double* Kel,double* fel,double* Kg,double* fg,int* dofs_el,int ndofs)
{

  int ipos, jpos;
  for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){
      ipos = dofs_el[i]-1;
      jpos = dofs_el[j]-1;
      Kg[ipos      *ndofs+      jpos] = Kg[ipos      *ndofs+      jpos] + Kel[i*6+j];
      //Kg[dofs_el[i]*ndofs+dofs_el[j]] = Kg[dofs_el[i]*ndofs+dofs_el[j]] + Kel[i*6+j];
    }
  }

  //int ipos;
  for(int i=0;i<6;i++){
    ipos = dofs_el[i]-1;
    //cout<<"pos="<<pos<<"fel[i]"<<fel[i]<<endl;
	      fg[ipos] = fg[ipos] + fel[i];
  }
  
}



void InitializeArray(double* A, int row, int col)
{
  
 	for(int i=0;i<row;i++){
	  for(int j=0;j<col;j++){
	    A[i*col+j]=0.0;
	  }
	}

}


double CompResNorm(double* A,int n)
{

      double sum=0.0;
      for (int i=0;i<n;i++){
      sum+=A[i]*A[i];
      }

      return sqrt(sum);
}


void AddArray(double* A, double* B, double* C, int row, int col, int op, int id)
{

  //(if op=0 --> add)

    if( id==0){ //id=0 --> vector
      	for(int i=0;i<row;i++){
	    if(op==0){
	      C[i]=A[i]+B[i];}
	    else{
	    C[i]=A[i]-B[i];}
	}
	
    }
    else{
 	for(int i=0;i<row;i++){
	  for(int j=0;j<col;j++){
	    if(op==0){
	      C[i*col+j]=A[i*col+j]+B[i*col+j];}
	    else{
	    C[i*col+j]=A[i*col+j]-B[i*col+j];}
	  }
	}
    }
}


void ApplyBCs(double* A,double* Af,int ndofs,int nfr,int* frdofs,int id)
{

  if( id==0){ //id=0 --> vector

    int jpos;
    for(int j=0;j<nfr;j++){
            jpos = frdofs[j]; 
	    Af[j] = A[jpos];
	    //cout<<j<<"---"<<Af[j]<<endl;
    }
    
    //int cnt=0;
    //for(int i=0;i<ndofs;i++){
    //for(int j=0;j<nfr;j++){
    //if(i==frdofs[j]){
    //Af[cnt]=A[i];
    //cnt++;}
    //}
    //}

  }
  else {
    int row = -1;
    int col;
    for(int i=0;i<ndofs;i++){
      for(int k=0;k<nfr;k++){
	if(i==frdofs[k]){
	  row++;
	  col=-1;
	for(int j=0;j<ndofs;j++){
	  for(int l=0;l<nfr;l++){
	    if(j==frdofs[l]){
	      col++;
	      Af[row*nfr+col]=A[i*ndofs+j];
	    }
	  }
	}
	}
      }
    }
  }
  
}



void Emodsix2three(double* A6, double* B3)
{  
	  int row=-1;
	  int col=-1;

	  double* A_inv6 = (double*) malloc(6*6*sizeof(double));
	  double* A_inv3 = (double*) malloc(3*3*sizeof(double));

	  invmatrix(A6,A_inv6,6);

	  cout<<"A_inv6 ="<<'\n';
          checkmat(A_inv6,6);
	  
	  for(int i=0;i<6;i++){
	    if(i<4 && i!=2){
	      row+=1;
	    for(int j=0;j<6;j++){
		 if(j<4 && j!=2){
		   col+=1;
		   A_inv3[row*3+col]=A_inv6[i*6+j];
		 }
	    }
	    }
	    col=-1;
	  }

	  cout<<"A_inv3 ="<<'\n';
          checkmat(A_inv3,3);
	  
	  //double* Ce3 = (double*) malloc(3*3*sizeof(double));
	  invmatrix(A_inv3,B3,3);

	  //free(Se);
	  //free(Se3);
	  free(A_inv6);
	  free(A_inv3);

	  //return Ce3;
}


void matrixmult(double* a, double* b, double* c, const int row, const int inner, const int col){

  int i,j,k;
  double sum;

  for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			sum = 0;
			for (k = 0; k < inner; k++) {

			  //sum = sum + a[i*row + k] * b[k*row + j];
				sum = sum + a[i*inner + k] * b[k*col + j];

			}
			c[i*col + j] = sum;
		}

  }

}


void matvecmult(double* a, double* b, double* c, int row, int col)
{


  int i,j;
  for (i = 0; i < row; i++) {
     double sum = 0;
                for (j = 0; j < col; j++) {
                  //sum = 0;                                                                           \
                                                                                                        
                        //for (k = 0; k < inner; k++) {                                                \
                                                                                                        

                  //sum = sum + a[i*row + k] * b[k*row + j];                                           \
                                                                                                        
                                sum = sum + a[i*col + j] * b[j];
                                //}                                                                    \
                                                                                                        
                                //c[i*row + j] = sum;                                                  \
                                                                                                        
                }
c[i] = sum;
  }

}





void invmatrix2(double *A,double *Ainv2,int n)
{
  double detA = A[0]*A[3] - A[1]*A[2];

  //Ainv2[0]=(1/det)*A[3];
  //Ainv2[1]=(-1/det)*A[1];
  //Ainv2[2]=(-1/det)*A[2];
  //Ainv2[3]=(1/det)*A[0];
      Ainv2[0] = A[3]/detA;
      Ainv2[1] = A[1]*(-1.0)/detA;
      Ainv2[2] = A[2]*(-1.0)/detA;
      Ainv2[3] = A[0]/detA;

}


double invmatrix(double *A,double *Ainv,int n)
{


  double* LL = (double*) malloc(n*n*sizeof(double));
  double* UU = (double*) malloc(n*n*sizeof(double));
  double* II = (double*) malloc(n*n*sizeof(double));
  double* PP = (double*) malloc(n*n*sizeof(double));
  
  double* Ainv0 = (double*) malloc(n*n*sizeof(double));

  LUdecomp(A,LL,UU,PP,II,n);

          cout<<"II ="<<'\n';
          checkmat(II,n);

          
          double* I_col = (double*) malloc(n*sizeof(double));

	  
          for (int i=0;i<n;i++){

            int jj=0;
            for(int j=i;j<n*n;j+=n){
              I_col[jj]=II[j];
              cout<<j<<"*"<<I_col[jj]<<'\n';
              jj+=1;
            }


	    double* A_col = (double*) malloc(n*sizeof(double));
            // matvecmult(P,I_col,I_col,3,3);                                                           
            SolveSysOEqs(LL,UU,I_col,A_col,n);

            //jj=0;                                                                                     
            //for(int j=0;j<n;j+=n){                                                                    
            //  cout<<j<<"#"<<A_col[jj]<<'\n';                                                          
            //}                                                                                         

            jj=0;
            for(int j=i;j<n*n;j+=n){
              Ainv0[j]=A_col[jj];
              jj++;
            }
            free(A_col);

          }

          cout<<"Ainv0 ="<<'\n';
          checkmat(Ainv0,n);

	  matrixmult(PP,Ainv0,Ainv,n,n,n);

	  cout<<"Ainv ="<<'\n';
          checkmat(Ainv,n);

	  free(Ainv0);
	  free(I_col);
	  free(LL);
	  free(UU);
	  free(PP);
	  free(II);
	  
          return 0;
}



void LUdecomp(double *A,double *L,double *U,double *P,double *I,int n)
{



  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if( i==j){
        L[i*n+j]=0;
        U[i*n+j]=A[i*n+j];
        P[i*n+j]=1;
        I[i*n+j]=1;}
      else{
      L[i*n+j]=0;
       U[i*n+j]=A[i*n+j];
        P[i*n+j]=0;
        I[i*n+j]=0;
      }
    }
  }





  //************************************************************                                        

for(int row = 0; row<n; row++){


    double* Acols = (double*) malloc((n-row)*sizeof(double));
    //int vsz = (n-row)*sizeof(double);                                                                  
    //vector<double> Acols(vsz);                                                                         
    int ii=0;
    for(int i = row; i<n; i++){
       Acols[ii] = abs(A[i*n+row]);
       //cout<<ii<<"  "<<Acols[ii]<<endl;                                                               
     ii++;
    }


    //double indx = 0;                                                                                   
    //int indx = max_element(Acols,Acols+n-row)-Acols;                                                  
    int indx = distance(Acols,max_element(Acols,Acols+n-row));
    indx = indx + row;//(row-1);                                                                        

    free(Acols);
    //cout<<row << indx <<'\n';                                                                         

    double tmp;
  for (int j=0; j<n ; j++){

    //if(row>0 && j<n-1){                                                                               
    if(j<n-1){
      tmp=L[row*n+j];
    L[row*n+j]=L[indx*n+j];
    L[indx*n+j]=tmp;}

    tmp=U[row*n+j];
    U[row*n+j]=U[indx*n+j];
    U[indx*n+j]=tmp;

    tmp=P[row*n+j];
    P[row*n+j]=P[indx*n+j];
    P[indx*n+j]=tmp;

  }


    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if( i==j){
        L[i*n+j]=1;}
      }
    }


  for (int col=row+1;col<n;col++){
    L[col*n+row] = U[col*n+row]/U[row*n+row];

      for (int k=0;k<n;k++){
        U[col*n+k]=U[col*n+k]-L[col*n+row]*U[row*n+k];
      }

  }


 }

  cout<<"L ="<<'\n';
checkmat(L,n);
  cout<<"U ="<<'\n';
 checkmat(U,n);
  cout<<"P ="<<'\n';
 checkmat(P,n);


}





  //************************************************************                                       \
                                                                                                        

void SolveSysOEqs(double* L,double* U,double* b,double* x, int n)
{


  double* y = (double*) malloc(n*sizeof(double));
  //double* x = (double*) malloc(n*sizeof(double));                                                        


  //y[0] = b[0]/L[0];                                                                                   
  //double sum;                                                                                          
  //for(int i=1;i<n;i++){                                                                               
  //  sum = 0;                                                                                          
  //  for(int j=1;j<i-1;j++){                                                                           
  //    sum+=L[i*n+j]*y[j];                                                                             
  //  }                                                                                                 
  //  y[i]=(1/L[i*n+i])*(b[i]-sum);                                                                     
  // }                                                                                                  


  //x[n-1] = y[n-1]/U[n*n-1];                                                                           
  //for(int i=n-2;i>=0;i--){                                                                            
  //  sum = 0;                                                                                          
  //  for(int j=i+1;j<n;j++){                                                                           
  //    sum+=U[i*n+j]*x[j];                                                                             
  //  }                                                                                                 
  //  x[i]=(1/U[i*n+i])*(y[i]-sum);                                                                     
  // }                                                                                                  




  //y[0] = b[0]/L[0];
  double sum;
  for(int i=0;i<n;i++){
    y[i]=b[i];
    for(int j=0;j<i;j++){
      y[i]-= (L[i*n+j]*(y[j]));
    }
    //y[i]/=(L[i*n+i]);
  }


  //x[n-1] = y[n-1]/U[n*n-1];
  for(int i=n-1;i>=0;i--){
    x[i]=y[i];
    for(int j=i+1;j<n;j++){
      x[i]-=U[i*n+j]*x[j];
    }
    x[i]/=(U[i*n+i]);
  }


  //double det=Determinant(a,n);                                                                        \
                                                                                                        
  //CoFactor(a,n,b);                                                                                   \
                                                                                                        
  //return Transpose(cofmat,n);                                                                        \



}




void checkmat(double* A, int n)
{
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){

      cout<< A[i*n+j] << " ";

       }
    cout<<'\n';
  }

}












/*
   Transpose of a square matrix, do it in place
*/
//void Transpose(double **a,int n)
void Transpose(double *a,double *at,int row, int col)  
{
   int i,j;
   //double tmp;

   
   //   for (i=1;i<n;i++) {
   //   for (j=0;j<i;j++) {
   //	//tmp = a[i][j];
   //	 tmp = a[i*n+j];
   //      //a[i][j] = a[j][i];
   //	 a[i*n+j] = a[j*n+i];
   //         //a[j][i] = tmp;
   //	 a[j*n+i] = tmp;
   //      }
   //   }


   for (i=0;i<col;i++) {
      for (j=0;j<row;j++) {
	at[i*row+j]=a[j*col+i];
      }
   }
       
   //return at;
     
}























