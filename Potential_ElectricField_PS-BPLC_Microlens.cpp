#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cfloat>
#include <ctime>
#include <cstdlib>
#include <cstring>

using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::clock;
using std::clock_t;
using std::ios;

// Macro or inline program
#define round(s) (s-floor(s)<0.5 ? int(floor(s)):int(floor(s)+1))

#define abs(x) ((x)<0 ? -(x):(x))

#define residul(x) ((2-((x)%2))>1 ? 2 : 1)

// physical constants
const double c0   = 299792458.0;              // speed of light in a vacuum in meters/second
const double c2   = c0*c0;                    // speed square of light in a vacuum in meters/second
const double mu0  = 1.256637061435917e-6;     // permeability of free space in henry/meter
const double epsf = 8.8541878176204e-12;      // permittivity of free space in farad/meter
const double pi   = 3.141592653589793;
const double Z0 = mu0*c0;                     // Impedance of free space
const double electron=1.602e-19;
const double mass=9.1e-31;
const double N_tautime=0.95e-15;
const double h_bar=(6.626e-34)/2.0/pi;

// simulation constants
const int nu = 2;
const int npmls=16;
const int nb=npmls;
const int nb1=npmls+1;
const int nb2=npmls+2;
const int nt=16384;
const double smax = 4.0e+8;
const double sqrt3=sqrt(3.0);
const double const4=-3.0*c0*(log(1e-8))/(2.0*npmls);
const double inverse_10201=1.0/double(10201);

////////////////////////////////// Functions ////////////////////////////////////

int power2(int, int);
void FileName_Choose(int);

////////////////////////////////////////////////////////////////////////////////

enum structure_type
 {
  non1,rod,air_hole
 }type;

enum form_type
 {
  non2,Triangular, Square
 }form;

enum edge_type
 {
  non3,plane,triangle
 }edge;

enum Waveguide_type
 {
  No_Defect,Single,Dual,Mach_Zehnder,Y_branch_I,Y_branch_II,DWDM,Grating
 }Waveguide;

enum Waveguide_circle
 {
  Original,move,change,change_defect
 }circle;

enum Input_type
 {
  non4,Odd_Input,Even_Input
 }input;

enum Output_type
 {
  non5,straight,slope
 }output;

enum form_lattice_type
 {
  non6,Large_a1,Large_a2
 }lattice;

enum Input_method
 {
  non7,Wave_AutoSpread,Waveguide_InputCoupling
 }coupling;

enum Perfect_type
 {
  Defect_0,Defect_nonzero,PC_separate,double_slit,Resonator,multilayer,grating,grating_PC,Wirelike,Wirelike_T,Periodic_layer
 }Perfect;

enum Grating_type
{
 Grating_Cycle_off,Grating_Cycle_on
}Grating_Cycle;

enum Defect_type
 {
  non8,single,dipole,Four,mutiple,Perioddefect_SingleWaveguide
 }Defect;

enum Random_type
 {
  Random_Defect_off,Random_Defect_on
 }Random_Defect;

enum Random_position
 {
  Random_inside,Random_outside,Random_both
 }Random_distribution;

enum Defect_reused_type
 {
  reuse_defect_No, reuse_defect_Yes
 }defect_reused;

enum Boundary_type
 {
  non9, Period_x, Period_y, Period_xy
 }Boundary;

enum eps_type
 {
  No, Yes
 }eps_reused;

enum H_type
 {
  H_No, H_Yes
 }H_record;

enum Material_type
 {
  non10, Dielectric, Metal
 }Material;

enum Metal_type
 {
  non11, Drude, Lorentz
 }Metal_eps;

enum Mode_type
 {
  non12, TE, TM, TETM
 }mode;

enum source_type
 {
  Pulse, pointwise, sine, Guassian_Beam, Plane_Wave, continuous_point
 }source;

enum source_position_type
 {
  source_No, source_Yes
 }source_default;

enum source_width_type
 {
  Width_No, Width_Yes
 }Width_default;

enum source_position_type_II
 {
  Out_No, Out_Yes
 }Out_default;

enum Fourier_Transform
 {
  Trans_No, Trans_Yes
 }Trans_default;

// ---主程式------------------------
int main()
{
  int i, j, k, m, n, z=0, ii, jj, kk, ad, i1=0, i2=0, j1, j2, nt, x1=0, y1;

  int ib, ie, jb, je, kb, ke, kp, k1, k2, NN, iix, kkz, nxmax, nymax, nzmax;

  int half_nx, top, air, core, first, second, substrate, check1=0, air2, tmax=1;

  int zmax, nz, first1, first2, core1, core2, air1, Perfect_grating, size_x_2=0;

  int substrate1, substrate2, second1, second2, count1, count2, sensor1, sensor2;

  int sensor1_time, sensor2_time, Source_Width1, Source_Width2, No_LC=0, Final_kk;            

  int Input_position, sensor3_time, sensor4_time, sensor5_time, sensor6_time;

  int y2=0, x3=0, x4=0, x5=0, x6=0, sensor5, sensor6, totalsize, enter_No, imaging_number;   

  double Refractive_Index, Epsilon, ne, no, Radius, Cellgap, V, R, Lambda, E, Es;

  double Wavelength, X_enter, E1, E2, EE, const1, dx, dy, dz, Radius2, Kerr_const;

  double const2, const3, const4, const5,const6, const7, const8, const9, x2=0, z1;

  double f, f3, alpha, beta, tan_thita, sin_thita, cos_thita, sin_phi, sin_phi2;

  double Director_Dr, Director_Dz, tan_thita0, Second_DDz_0, EE_inv, dz_step, z2;

  double coeff_k1, coeff_k2, coeff_k3, coeff_k4, coeff_l1, coeff_l2, coeff_l3;

  double sin_Gamma, cos_Gamma, tan_Gamma, Gamma, cos_thita2, cos_phi, k_initial;

  double Delta_n, Dalpha_dr, Dalpha_dz, no_inv, Kerr_sign, Dn_dr, Dn_dz, ne_inv;

  double X_enter_start, X_enter_interval, X_enter_final, Vc, const10, const11;

  char *FileNameOut="PlotFile.m";

  char *FileNameOut0="Voltage.dat";

  char *FileNameOut1="E_radius.dat";

  char *FileNameOut2="Ez.dat";

  char *FileNameOut3="Director.dat";

  char *FileNameOut4="E_radius_Dr.dat";

  char *FileNameOut5="E_radius_Dz.dat";

  char *FileNameOut6="Ez_Dr.dat";

  char *FileNameOut7="Ez_Dz.dat";

///////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Decimal point display //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

  cout.precision(15);

///////////////////////////////////////////////////////////////////////////////////
////////////////////////////// File opening operation /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

  ofstream FileOutput;
  FileOutput.open( FileNameOut, ios::out );

  FileOutput<<"clear all;                  "<<endl;
  FileOutput<<"                            "<<endl;
  
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Parameter input /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

  cout<<"                                                               "<<endl;
  cout<<"Please enter the ordinary refractive index of Blue-phase Liquid Crystal when there is no external potential: no= ";
  cin>>no;
  cout<<endl;

  cout<<"                                                               "<<endl;
  cout<<"Please enter the extra-ordinary refractive index of Blue-phase Liquid Crystal when there is no external potential: ne= ";
  cin>>ne;
  cout<<endl;

  cout<<"                                                               "<<endl;
  cout<<"Please enter the dielectric constant when Blue-phase Liquid Crystal is driven by AC voltage: ";
  cin>>Epsilon;
  cout<<endl;

  Refractive_Index=(2*no+ne)/3.0;

  cout<<"                                                               "<<endl;
  cout<<"Please enter the potential difference applied to Blue-phase Liquid Crystal: ";
  cin>>V;
  cout<<endl;

  cout<<"                                                               "<<endl;
  cout<<"Please enter the cell gap of Blue-phase Liquid Crystal (in units of micrometer): ";
  cin>>Cellgap;
  cout<<endl;

  E1=V/Cellgap;                                 // in unit of V/um 
  E2=0.0;
  const1=2.0*(Epsilon*E1+E2)/(1+Epsilon)/pi;    // It is A in Eqs. (5) and (6)!!

  cout<<"                                                               "<<endl;
  cout<<"Please enter the Kerr constant of Blue-Phase Liquid Crystal (in units of nm/VV): ";
  cin>>Kerr_const;
  cout<<endl;

  cout<<"                                                               "<<endl;
  cout<<"Please enter the saturation electric field of Blue-Phase Liquid Crystal (in units of V/um): ";
  cin>>Es;
  cout<<endl;

  cout<<"                                                               "<<endl;
  cout<<"Please enter the Threshold voltage that causes Blue-Phase Liquid Crystal to perform Kerr Effect (>=0): ";
  cin>>Vc;
  cout<<endl;

  cout<<"                                                               "<<endl;
  cout<<"Please enter the sign to make the Kerr effect of Blue-Phase Liquid Crystal produce Birefringence (1 or -1): ";
  cin>>Kerr_sign;
  cout<<endl;
  
  cout<<"                                                               "<<endl;
  cout<<"Please enter the number of additional mirror potentials required to satisfy the upper and lower plate boundary conditions: ";            
  cin>>imaging_number;
  cout<<endl; 

  cout<<"                                                               "<<endl;
  cout<<"Please enter the opening radius of the Blue-phase Liquid Crystal lens (in units of micrometer): ";
  cin>>Radius;
  cout<<endl;

  Radius2=Radius*Radius;

  cout<<"                                                               "<<endl;
  cout<<"Please enter the wavelength of the incident light (in units of nanometer): ";
  cin>>Wavelength;
  cout<<endl;

  cout<<"                                                "<<endl;
  cout<<"Please enter the grid spacing in the x direction (radial direction) (in units of micrometer): ";
  cin>>dx;
  cout<<endl;

  nxmax=round(ceil((Radius)/dx));        // nxmax=round(ceil((2*Radius)/dx+1));
  //half_nx=round(floor(0.5*nxmax));

  cout<<"                                                "<<endl;
  cout<<"Please enter the grid spacing in the z direction (radius direction) (in units of micrometer): ";
  cin>>dz;
  cout<<endl;

  nzmax=round(ceil((Cellgap)/dz));

  cout<<"no="<<no<<";                                     "<<endl;
  cout<<"ne="<<ne<<";                                     "<<endl;
  cout<<"Average Refractive Index="<<Refractive_Index<<"; "<<endl;
  cout<<"Average dielectric constant="<<Epsilon<<";       "<<endl;
  cout<<"V="<<V<<";                                       "<<endl;
  cout<<"Cellgap="<<Cellgap<<";                           "<<endl;
  cout<<"E1="<<E1<<";                                     "<<endl;
  cout<<"E2="<<E2<<";                                     "<<endl;
  cout<<"const1="<<const1<<";                             "<<endl;
  cout<<"Radius="<<Radius<<";                             "<<endl;
  cout<<"Radius2="<<Radius2<<";                           "<<endl;
  cout<<"Wavelength="<<Wavelength<<";                     "<<endl;
  cout<<"Kerr_const="<<Kerr_const<<";                     "<<endl;
  cout<<"Es="<<Es<<";                                     "<<endl;
  cout<<"Vc="<<Vc<<";                                     "<<endl;
  cout<<"imaging_number="<<imaging_number<<";             "<<endl;      
  cout<<"dx="<<dx<<";                                     "<<endl;
  cout<<"dz="<<dz<<";                                     "<<endl;
  cout<<"nxmax="<<nxmax<<";                               "<<endl;
  cout<<"nzmax="<<nzmax<<";                               "<<endl;
  cout<<"                                                 "<<endl;

  FileOutput<<"no="<<no<<";                               "<<endl;
  FileOutput<<"ne="<<ne<<";                               "<<endl;
  FileOutput<<"Refractive_Index="<<Refractive_Index<<";   "<<endl;
  FileOutput<<"Epsilon="<<Epsilon<<";                     "<<endl;
  FileOutput<<"V="<<V<<";                                 "<<endl;
  FileOutput<<"Cellgap="<<Cellgap<<";                     "<<endl;
  FileOutput<<"E1="<<E1<<";                               "<<endl;
  FileOutput<<"E2="<<E2<<";                               "<<endl;
  FileOutput<<"const1="<<const1<<";                       "<<endl;
  FileOutput<<"Radius="<<Radius<<";                       "<<endl;
  FileOutput<<"Radius2="<<Radius2<<";                     "<<endl;
  FileOutput<<"Wavelength="<<Wavelength<<";               "<<endl;
  FileOutput<<"Kerr_const="<<Kerr_const<<";               "<<endl;
  FileOutput<<"Es="<<Es<<";                               "<<endl;
  FileOutput<<"Vc="<<Vc<<";                               "<<endl;
  FileOutput<<"imaging_number="<<imaging_number<<";       "<<endl;     
  FileOutput<<"dx="<<dx<<";                               "<<endl;
  FileOutput<<"dz="<<dz<<";                               "<<endl;
  FileOutput<<"nxmax="<<nxmax<<";                         "<<endl;
  FileOutput<<"nzmax="<<nzmax<<";                         "<<endl;
  FileOutput<<"                                           "<<endl;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration Voltage ////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double **Voltage = new double *[nxmax+1];

  for(i=0; i<=nxmax; i++)
      Voltage[i] = new double [nzmax+1];

  for(i=0; i<=nxmax; i++)
      for(k=0; k<=nzmax; k++)
          Voltage[i][k] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration E_radius ///////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double **E_radius = new double *[nxmax+1];

  for(i=0; i<=nxmax; i++)
      E_radius[i] = new double [nzmax+1];

  for(i=0; i<=nxmax; i++)
      for(k=0; k<=nzmax; k++)
          E_radius[i][k] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration Ez /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double **Ez = new double *[nxmax+1];

  for(i=0; i<=nxmax; i++)
      Ez[i] = new double [nzmax+1];

  for(i=0; i<=nxmax; i++)
      for(k=0; k<=nzmax; k++)
          Ez[i][k] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration E_radius_Dr ////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double **E_radius_Dr = new double *[nxmax+1];

  for(i=0; i<=nxmax; i++)
      E_radius_Dr[i] = new double [nzmax+1];

  for(i=0; i<=nxmax; i++)
      for(k=0; k<=nzmax; k++)
          E_radius_Dr[i][k] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration E_radius_Dz ////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double **E_radius_Dz = new double *[nxmax+1];

  for(i=0; i<=nxmax; i++)
      E_radius_Dz[i] = new double [nzmax+1];

  for(i=0; i<=nxmax; i++)
      for(k=0; k<=nzmax; k++)
          E_radius_Dz[i][k] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration Ez_Dr //////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double **Ez_Dr = new double *[nxmax+1];

  for(i=0; i<=nxmax; i++)
      Ez_Dr[i] = new double [nzmax+1];

  for(i=0; i<=nxmax; i++)
      for(k=0; k<=nzmax; k++)
          Ez_Dr[i][k] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration Ez_Dz //////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double **Ez_Dz = new double *[nxmax+1];

  for(i=0; i<=nxmax; i++)
      Ez_Dz[i] = new double [nzmax+1];

  for(i=0; i<=nxmax; i++)
      for(k=0; k<=nzmax; k++)
          Ez_Dz[i][k] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration E_radius_basis /////////////
/////////////////////////////////////////////////////////////////////////////////////

  double *E_radius_basis = new double [nxmax+1];

  for(i=0; i<=nxmax; i++)
      E_radius_basis[i] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration Ez_basis ///////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double *Ez_basis = new double [nxmax+1];

  for(i=0; i<=nxmax; i++)
      Ez_basis[i] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration E_radius_Dr_basis //////////
/////////////////////////////////////////////////////////////////////////////////////

  double *E_radius_Dr_basis = new double [nxmax+1];

  for(i=0; i<=nxmax; i++)
      E_radius_Dr_basis[i] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration E_radius_Dz_basis //////////
/////////////////////////////////////////////////////////////////////////////////////

  double *E_radius_Dz_basis = new double [nxmax+1];

  for(i=0; i<=nxmax; i++)
      E_radius_Dz_basis[i] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration Ez_Dr_basis ////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double *Ez_Dr_basis = new double [nxmax+1];

  for(i=0; i<=nxmax; i++)
      Ez_Dr_basis[i] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Dynamic Memory Configuration Ez_Dz_basis ////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double *Ez_Dz_basis = new double [nxmax+1];

  for(i=0; i<=nxmax; i++)
      Ez_Dz_basis[i] = 0.0;

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Voltage Calculation ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<=nxmax; i++)
     {
      x2=(i*dx)*(i*dx);

      for(k=0; k<=nzmax; k++)
         {
          z1=(k-nzmax)*dz;           
          z2=z1*z1;

          Lambda=(x2+z2)/Radius2-1.0;
          R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

          Voltage[i][k] = V + E1*z1 -const1*(Radius*sqrt(0.5*(R-Lambda))+z1*atan(sqrt(2.0/(R+Lambda))));     

          for(m=1; m<=imaging_number; m++)                     
             {
              z1=-2*m*Cellgap-(k-nzmax)*dz;                
              z2=z1*z1;

              Lambda=(x2+z2)/Radius2-1.0;
              R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

              Voltage[i][k] = Voltage[i][k] + const1*(Radius*sqrt(0.5*(R-Lambda))+z1*atan(sqrt(2.0/(R+Lambda))));     
			 
			  if(m<imaging_number)
			    {
                 z1=(k-nzmax)*dz-2*m*Cellgap;               
                 z2=z1*z1;

                 Lambda=(x2+z2)/Radius2-1.0;
                 R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

                 Voltage[i][k] = Voltage[i][k] - const1*(Radius*sqrt(0.5*(R-Lambda))+z1*atan(sqrt(2.0/(R+Lambda))));     
				}
			 }
         }
     }

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// E_radius Calculation /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<=nxmax; i++)
     {
      x1=i*dx;
      x2=(i*dx)*(i*dx);

      for(k=0; k<=nzmax; k++)
         {
          z1=(k-nzmax)*dz;
          z2=z1*z1;

          Lambda=(x2+z2)/Radius2-1.0;
          R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

          E_radius[i][k] = -2.0*const1*(x1/(R*Radius))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0));

          for(m=1; m<=imaging_number; m++)                     
             {
              z1=-2*m*Cellgap-(k-nzmax)*dz;                
              z2=z1*z1;

              Lambda=(x2+z2)/Radius2-1.0;
              R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

              E_radius[i][k] = E_radius[i][k] + 2.0*const1*(x1/(R*Radius))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0));
			 
			  if(m<imaging_number)
			    {
                 z1=(k-nzmax)*dz-2*m*Cellgap;              
                 z2=z1*z1;

                 Lambda=(x2+z2)/Radius2-1.0;
                 R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

                 E_radius[i][k] = E_radius[i][k] - 2.0*const1*(x1/(R*Radius))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0));
				}
			 }
         }
     }

////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Ez Calculation ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<=nxmax; i++)
     {
      x2=(i*dx)*(i*dx);

      for(k=0; k<=nzmax; k++)
         {
          z1=(k-nzmax)*dz;
          z2=z1*z1;

          Lambda=(x2+z2)/Radius2-1.0;
          R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

          Ez[i][k] = -E1 + const1*(z1*(Lambda-R+2.0)/(Radius*R*sqrt(2.0*(R-Lambda)))+
                                   atan(sqrt(2.0/(R+Lambda)))-
                                   2.0*z2/(Radius2*R*sqrt(2.0*(R+Lambda))));

          for(m=1; m<=imaging_number; m++)                     
             {
              z1=-2*m*Cellgap-(k-nzmax)*dz;                
              z2=z1*z1;

              Lambda=(x2+z2)/Radius2-1.0;
              R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

              Ez[i][k] = Ez[i][k] - const1*(z1*(Lambda-R+2.0)/(Radius*R*sqrt(2.0*(R-Lambda)))+
                                            atan(sqrt(2.0/(R+Lambda)))-
                                            2.0*z2/(Radius2*R*sqrt(2.0*(R+Lambda))));
			 
			  if(m<imaging_number)
			    {
                 z1=(k-nzmax)*dz-2*m*Cellgap;              
                 z2=z1*z1;

                 Lambda=(x2+z2)/Radius2-1.0;
                 R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

                 Ez[i][k] = Ez[i][k] + const1*(z1*(Lambda-R+2.0)/(Radius*R*sqrt(2.0*(R-Lambda)))+
                                               atan(sqrt(2.0/(R+Lambda)))-
                                               2.0*z2/(Radius2*R*sqrt(2.0*(R+Lambda))));
				}
			 }
         }
     }

////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// E_radius_Dr Calculation ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<=nxmax; i++)
     {
      x2=(i*dx)*(i*dx);

      for(k=0; k<=nzmax; k++)
         {
          z1=(k-nzmax)*dz;
          z2=z1*z1;

          Lambda=(x2+z2)/Radius2-1.0;
          R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

          E_radius_Dr[i][k] = ((2.0)/(R*Radius))*((1.0-2.0*x2*Lambda/(R*R*Radius2))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0))+
                                                  (x2/Radius2)*(-0.5*sqrt(0.5*(R-Lambda))/R-(z1/Radius)*(R+Lambda-2.0)*sqrt(0.5*(R+Lambda))/(R*(R+Lambda+2.0)*(R+Lambda+2.0))));   // 1.0e+6 (不需要額外添加，因為dz單位為um，會消去此值!!)
                              
          for(m=1; m<=imaging_number; m++)                   
             {
              z1=-2*m*Cellgap-(k-nzmax)*dz;                
              z2=z1*z1;

              Lambda=(x2+z2)/Radius2-1.0;
              R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

              E_radius_Dr[i][k] = E_radius_Dr[i][k] - ((2.0)/(R*Radius))*((1.0-2.0*x2*Lambda/(R*R*Radius2))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0))+
                                                                          (x2/Radius2)*(-0.5*sqrt(0.5*(R-Lambda))/R-(z1/Radius)*(R+Lambda-2.0)*sqrt(0.5*(R+Lambda))/(R*(R+Lambda+2.0)*(R+Lambda+2.0))));   // 1.0e+6 (不需要額外添加，因為dz單位為um，會消去此值!!)
			  if(m<imaging_number)
			    {
                 z1=(k-nzmax)*dz-2*m*Cellgap;              
                 z2=z1*z1;

                 Lambda=(x2+z2)/Radius2-1.0;
                 R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

                 E_radius_Dr[i][k] = E_radius_Dr[i][k] + ((2.0)/(R*Radius))*((1.0-2.0*x2*Lambda/(R*R*Radius2))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0))+
                                                                             (x2/Radius2)*(-0.5*sqrt(0.5*(R-Lambda))/R-(z1/Radius)*(R+Lambda-2.0)*sqrt(0.5*(R+Lambda))/(R*(R+Lambda+2.0)*(R+Lambda+2.0))));   // 1.0e+6 (不需要額外添加，因為dz單位為um，會消去此值!!)
				}
			 }
         }                                                                      
     } 

////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// E_radius_Dz Calculation ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<=nxmax; i++)
     {
      x1=i*dx;
      x2=(i*dx)*(i*dx);

      for(k=0; k<=nzmax; k++)
         {
          z1=(k-nzmax)*dz;
          z2=z1*z1;

          Lambda=(x2+z2)/Radius2-1.0;
          R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

          E_radius_Dz[i][k] = ((2.0)*x1/(R*Radius2))*((-2.0*z1*(Lambda+2.0)/(R*R*Radius))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0))+
                                                        0.5*z1*(Lambda-R+2.0)/(R*Radius*sqrt(2.0*(R-Lambda)))+
                                                        sqrt(0.5*(R+Lambda))/(R+Lambda+2.0)-
                                                        (z2/Radius2)*(R+Lambda-2.0)/(R*(R+Lambda+2.0)*sqrt(2.0*(R+Lambda))));   
                              
          for(m=1; m<=imaging_number; m++)                   
             {
              z1=-2*m*Cellgap-(k-nzmax)*dz;              
              z2=z1*z1;

              Lambda=(x2+z2)/Radius2-1.0;
              R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

              E_radius_Dz[i][k] = E_radius_Dz[i][k] - ((2.0)*x1/(R*Radius2))*((-2.0*z1*(Lambda+2.0)/(R*R*Radius))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0))+
                                                                              0.5*z1*(Lambda-R+2.0)/(R*Radius*sqrt(2.0*(R-Lambda)))+
                                                                              sqrt(0.5*(R+Lambda))/(R+Lambda+2.0)-
                                                                              (z2/Radius2)*(R+Lambda-2.0)/(R*(R+Lambda+2.0)*sqrt(2.0*(R+Lambda))));
			 
			  if(m<imaging_number)
			    {
                 z1=(k-nzmax)*dz-2*m*Cellgap;            
                 z2=z1*z1;

                 Lambda=(x2+z2)/Radius2-1.0;
                 R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

                 E_radius_Dz[i][k] = E_radius_Dz[i][k] + ((2.0)*x1/(R*Radius2))*((-2.0*z1*(Lambda+2.0)/(R*R*Radius))*(0.5*sqrt(0.5*(R-Lambda))+(z1/Radius)*sqrt(0.5*(R+Lambda))/(R+Lambda+2.0))+
                                                                                 0.5*z1*(Lambda-R+2.0)/(R*Radius*sqrt(2.0*(R-Lambda)))+
                                                                                 sqrt(0.5*(R+Lambda))/(R+Lambda+2.0)-
                                                                                 (z2/Radius2)*(R+Lambda-2.0)/(R*(R+Lambda+2.0)*sqrt(2.0*(R+Lambda))));
				}
			 }
         }                                               
     }

////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Ez_Dr Calculation /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<=nxmax; i++)
     {
      x1=i*dx;
      x2=(i*dx)*(i*dx);

      for(k=0; k<=nzmax; k++)
         {
          z1=(k-nzmax)*dz;
          z2=z1*z1;

          Lambda=(x2+z2)/Radius2-1.0;
          R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

          Ez_Dr[i][k] = ((2.0)*x1/(Radius2))*(0.5*(z1/Radius)*(2.0*R*(R-Lambda)-(R-2.0*Lambda)*(R-Lambda-2.0))/(sqrt(2.0*(R-Lambda))*R*R*R)-
                                                sqrt(0.5*(R+Lambda))/(R*(R+Lambda+2.0))+
                                                z2*(2.0*Lambda+R)/(Radius2*R*R*R*sqrt(2.0*(R+Lambda))));  
                        
          for(m=1; m<=imaging_number; m++)                     
             {
              z1=-2*m*Cellgap-(k-nzmax)*dz;                
              z2=z1*z1;

              Lambda=(x2+z2)/Radius2-1.0;
              R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

              Ez_Dr[i][k] = Ez_Dr[i][k] - ((2.0)*x1/(Radius2))*(0.5*(z1/Radius)*(2.0*R*(R-Lambda)-(R-2.0*Lambda)*(R-Lambda-2.0))/(sqrt(2.0*(R-Lambda))*R*R*R)-
                                                                sqrt(0.5*(R+Lambda))/(R*(R+Lambda+2.0))+
                                                                z2*(2.0*Lambda+R)/(Radius2*R*R*R*sqrt(2.0*(R+Lambda))));
			 
			  if(m<imaging_number)
			    {
                 z1=(k-nzmax)*dz-2*m*Cellgap;              
                 z2=z1*z1;

                 Lambda=(x2+z2)/Radius2-1.0;
                 R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

                 Ez_Dr[i][k] = Ez_Dr[i][k] + ((2.0)*x1/(Radius2))*(0.5*(z1/Radius)*(2.0*R*(R-Lambda)-(R-2.0*Lambda)*(R-Lambda-2.0))/(sqrt(2.0*(R-Lambda))*R*R*R)-
                                                                   sqrt(0.5*(R+Lambda))/(R*(R+Lambda+2.0))+
                                                                   z2*(2.0*Lambda+R)/(Radius2*R*R*R*sqrt(2.0*(R+Lambda))));
				}
			 }
         }
     }

////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Ez_Dz Calculation /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<=nxmax; i++)
     {
      x2=(i*dx)*(i*dx);

      for(k=0; k<=nzmax; k++)
         {
          z1=(k-nzmax)*dz;
          z2=z1*z1;

          Lambda=(x2+z2)/Radius2-1.0;
          R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

          Ez_Dz[i][k] = ((1.0)/(Radius*R))*((Lambda-R+2.0)/sqrt(2.0*(R-Lambda))+
                                            (z2/Radius2)*(R-Lambda+2.0)*(2.0*(R-Lambda)*(R+Lambda+2.0)-R*(R-Lambda-2.0))/(sqrt(2.0*(R-Lambda)*(R-Lambda)*(R-Lambda))*R*R)-
                                            (6.0*z1)/(Radius*sqrt(2.0*(R+Lambda)))+
                                            (4.0*z1*z2/(Radius*Radius2*sqrt(2.0*(R+Lambda))))*(0.5*(2.0*Lambda+R+4.0)/(R*R)+1.0/(R*(R+Lambda))));   

          for(m=1; m<=imaging_number; m++)                     
             {
              z1=-2*m*Cellgap-(k-nzmax)*dz;                
              z2=z1*z1;

              Lambda=(x2+z2)/Radius2-1.0;
              R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

              Ez_Dz[i][k] = Ez_Dz[i][k] - ((1.0)/(Radius*R))*((Lambda-R+2.0)/sqrt(2.0*(R-Lambda))+
                                                              (z2/Radius2)*(R-Lambda+2.0)*(2.0*(R-Lambda)*(R+Lambda+2.0)-R*(R-Lambda-2.0))/(sqrt(2.0*(R-Lambda)*(R-Lambda)*(R-Lambda))*R*R)-
                                                              (6.0*z1)/(Radius*sqrt(2.0*(R+Lambda)))+
                                                              (4.0*z1*z2/(Radius*Radius2*sqrt(2.0*(R+Lambda))))*(0.5*(2.0*Lambda+R+4.0)/(R*R)+1.0/(R*(R+Lambda))));
			 
			  if(m<imaging_number)
			    {
                 z1=(k-nzmax)*dz-2*m*Cellgap;              
                 z2=z1*z1;

                 Lambda=(x2+z2)/Radius2-1.0;
                 R=sqrt(Lambda*Lambda+4.0*z2/Radius2);

                 Ez_Dz[i][k] = Ez_Dz[i][k] + ((1.0)/(Radius*R))*((Lambda-R+2.0)/sqrt(2.0*(R-Lambda))+
                                                                 (z2/Radius2)*(R-Lambda+2.0)*(2.0*(R-Lambda)*(R+Lambda+2.0)-R*(R-Lambda-2.0))/(sqrt(2.0*(R-Lambda)*(R-Lambda)*(R-Lambda))*R*R)-
                                                                 (6.0*z1)/(Radius*sqrt(2.0*(R+Lambda)))+
                                                                 (4.0*z1*z2/(Radius*Radius2*sqrt(2.0*(R+Lambda))))*(0.5*(2.0*Lambda+R+4.0)/(R*R)+1.0/(R*(R+Lambda))));
				}
			 }
         }
     }

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Record Voltage ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  if(i>=0)
    {
     ofstream FileOutput;
	 FileOutput.open(FileNameOut0, ios::out);

     for(i=0; i<=nxmax; i++)
        {
         for(k=0; k<=nzmax; k++)
            {
             if(k==nzmax)
                FileOutput<<Voltage[i][k]<<endl;
           else
                FileOutput<<Voltage[i][k]<<" ";
            }
        }

     FileOutput.close();
    }

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Record E_radius //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  if(i>=0)
    {
     ofstream FileOutput;
	 FileOutput.open(FileNameOut1, ios::out);

     for(i=0; i<=nxmax; i++)
        {
         for(k=0; k<=nzmax; k++)
            {
             if(k==nzmax)
                FileOutput<<E_radius[i][k]<<endl;
           else
                FileOutput<<E_radius[i][k]<<" ";
            }
        }

     FileOutput.close();
    }

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Record Ez ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  if(i>=0)
    {
     ofstream FileOutput;
	 FileOutput.open(FileNameOut2, ios::out);

     for(i=0; i<=nxmax; i++)
        {
         for(k=0; k<=nzmax; k++)
            {
             if(k==nzmax)
                FileOutput<<Ez[i][k]<<endl;
           else
                FileOutput<<Ez[i][k]<<" ";
            }
        }

     FileOutput.close();
    }

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Record Director //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  if(i>=0)
    {
     ofstream FileOutput;
	 FileOutput.open(FileNameOut3, ios::out);

     for(i=0; i<=nxmax; i++)
        {
         for(k=0; k<=nzmax; k++)
            {
             if(k==nzmax)
               {
                if(E_radius[i][k]==0)
                   FileOutput<<0.5*pi<<endl;
              else
                   FileOutput<<atan(Ez[i][k]/E_radius[i][k])<<endl;
               }
           else
               {
                if(E_radius[i][k]==0)
                   FileOutput<<0.5*pi<<" "<<endl;
               else
                   FileOutput<<atan(Ez[i][k]/E_radius[i][k])<<" ";
               }
            }
        }

     FileOutput.close();
    }

///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Release Memory Voltage /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

  // Release Memory Voltage
  for(i=0; i<=nxmax; i++)
      delete [] Voltage[i];
  delete [] Voltage;

///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Release Memory E_radius Ez /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

  // Release Memory E_radius
  for(i=0; i<=nxmax; i++)
      delete [] E_radius[i];
  delete [] E_radius;

  delete [] E_radius_basis; 

  // Release Memory Ez
  for(i=0; i<=nxmax; i++)
      delete [] Ez[i];
  delete [] Ez;

  delete [] Ez_basis;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Plotting & Files Saving /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

FileOutput<<"ShowDistribution=1;                                                "<<endl;
FileOutput<<"color_low=input('Please input the lowest value of the Color bar:');"<<endl;
FileOutput<<"color_top=input('Please input the topest value of the Color bar:');"<<endl;

FileOutput<<"Voltage=zeros(nxmax+1,nzmax+1);                                    "<<endl;
FileOutput<<"E_radius=zeros(nxmax+1,nzmax+1);                                   "<<endl;
FileOutput<<"Ez=zeros(nxmax+1,nzmax+1);                                         "<<endl;
FileOutput<<"Ray_radius=zeros(tmax+1,1);                                        "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"fid=fopen('"<<FileNameOut0<<"', 'r');                              "<<endl;
FileOutput<<"Voltage=fscanf(fid, '%g', [nzmax+1 nxmax+1]);                      "<<endl;
FileOutput<<"fclose(fid);                                                       "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"fid=fopen('"<<FileNameOut1<<"', 'r');                              "<<endl;
FileOutput<<"E_radius=fscanf(fid, '%g', [nzmax+1 nxmax+1]);                     "<<endl;
FileOutput<<"fclose(fid);                                                       "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"figure(3);                                                         "<<endl;
FileOutput<<"fid=fopen('"<<FileNameOut2<<"', 'r');                              "<<endl;
FileOutput<<"Ez=fscanf(fid, '%g', [nzmax+1 nxmax+1]);                           "<<endl;
FileOutput<<"fclose(fid);                                                       "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"pause;                                                             "<<endl;
FileOutput<<"count=1;                                                           "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"    clf;                                                           "<<endl;
FileOutput<<"    figure(1);                                                     "<<endl;
FileOutput<<"    pcolor(Voltage);                                               "<<endl;
FileOutput<<"    t3=['Voltage'];                                                "<<endl;
FileOutput<<"    axis('equal');                                                 "<<endl;
FileOutput<<"    axis([1 nxmax+1 1 nzmax+1]);                                   "<<endl;
FileOutput<<"    caxis([color_low color_top]);                                  "<<endl;
FileOutput<<"    shading interp;                                                "<<endl;
FileOutput<<"    colorbar;                                                      "<<endl;
FileOutput<<"    title(t3);                                                     "<<endl;
FileOutput<<"    hold;                                                          "<<endl;
FileOutput<<"    xlabel('x coordinate (um)');                                    "<<endl;
FileOutput<<"    ylabel('z coordinate (um)');                                    "<<endl;
FileOutput<<"    axis('equal');                                                  "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"    clf;                                                           "<<endl;
FileOutput<<"    figure(2);                                                     "<<endl;
FileOutput<<"    pcolor(E_radius);                                              "<<endl;
FileOutput<<"    t3=['E_radius'];                                               "<<endl;
FileOutput<<"    axis('equal');                                                 "<<endl;
FileOutput<<"    axis([1 nxmax+1 1 nzmax+1]);                                   "<<endl;
FileOutput<<"    caxis([color_low color_top]);                                  "<<endl;
FileOutput<<"    shading interp;                                                "<<endl;
FileOutput<<"    colorbar;                                                      "<<endl;
FileOutput<<"    title(t3);                                                     "<<endl;
FileOutput<<"    hold;                                                          "<<endl;
FileOutput<<"    xlabel('x coordinate (um)');                                    "<<endl;
FileOutput<<"    ylabel('z coordinate (um)');                                    "<<endl;
FileOutput<<"    axis('equal');                                                  "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"    clf;                                                           "<<endl;
FileOutput<<"    figure(3);                                                     "<<endl;
FileOutput<<"    pcolor(Ez);                                                    "<<endl;
FileOutput<<"    t3=['Ez'];                                                     "<<endl;
FileOutput<<"    axis('equal');                                                 "<<endl;
FileOutput<<"    axis([1 nxmax+1 1 nzmax+1]);                                   "<<endl;
FileOutput<<"    caxis([color_low color_top]);                                  "<<endl;
FileOutput<<"    shading interp;                                                "<<endl;
FileOutput<<"    colorbar;                                                      "<<endl;
FileOutput<<"    title(t3);                                                     "<<endl;
FileOutput<<"    hold;                                                          "<<endl;
FileOutput<<"    xlabel('x coordinate (um)');                                    "<<endl;
FileOutput<<"    ylabel('z coordinate (um)');                                    "<<endl;
FileOutput<<"    axis('equal');                                                  "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"    t3=['Equi-Voltage'];                                           "<<endl;
FileOutput<<"    V_record=zeros(11,nxmax+1);                                    "<<endl;
FileOutput<<"    count=1;                                                       "<<endl;
FileOutput<<"    for V1=90:-8:10                                                "<<endl;
FileOutput<<"        for k=1:nxmax+1                                            "<<endl;
FileOutput<<"            for i=1:nzmax                                          "<<endl;
FileOutput<<"                if Voltage(i,k)<=V1 && Voltage(i+1,k)>=V1          "<<endl;
FileOutput<<"                   V_record(count,k)=i+(V-Voltage(i,k))/(Voltage(i+1,k)-Voltage(i,k)); "<<endl;
FileOutput<<"               end                                                 "<<endl;
FileOutput<<"            end                                                    "<<endl;
FileOutput<<"        end                                                        "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"        count=count+1;                                             "<<endl;
FileOutput<<"    end                                                            "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"    figure(4)                                                      "<<endl;
FileOutput<<"    i=1:1:11                                                       "<<endl;
FileOutput<<"    plot((1:nxmax+1)*dx,V_record(i,1:nxmax+1)*dz,'b.');            "<<endl;
FileOutput<<"    hold on;                                                       "<<endl;
FileOutput<<"    title(t3);                                                     "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"                                                                   "<<endl;


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Plot Director /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

FileOutput<<"%%This area can display the Director size and distribution in Liquid Crystal!! "<<endl;
FileOutput<<"fid=fopen('"<<FileNameOut3<<"', 'r');                              "<<endl;
FileOutput<<"Directo=fscranf(fid, '%g', ["<<nxmax+1<<" "<<nzmax+1<<"]);         "<<endl;
FileOutput<<"fclose(fid);                                                       "<<endl;
FileOutput<<"    figure(9)                                                      "<<endl;
FileOutput<<"for i=1:1:"<<nxmax+1<<"                                            "<<endl;
FileOutput<<"    for k=1:"<<nzmax+1<<"                                          "<<endl;
FileOutput<<"        x=(1e-4)*E_radius(i,k);                                    "<<endl;
FileOutput<<"        z=(1e-4)*Ez(i,k);                                          "<<endl;
FileOutput<<"        quiver(i,k,x,z,'g');                                       "<<endl;
FileOutput<<"        hold on;                                                   "<<endl;
FileOutput<<"    end                                                            "<<endl;
FileOutput<<"end                                                                "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"                                                                   "<<endl;
FileOutput<<"                                                                   "<<endl;

FileOutput.close();

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Record E_radius_Dr //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

  if(i>=0)
    {
     ofstream FileOutput;
	 FileOutput.open(FileNameOut4, ios::out);

     for(i=0; i<=nxmax; i++)
        {
         for(k=0; k<=nzmax; k++)
            {
             if(k==nzmax)
                FileOutput<<E_radius_Dr[i][k]<<endl;
           else
                FileOutput<<E_radius_Dr[i][k]<<" "<<endl;
            }
        }

     FileOutput.close();
    }

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Record E_radius_Dz //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

  if(i>=0)
    {
     ofstream FileOutput;
	 FileOutput.open(FileNameOut5, ios::out);

     for(i=0; i<=nxmax; i++)
        {
         for(k=0; k<=nzmax; k++)
            {
             if(k==nzmax)
                FileOutput<<E_radius_Dz[i][k]<<endl;
           else
                FileOutput<<E_radius_Dz[i][k]<<" "<<endl;
            }
        }

     FileOutput.close();
    }

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Record Ez_Dr ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

  if(i>=0)
    {
     ofstream FileOutput;
	 FileOutput.open(FileNameOut6, ios::out);

     for(i=0; i<=nxmax; i++)
        {
         for(k=0; k<=nzmax; k++)
            {
             if(k==nzmax)
                FileOutput<<Ez_Dr[i][k]<<endl;
           else
                FileOutput<<Ez_Dr[i][k]<<" "<<endl;
            }
        }

     FileOutput.close();
    }

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Record Ez_Dz ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

  if(i>=0)
    {
     ofstream FileOutput;
	 FileOutput.open(FileNameOut7, ios::out);

     for(i=0; i<=nxmax; i++)
        {
         for(k=0; k<=nzmax; k++)
            {
             if(k==nzmax)
                FileOutput<<Ez_Dz[i][k]<<endl;
           else
                FileOutput<<Ez_Dz[i][k]<<" "<<endl;
            }
        }

     FileOutput.close();
    }

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Release Memory //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

  // Release Memory E_radius_Dr
  for(i=0; i<=nxmax; i++)
      delete [] E_radius_Dr[i];
  delete [] E_radius_Dr;

  // Release Memory E_radius_Dz
  for(i=0; i<=nxmax; i++)
      delete [] E_radius_Dz[i];
  delete [] E_radius_Dz;

  // Release Mmeory Ez_Dr
  for(i=0; i<=nxmax; i++)
      delete [] Ez_Dr[i];
  delete [] Ez_Dr;

  // Release Memory Ez_Dz
  for(i=0; i<=nxmax; i++)
      delete [] Ez_Dz[i];
  delete [] Ez_Dz;

  delete [] E_radius_Dr_basis;

  delete [] E_radius_Dz_basis;

  delete [] Ez_Dr_basis;

  delete [] Ez_Dz_basis;

/*  end=clock();
  cout<<"CLK_TCK="<<CLK_TCK<<endl;
  cout<<" The total time is "<<(end-start)/(CLK_TCK)<<" seconds:";
  cout<<endl;
*/

  cout<<"Please input any number to end this program:";
  cin>>i;
  return 0;
}

int power2(int nmax, int n)
{
 if(n>1)
    nmax=2*power2(nmax, n/2);

 return nmax;
}



