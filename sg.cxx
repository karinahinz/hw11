#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);
void step(cmplx* const psi1,cmplx* const psi0,const double omega,const int Nx, const double dx,const double dt, 
	  const double xmin);

void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
  const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt = 0.1*dx  ;
  double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
  const double omega = 0.2;
  const double alpha =sqrt(omega);

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h = new cmplx[Nx];
	
	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
		  step(psi1,psi0,omega,Nx,dx,dt,xmin);
		  h=psi0;
		  psi0=psi1;
		  psi1=h;
	
         t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
delete[] psi0;
delete[] psi1;
delete[] h;
	return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = -omega*t/2 - xi * xil * sin(omega*t);
    h3 =  0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}

void step(cmplx*psi1,cmplx* psi0,const double omega,const int Nx, const double dx,const double dt, 
	  const double xmin)
{	cmplx* d = new cmplx [Nx];
	cmplx* a= new cmplx[Nx];
	cmplx* cond = new cmplx [Nx];
	cmplx* cona= new cmplx[Nx];
	cmplx* Apsi= new cmplx[Nx];
	double x;
	
	for(int i=1;i<Nx;i++){
	x=xmin+i*dx;
	d[i]=cmplx(1.0,dt/(2*dx*dx)+dt*0.25*omega*omega*x*x);
	cond[i]=conj(d[i]);
	
	a[i]=cmplx(0.0,-dt/4*dx*dx);
	cona[i]=conj(a[i]);
	
	}
	
	Apsi[0]=cond[0]*psi0[0]+cona[1]*psi0[1];
	
	for(int i=0;i<Nx-1;i++){
	  Apsi[i]=cona[i]*psi0[i-1]+cond[i]*psi0[i]+cona[i+1]*psi0[i+1];
	}
	
	Apsi[Nx-1]=cona[Nx-2]*psi0[Nx-2]+cond[Nx-1]*psi0[Nx-1];
	
	for(int i=1;i<Nx-1;i++){
	  cmplx c=a[i]/d[i-1];
	  d[i]=d[i]-c*a[i];
	  Apsi[i]=Apsi[i]-a[i]/d[i-1]*Apsi[i-1];
	}
	
	for(int i=1;i<Nx;i++){
	  a[i]/=d[i];
	  Apsi[i]/=d[i];
	}
	psi1[Nx-1]=Apsi[Nx-1];
	
	for(int i=Nx-2;i>=0;i--){
	psi1[i]=Apsi[i]-a[i]*psi1[i+1];
	}
	delete[] d;
	delete[] a;
	delete[] cond;
	delete[] cona;
	delete[] Apsi;
}