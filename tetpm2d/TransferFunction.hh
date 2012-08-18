#ifndef __TRANSFERFUNCTION_HH
#define __TRANSFERFUNCTION_HH

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


/* Convenience from Numerical Recipes in C, 2nd edition */
extern double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
extern double cubearg;
#define CUBE(a) ((cubearg=(a)) == 0.0 ? 0.0 : cubearg*cubearg*cubearg)
extern double pow4arg;
#define POW4(a) ((pow4arg=(a)) == 0.0 ? 0.0 : pow4arg*pow4arg*pow4arg*pow4arg)

extern double
g_H0,       // Hubble constant
g_Omega_b,      // baryon content
g_Omega_m,//0.276,       // matter content
g_Omega_L,//0.724,       // dark energy content
g_sigma_8,//0.811,        // spectrum normalization
g_nspect,//0.961,        // spectral index
g_zstart;


extern const double
g_deltac,      // critical density for spherical collapse
g_G, // gravitational constant in Mpc km^2/s^2 1.0e+10Msun
g_rhoc; // critical density




//! Abstract base class for transfer functions
/*!
    This class implements a purely virtual interface that can be
    used to derive instances implementing various transfer functions.
*/ 
class TransferFunction{
public:
  TransferFunction(){ };
  virtual double compute( double k ) = 0;
  virtual ~TransferFunction(){ };
};


enum tf_type{
	total, cdm, baryon, vcdm, vbaryon
};



class TransferFunction_MUSIC : public TransferFunction
{
public:
    std::string m_filename_Tk;
    std::vector<double> m_tab_k, m_tab_Tk_tot, m_tab_Tk_cdm, m_tab_Tk_baryon, m_tab_Tvk_cdm, m_tab_Tvk_baryon;
    gsl_interp_accel *acc_dtot, *acc_dcdm, *acc_dbaryon, *acc_vcdm, *acc_vbaryon;
	gsl_spline *spline_dtot, *spline_dcdm, *spline_dbaryon, *spline_vcdm, *spline_vbaryon;
	
	void read_table( void ){
		std::cerr 
			<< " - reading tabulated transfer function data from file \n"
			<< "    \'" << m_filename_Tk << "\'\n";
			
			std::string line;
			std::ifstream ifs( m_filename_Tk.c_str() );
			
			if(! ifs.good() )
				throw std::runtime_error("Could not find transfer function file \'"+m_filename_Tk+"\'");
			
			m_tab_k.clear();
			m_tab_Tk_tot.clear();
			m_tab_Tk_cdm.clear();
			m_tab_Tk_baryon.clear();
			m_tab_Tvk_cdm.clear();
			m_tab_Tvk_baryon.clear();
			
			while( !ifs.eof() ){
				getline(ifs,line);
				
				if(ifs.eof()) break;
				
				std::stringstream ss(line);
				
				double k, Tkc, Tkb, Tktot, Tkvc, Tkvb;
				ss >> k;
				ss >> Tktot;
				ss >> Tkc;
				ss >> Tkb;
				ss >> Tkvc;
				ss >> Tkvb;
				
				
				m_tab_k.push_back( log10(k) );
				
				m_tab_Tk_tot.push_back( log10(Tktot) );
				m_tab_Tk_baryon.push_back( log10(Tkb) );
				m_tab_Tk_cdm.push_back( log10(Tkc) );
				m_tab_Tvk_cdm.push_back( log10(Tkvc) );
				m_tab_Tvk_baryon.push_back( log10(Tkvb) );
			}
			
			ifs.close();					
	}
	
public:
	TransferFunction_MUSIC( std::string fname )
	: m_filename_Tk( fname )
	{
		read_table( );
		
		acc_dtot = gsl_interp_accel_alloc();
		acc_dcdm = gsl_interp_accel_alloc();
		acc_dbaryon = gsl_interp_accel_alloc();
		acc_vcdm = gsl_interp_accel_alloc();
		acc_vbaryon = gsl_interp_accel_alloc();
		
		spline_dtot = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_dcdm = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_dbaryon = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_vcdm = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_vbaryon = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		
		gsl_spline_init (spline_dtot, &m_tab_k[0], &m_tab_Tk_tot[0], m_tab_k.size() );
		gsl_spline_init (spline_dcdm, &m_tab_k[0], &m_tab_Tk_cdm[0], m_tab_k.size() );
		gsl_spline_init (spline_dbaryon, &m_tab_k[0], &m_tab_Tk_baryon[0], m_tab_k.size() );
		gsl_spline_init (spline_vcdm, &m_tab_k[0], &m_tab_Tvk_cdm[0], m_tab_k.size() );
		gsl_spline_init (spline_vbaryon, &m_tab_k[0], &m_tab_Tvk_baryon[0], m_tab_k.size() );
		
	}
	
	~TransferFunction_MUSIC()
	{
		gsl_spline_free (spline_dtot);
		gsl_spline_free (spline_dcdm);
		gsl_spline_free (spline_dbaryon);
		gsl_spline_free (spline_vcdm);
		gsl_spline_free (spline_vbaryon);
		
		gsl_interp_accel_free (acc_dtot);
		gsl_interp_accel_free (acc_dcdm);
		gsl_interp_accel_free (acc_dbaryon);
		gsl_interp_accel_free (acc_vcdm);
		gsl_interp_accel_free (acc_vbaryon);
	}
	
	inline double extrap_left( double k, const tf_type& type ) 
	{
		if( k<1e-8 )
			return 1.0;
		
		double v1(1.0), v2(1.0);
		switch( type )
		{
			case cdm:
				v1 = m_tab_Tk_cdm[0];
				v2 = m_tab_Tk_cdm[1];
				break;
			case baryon:
				v1 = m_tab_Tk_baryon[0];
				v2 = m_tab_Tk_baryon[1];
				break;
			case vcdm:
				v1 = m_tab_Tvk_cdm[0];
				v2 = m_tab_Tvk_cdm[1];
				break;
			case vbaryon:
				v1 = m_tab_Tvk_baryon[0];
				v2 = m_tab_Tvk_baryon[1];
				break;
			case total: 
				v1 = m_tab_Tk_tot[0];
				v2 = m_tab_Tk_tot[1];
				break;
				
			default:
				throw std::runtime_error("Invalid type requested in transfer function evaluation");
		}
		
		double lk = log10(k);
		double dk = m_tab_k[1]-m_tab_k[0];
		double delk = lk-m_tab_k[0];
		
		//double xi = (v2-v1)/dk;
		return pow(10.0,(v2-v1)/dk*(delk)+v1);
	}
	
	inline double extrap_right( double k, const tf_type& type ) 
	{
		double v1(1.0), v2(1.0);
		
		int n=m_tab_k.size()-1, n1=n-1;
		switch( type )
		{
			case cdm:
				v1 = m_tab_Tk_cdm[n1];
				v2 = m_tab_Tk_cdm[n];
				break;
			case baryon:
				v1 = m_tab_Tk_baryon[n1];
				v2 = m_tab_Tk_baryon[n];
				break;
			case vcdm:
				v1 = m_tab_Tvk_cdm[n1];
				v2 = m_tab_Tvk_cdm[n];
				break;
			case vbaryon:
				v1 = m_tab_Tvk_baryon[n1];
				v2 = m_tab_Tvk_baryon[n];
				break;
			case total: 
				v1 = m_tab_Tk_tot[n1];
				v2 = m_tab_Tk_tot[n];
				break;
				
			default:
				throw std::runtime_error("Invalid type requested in transfer function evaluation");
		}
		
		double lk = log10(k);
		double dk = m_tab_k[n]-m_tab_k[n1];
		double delk = lk-m_tab_k[n];
		
		//double xi = (v2-v1)/dk;
		return pow(10.0,(v2-v1)/dk*(delk)+v2);
	}
	
	inline double compute( double k, tf_type type ){
		
		double lk = log10(k);
		
		//if( lk<m_tab_k[1])
		//	return 1.0;
		
		//if( lk>m_tab_k[m_tab_k.size()-2] );
		//	return m_tab_Tk_cdm[m_tab_k.size()-2]/k/k;
		
		if( k<get_kmin() )
			return extrap_left(k, type );
		
		if( k>get_kmax() )
			return extrap_right(k,type );
		
		
		switch( type )
		{
			case cdm:
				return pow(10.0, gsl_spline_eval (spline_dcdm, lk, acc_dcdm) );
			case baryon:
				return pow(10.0, gsl_spline_eval (spline_dbaryon, lk, acc_dbaryon) );
			case vcdm:
				return pow(10.0, gsl_spline_eval (spline_vcdm, lk, acc_vcdm) );
			case vbaryon:
				return pow(10.0, gsl_spline_eval (spline_vbaryon, lk, acc_vbaryon) );
			case total: 
				return pow(10.0, gsl_spline_eval (spline_dtot, lk, acc_dtot) );
				
			default:
				throw std::runtime_error("Invalid type requested in transfer function evaluation");
		}
		
		return 1.0;
	}
    
    inline double compute( double k )
    {
        return this->compute(k, total);
    }
	
	inline double get_kmin( void ){
		return pow(10.0,m_tab_k[0]);
	}
	
	inline double get_kmax( void ){
		return pow(10.0,m_tab_k[m_tab_k.size()-1]);
	}
};



//! Implementation of abstract base class TransferFunction for the BBKS transfer function 
/*!
    This class implements the analytical fit to the matter transfer
    function by Bardeen, Bond, Kaiser & Szalay (BBKS).
    ( see Bardeen et al. (1986) )
*/
class TransferFunction_BBKS : public TransferFunction{
private:
  double      m_Gamma;

public:
  //! Constructor
  /*!
    \param aCosm Structure of type Cosmology carrying the cosmological parameters
  */
  TransferFunction_BBKS( bool bSugiyama=true )
  {  
    /*m_Gamma = g_Omega_m*g_H0*0.01
      *exp(-g_Omega_b-sqrt(2.0*g_H0*0.01)
	   * g_Omega_b/(g_Omega_m+g_Omega_b+g_Omega_L));*/

    double Omega0 = g_Omega_m;//-g_Omega_b;
    //m_Gamma = Omega0*0.01*g_H0;
    m_Gamma = 0.25;
    
    if( bSugiyama )
      m_Gamma *= exp(-g_Omega_b*(1.0+sqrt(2.0*0.01*g_H0)/Omega0));
  }

  //! computes the value of the BBKS transfer function for mode k (in h/Mpc)
  virtual inline double compute( double k ){
    double q, f1, f2;
    
    q = k/(m_Gamma);
    f1 = log(1.0 + 2.34*q)/(2.34*q);
    f2 = 1.0 + q*(3.89 + q*(259.21 + q*(162.771336 + q*2027.16958081)));
    
    return f1/sqrt(sqrt(f2));

  }
};



//! Implementation of abstract base class TransferFunction for the Eisenstein & Hu transfer function 
/*!
    This class implements the analytical fit to the matter transfer
    function by Eisenstein & Hu (1999).
*/
class TransferFunction_Eisenstein : public TransferFunction{
private:
  double  m_h0;
  double  omhh,		/* Omega_matter*h^2 */
    obhh,		/* Omega_baryon*h^2 */
    theta_cmb,	/* Tcmb in units of 2.7 K */
    z_equality,	/* Redshift of matter-radiation equality, really 1+z */
    k_equality,	/* Scale of equality, in Mpc^-1 */
    z_drag,		/* Redshift of drag epoch */
    R_drag,		/* Photon-baryon ratio at drag epoch */
    R_equality,	/* Photon-baryon ratio at equality epoch */
    sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
    k_silk,		/* Silk damping scale, in Mpc^-1 */
    alpha_c,	/* CDM suppression */
    beta_c,		/* CDM log shift */
    alpha_b,	/* Baryon suppression */
    beta_b,		/* Baryon envelope shift */
    beta_node,	/* Sound horizon shift */
    k_peak,		/* Fit to wavenumber of first peak, in Mpc^-1 */
    sound_horizon_fit,	/* Fit to sound horizon, in Mpc */
    alpha_gamma;	/* Gamma suppression in approximate TF */

  //! private member function: sets internal quantities for Eisenstein & Hu fitting
  void TFset_parameters(float omega0hh, double f_baryon, double Tcmb)
  /* Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
  /* Input: omega0hh -- The density of CDM and baryons, in units of critical dens,
     multiplied by the square of the Hubble constant, in units
     of 100 km/s/Mpc */
  /* 	  f_baryon -- The fraction of baryons to CDM */
  /*        Tcmb -- The temperature of the CMB in Kelvin.  Tcmb<=0 forces use
	    of the COBE value of  2.728 K. */
  /* Output: Nothing, but set many global variables used in TFfit_onek(). 
     You can access them yourself, if you want. */
  /* Note: Units are always Mpc, never h^-1 Mpc. */
  {
    double z_drag_b1, z_drag_b2;
    double alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;
    
    if (f_baryon<=0.0 || omega0hh<=0.0) {
      fprintf(stderr, "TFset_parameters(): Illegal input.\n");
      exit(1);
    }
    omhh = omega0hh;
    obhh = omhh*f_baryon;
    if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
    theta_cmb = Tcmb/2.7;
    
    z_equality = 2.50e4*omhh/POW4(theta_cmb);  /* Really 1+z */
    k_equality = 0.0746*omhh/SQR(theta_cmb);

    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag = 1291*pow(omhh,0.251)/(1+0.659*pow(omhh,0.828))*
		(1+z_drag_b1*pow(obhh,z_drag_b2));
    
    R_drag = 31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));
    R_equality = 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);

    sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
	    log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));

    k_silk = 1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1+pow(10.4*omhh,-0.95));

    alpha_c_a1 = pow(46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
    alpha_c_a2 = pow(12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
    alpha_c = pow(alpha_c_a1,-f_baryon)*
		pow(alpha_c_a2,-CUBE(f_baryon));
    
    beta_c_b1 = 0.944/(1+pow(458*omhh,-0.708));
    beta_c_b2 = pow(0.395*omhh, -0.0266);
    beta_c = 1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));

    y = z_equality/(1+z_drag);
    alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    alpha_b = 2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;

    beta_node = 8.41*pow(omhh, 0.435);
    beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);

    k_peak = 2.5*3.14159*(1+0.217*omhh)/sound_horizon;
    sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*pow(obhh,0.75));

    alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*
		SQR(f_baryon);
    
    return;
  }

  //! private member function: computes transfer function for mode k (k in Mpc)
  inline double TFfit_onek(double k, double *tf_baryon, double *tf_cdm)
  /* Input: k -- Wavenumber at which to calculate transfer function, in Mpc^-1.
   *tf_baryon, *tf_cdm -- Input value not used; replaced on output if
   the input was not NULL. */
  /* Output: Returns the value of the full transfer function fitting formula.
     This is the form given in Section 3 of Eisenstein & Hu (1997).
     *tf_baryon -- The baryonic contribution to the full fit.
     *tf_cdm -- The CDM contribution to the full fit. */
  /* Notes: Units are Mpc, not h^-1 Mpc. */
  {
    double T_c_ln_beta, T_c_ln_nobeta, T_c_C_alpha, T_c_C_noalpha;
    double q, xx, xx_tilde;//, q_eff;
    double T_c_f, T_c, s_tilde, T_b_T0, T_b, f_baryon, T_full;
    //double T_0_L0, T_0_C0, T_0, gamma_eff; 
    //double T_nowiggles_L0, T_nowiggles_C0, T_nowiggles;

    k = fabs(k);	/* Just define negative k as positive */
    if (k==0.0) {
	if (tf_baryon!=NULL) *tf_baryon = 1.0;
	if (tf_cdm!=NULL) *tf_cdm = 1.0;
	return 1.0;
    }

    q = k/13.41/k_equality;
    xx = k*sound_horizon;

    T_c_ln_beta = log(2.718282+1.8*beta_c*q);
    T_c_ln_nobeta = log(2.718282+1.8*q);
    T_c_C_alpha = 14.2/alpha_c + 386.0/(1+69.9*pow(q,1.08));
    T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));

    T_c_f = 1.0/(1.0+POW4(xx/5.4));
    T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*SQR(q)) +
	    (1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*SQR(q));
    
    s_tilde = sound_horizon*pow(1+CUBE(beta_node/xx),-1./3.);
    xx_tilde = k*s_tilde;

    T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*SQR(q));
    T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1+SQR(xx/5.2))+
		alpha_b/(1+CUBE(beta_b/xx))*exp(-pow(k/k_silk,1.4)));
    
    f_baryon = obhh/omhh;
    T_full = f_baryon*T_b + (1-f_baryon)*T_c;

    /* Now to store these transfer functions */
    if (tf_baryon!=NULL) *tf_baryon = T_b;
    if (tf_cdm!=NULL) *tf_cdm = T_c;
    return T_full;
}

public:
  //! Constructor for Eisenstein & Hu fitting for transfer function
  /*!
    \param aCosm structure of type Cosmology carrying the cosmological parameters
    \param Tcmb mean temperature of the CMB fluctuations (defaults to
    Tcmb = 2.726 if not specified)
  */
  TransferFunction_Eisenstein( double Tcmb = 2.726 )
    :  m_h0( g_H0*0.01 ){
    TFset_parameters( (g_Omega_m)*g_H0*g_H0*(0.01*0.01), 
		      g_Omega_b/(g_Omega_m-g_Omega_b),//-aCosm.Omega_b), 
		      Tcmb);
}

  //! Computes the transfer function for k in Mpc/h by calling TFfit_onek
  virtual inline double compute( double k ){
    double tfb, tfcdm, fb, fc;
    TFfit_onek( k*m_h0, &tfb, &tfcdm );
    
    fb = g_Omega_b/(g_Omega_m);
    fc = (g_Omega_m-g_Omega_b)/(g_Omega_m) ;
    
    return fb*tfb+fc*tfcdm;
  }
  

};



#endif
