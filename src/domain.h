#ifndef DOM
#define DOM

#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "matrix.h"

using namespace std;


struct node_data
{
	double V;
	complex<double> psi;
};


class domain
{
	typedef std::complex<double> comp_t;
	typedef dat<complex<double> > dat_t;

public:

	domain(int rows, int cols, double dx, double dy);

	~domain()
	{
		//delete cn_A;
		//delete cn_b;
		//delete cn_x;
	}

	//utility
	void		allocate_nodes();
	void		dump(ostream &str);
	void		dump_separate(string prefix);
	inline int	idx(int i, int j){	return i*cols_+j; }
	void		copy_to_last();
	
	//potential
	void make_barrier(double V0, double x0, double x1, double y0, double y1);
	void make_wall_barrier(double a, double c);

	void gauss_barrier_y(	double y0, double x0, double x1, 
								double a, double b, double c);
	void gauss_barrier_x(	double x0, double y0, double y1, 
								double a, double b, double c);
	
	//integrators
	void timestep_euler(double dt);
	void timestep_cranknicholson(double dt);
	void cn_init(double dt);
	
	//other
	double	psi_sum();
	void	init_psi(double kx, double ky, double sigma, double x0, double y0);
	
	
	//members
	double		dx_;
	double		dy_;
	int			rows_; //diff. by dy
	int			cols_; //diff. by dx
	node_data**	d_; 
	node_data** d_last_;

	cmat< complex<double> >*		cn_b;
	cmat< complex<double> >*		cn_x;
	cmat< complex<double> >*		cn_A;
		

};



















#endif
