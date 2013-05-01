#include "domain.h"

domain::domain(int rows, int cols, double dx, double dy):
	  rows_(rows),cols_(cols),dx_(dx),dy_(dy)
	{	
	  
		cn_A = NULL;	//new cmat(rows*cols,rows*cols);
		cn_b = NULL;	//new cmat(rows*cols,1);
		cn_x = NULL;	//new cmat(rows*cols,1);  
	  
	 }

void domain::allocate_nodes()
{
	d_ =		new node_data*[rows_*cols_];
	d_last_ =	new node_data*[rows_*cols_];
		
	for(int i = 0; i < rows_*cols_; i++)
	{
		d_[i] =			new node_data;
		d_[i]->psi =	comp_t(0.,0.);
		d_[i]->V =		0.;		

		d_last_[i] =		new node_data;
		d_last_[i]->psi =	comp_t(0.,0.);
		d_last_[i]->V =		0.;
	}

	printf("Domain allocated.  Mem. required is about: %i bytes",
		2*((rows_*cols_)*(sizeof(void*) + sizeof(node_data))) );
}

void domain::dump(ostream &str)
{
	cout << "Printing Psi: \n";
	for(int i = 0; i < rows_; i++){
		for(int j = 0; j < cols_; j++){
			str << setw(25) << d_[i*cols_+j]->psi;
		}
		str << "\n";
	}

	str << "Printing V: \n";
	for(int i = 0; i < rows_; i++){
		for(int j = 0; j < cols_; j++){
			str << setw(15) << d_[i*cols_+j]->V;
		}
		str << "\n";
	}
}

void domain::dump_separate(string prefix)
{

	string v_name = prefix + "_potential.txt";
	string psi_name_r = prefix + "_psi_real.txt";
	string psi_name_im = prefix + "_psi_imag.txt";
	string psi_mag_name = prefix + "_psi_mag.txt";

	ofstream f;
	f.open(v_name.c_str() );
		
	for(int i = 0; i < rows_; i++){
		for(int j = 0; j < cols_; j++){
			f << setw(20) << d_[i*cols_+j]->V;
		}
		f << "\n";
	}
		
	f.close();

	f.open(psi_name_r.c_str());
	for(int i = 0; i < rows_; i++){
		for(int j = 0; j < cols_; j++){
			f << setw(20) << d_[i*cols_+j]->psi.real();
		}
		f << "\n";
	}

	f.close();
		
	f.open(psi_name_im.c_str());
	for(int i = 0; i < rows_; i++){
		for(int j = 0; j < cols_; j++){
			f << setw(20) << d_[i*cols_+j]->psi.imag();
		}
		f << "\n";
	}

	f.close();


	f.open(psi_mag_name.c_str());
	for(int i = 0; i < rows_; i++){
		for(int j = 0; j < cols_; j++){
			f << setw(20) << abs(d_[i*cols_+j]->psi)*abs(d_[i*cols_+j]->psi);
		}
		f << "\n";
	}
	f.close();

}

void domain::make_barrier(double V0, double x0, double x1, double y0, double y1)
{		
	//set up specified potential barrier
	//(note the addition)
	for(int i = 0; i < rows_; i++){
		for(int j = 0; j < cols_; j++){
		
		if( j*dx_ < x1 && j*dx_ > x0 && i*dy_ < y1 && i*dy_ > y0)
			d_[i*cols_+j]->V += V0;

		}
	}
}

void domain::make_wall_barrier(double a, double c)
{
		
	for(int i = 0; i < rows_; i++){
		for(int j = 0; j < cols_; j++){
		
			d_[i*cols_+j]->V += c * (	exp( a*(i*dy_ - dy_*(rows_-1))) +//				;// +											
										exp( a*(j*dx_ - dx_*(cols_-1))) +//; //+// +
										exp( -1.*a*(j*dx_))  +//  +//;
										exp( -1.*a*(i*dy_)) ); 
		}
	}
}

void domain::copy_to_last()
{
	for(int i = 0; i < rows_*cols_; i++)
	{
		d_last_[i]->V =		d_[i]->V;
		d_last_[i]->psi =	d_[i]->psi;
	}
}
void domain::timestep_euler(double dt)
{
	//alias swap?
	std::swap(d_last_,d_);
	
	comp_t a = -1.*(dt/ (2.*comp_t(0,1.)*dx_*dx_) );
	comp_t b = -1.*(dt/ (2.*comp_t(0,1.)*dy_*dy_) );		
	comp_t c = dt/(2.*comp_t(0,1.));
		
	int ilast, inext, jlast, jnext;

	for(int i = 0; i < rows_; i++)
	{
		ilast = i-1; inext = i+1;
		for(int j = 0; j < cols_; j++)
		{			
			jlast = j-1; jnext = j+1;
				
			if(i == 0 || j ==0 || i == rows_-1 || j == cols_-1)
				d_[idx(i,j)]->psi = 0.;					//cn_A->mat_data_.push_back( dat_t(comp_t(1.,0),i*cols_ + j) );
			//fill in matrix bands
			else
			{		

				d_[idx(i,j)]->psi = b*d_last_[idx(ilast,j)]->psi + 
									b*d_last_[idx(inext,j)]->psi +
									a*d_last_[idx(i,jlast)]->psi + 
									a*d_last_[idx(i,jnext)]->psi +
									2.*(-a + -b + c*d_last_[idx(i,j)]->V + 0.5)*d_last_[idx(i,j)]->psi;	
				
				//d_[idx(i,j)]->psi = d_last_[idx(i,j)]->psi;	

			}
		}
	}	
}

void domain::timestep_cranknicholson(double dt)
{

	//alias swap?
	std::swap(d_last_,d_);
	
	comp_t a = -1.*(dt/ (4.*comp_t(0,1.)*dx_*dx_) );
	comp_t b = -1.*(dt/ (4.*comp_t(0,1.)*dy_*dy_) );		
	comp_t c = dt/(2.*comp_t(0,1.));	
	///////////////////////////
	//set top and bottom of b!!
	////////////////////////////
	
	//malloc x and b
	comp_t* x_vect = new comp_t[rows_*cols_];
	comp_t* b_vect = new comp_t[rows_*cols_];
	
	//create b
	for(int i = 0; i < rows_; i++)
	{
		for(int j = 0; j < cols_; j++)
		{

			if(i == 0 || j ==0 || i == rows_-1 || j == cols_-1)
				b_vect[i*cols_ + j] = d_last_[idx(i,j)]->psi;
			else
			//we are also missing some rows of b here...
			//+ or -?
			b_vect[i*cols_ + j] =			-1.*(	b*d_last_[idx(i-1,j)]->psi + b*d_last_[idx(i+1,j)]->psi +
													a*d_last_[idx(i,j-1)]->psi + a*d_last_[idx(i,j+1)]->psi +
													(-2.*a + -2.*b + c*d_last_[idx(i,j)]->V + 1.)*d_last_[idx(i,j)]->psi	);

		}
	}
	
	//make guess from d_last
	for(int i = 0; i < rows_*cols_; i++)
		x_vect[i] = d_last_[i]->psi;
		
	//jacobi solve
	//cn_A->make_diag_table();
	//cout << endl;
	//cn_A->jacobi(x_vect,b_vect,200);	
	cn_A->cg(x_vect,b_vect,200);
	//copy x to d_
	for(int i = 0; i < rows_*cols_; i++)
		d_[i]->psi = x_vect[i];
		
	//cleanup (we should just make these class members)
	delete[] x_vect;
	delete[] b_vect;
		
}

void domain::cn_init(double dt)
{
	comp_t a = -1.*(dt/ (4.*comp_t(0,1.)*dx_*dx_) );
	comp_t b = -1.*(dt/ (4.*comp_t(0,1.)*dy_*dy_) );		
	comp_t c = dt/(2.*comp_t(0,1.));
	
	comp_t self_term;

	cn_A = new cmat<comp_t>();
	
	for(int i = 0; i < rows_; i++)
	{
		for(int j = 0; j < cols_; j++)
		{
			//add a new sparse row
			cn_A->new_row();
			
			//enforce d_next = d on edges
			if(i == 0 || j ==0 || i == rows_-1 || j == cols_-1)
				cn_A->mat_data_.push_back( dat_t(comp_t(1.,0),i*cols_ + j) );
			//fill in matrix bands
			else
			{
				self_term = (-2.*a + -2.*b + c*d_[i*cols_+j]->V - 1.);

				cn_A->mat_data_.push_back( dat_t(self_term,i*cols_+j) );
				
				cn_A->mat_data_.push_back( dat_t(a,i*cols_+j +1) );
				cn_A->mat_data_.push_back( dat_t(a,i*cols_+j -1) );
				
				cn_A->mat_data_.push_back( dat_t(b,i*cols_+j + cols_) );
				cn_A->mat_data_.push_back( dat_t(b,i*cols_+j - cols_) );

			}			
		}
	}
}

double domain::psi_sum()
{
	double sum = 0;

	for(int i = 0; i < rows_; i++)
	{
		for(int j = 0; j < cols_; j++)
		{
			sum += abs(d_[idx(i,j)]->psi)*abs(d_[idx(i,j)]->psi)*dx_*dy_;
		}
	}

	//cout << "sum is: " << sum << endl;
	return sum;
}

void domain::init_psi(double kx, double ky, double sigma, double x0, double y0)
{

	double x,y,dist2,kr;
	for(int i = 0; i < rows_; i++){
			
		y = i*dy_;				
		for(int j = 0; j < cols_; j++){
				
			x =		j*dx_;
			dist2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);
			kr =	kx*x + ky*y;				
			d_[i*cols_+j]->psi = 5*exp( -1.*dist2 /(2*sigma*sigma) )*exp(comp_t(0.,kr));
			
		}
	}
}

void domain::gauss_barrier_y(	double y0, double x0, double x1, 
								double a, double b, double c)
{
	//cols_ .. x
	//rows_ .. y



	int x_start = int( x0 / dx_ );
	int x_end = ceil((x1 / dx_));

	for(int i = x_start; i < x_end; i++)
	{
		for(int j = 0; j < cols_; j++)
		{
			d_[idx(i,j)]->V += a*exp( -( (j*dy_-b)*(j*dy_-b) / (2.*c*c) ) );
		}
	}
}

void domain::gauss_barrier_x(	double x0, double y0, double y1, 
								double a, double b, double c)
{
	//cols_ .. x
	//rows_ .. y

	int y_start = int( y0 / dy_ );
	int y_end = ceil((y1 / dy_));

	for(int i = y_start; i < y_end; i++) //dy
	{
		for(int j = 0; j < cols_; j++)
		{
			d_[idx(i,j)]->V += a*exp( -( (j*dx_-b)*(j*dx_-b) / (2.*c*c) ) );
		}
	}
}


