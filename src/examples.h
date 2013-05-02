#ifndef UTIL_1
#define	UTIL_1

#include "examples.h"
#include "matrix.h"
#include "domain.h"


namespace examples
{


	void matrix_test()
	{
		typedef cmat<complex<double> > cmatcd;
		typedef dat<complex<double> > datcd;
		/////////////
		//test for cmat class
		/////////////
		cmat<complex<double> > A;

		A.row_idx_.push_back(0);
		A.mat_data_.push_back( datcd(3.2,0) );
		A.mat_data_.push_back( datcd(0.2,2) );
		A.mat_data_.push_back( datcd(0.2,4) );


		A.new_row();
		A.mat_data_.push_back( datcd(0.5,0) );
		A.mat_data_.push_back( datcd(2.2,1) );
		A.mat_data_.push_back( datcd(0.2,3) );

		A.new_row();
		A.mat_data_.push_back( datcd(0.8,1) );
		A.mat_data_.push_back( datcd(1.2,2) );

		A.new_row();
		A.mat_data_.push_back( datcd(0.2,1) );
		A.mat_data_.push_back( datcd(5.2,3) );

		A.new_row();
		A.mat_data_.push_back( datcd(0.6,3) );
		A.mat_data_.push_back( datcd(8.2,4) );
		

		A.print();
		A.print_full(cout);
		complex<double>* b = new complex<double>[A.m()];
		complex<double>* x = new complex<double>[A.m()];

		//x is guess
		for(int i =0; i < A.m(); i++)
		{
			x[i] = complex<double>(1,0);
			b[i] = complex<double>(1,0);
		}

		A.make_diag_table();
		A.cg(x,b,20);

		A.print_dense_vector(x,A.m());
		
		A.mv(x,b);
		A.print_dense_vector(b,A.m());
	}


	void matrix_test2()
	{
		typedef cmat<complex<double> > cmatcd;
		typedef dat<complex<double> > datcd;
		typedef std::complex<double> comp_t;
		typedef dat<complex<double> > dat_t;
		
		int n = 2000;
		/////////////
		//test for cmat class
		/////////////
		cmat<complex<double> > A;
		
		complex<double> a = 2.65;
		complex<double> bb = 0.5;
		complex<double> c = 0.1;
	
		//
		for(int i = 0; i < n ; i++)
		{
			for(int j = 0; j < n ; j++)
			{
				//add a new sparse row
				A.new_row();			
			
				if(i == 0 || j ==0 || i == n-1 || j == n-1)
					A.mat_data_.push_back( dat_t(comp_t(1.,0),i*n + j) );
				//fill in matrix bands
				else
				{			
					A.mat_data_.push_back( dat_t(a,i*n+j) );
					
					A.mat_data_.push_back( dat_t(bb,i*n+j +1) );
					A.mat_data_.push_back( dat_t(bb,i*n+j -1) );
					
					A.mat_data_.push_back( dat_t(c,i*n+j + n) );
					A.mat_data_.push_back( dat_t(c,i*n+j - n) );
				}
			}			
		}
		
		//A.print();
		//A.print_full();
		
		complex<double>* b = new complex<double>[A.m()];
		complex<double>* x = new complex<double>[A.m()];

		//x is guess
		for(int i =0; i < A.m(); i++)
		{
			x[i] = complex<double>(1,0);
			b[i] = complex<double>(1,0);
		}

		A.make_diag_table();
		//A.jacobi(x,b,50);
		A.cg(x,b,50);
		//A.print_dense_vector(x,A.m());
		
		A.mv(x,b);
		//A.print_dense_vector(b,A.m());
	}



	void setup_potential(){}


	void run_cn()
	{
	
		double dx =		.005;
		double dy =		.005;
		double x_max =	1;
		double y_max =	1.5;
		double dt =		1e-5;
		double n =		200;
		ofstream mat_out;
		vector<double> norms;
		double normloss;
		mat_out.open("cn_matrix.txt");

		domain* mydom = new domain(int(x_max/dx),int(y_max/dy),dx,dy);
		mydom->allocate_nodes();
		mydom->init_psi(0,20,.1,x_max/2,y_max/2);
		
		//exp(1x) * 10
		mydom->make_wall_barrier(10.,1);
		mydom->cn_init(dt);
	
		//mydom->cn_A->print_full(mat_out);

		mydom->copy_to_last();
	
		for(int i = 0; i < n; i++)
		{			
			//mydom->dump(cout);
			mydom->timestep_cranknicholson(dt);
			//mydom->psi_sum();				
			if(i % 10 == 0)
			{				
				cout << "Step: " << i << endl;
				norms.push_back( mydom->psi_sum());
				normloss = ( (norms[norms.size()-1] - norms[0]) / norms[0])*100;
				cout << "Norm Loss is: " << normloss << endl;

				mydom->dump_separate( string( "output_step" + to_string( (long long) i) +"_") );
			}
		}
	}

	void run_cn2()
	{
	
		/////////////////////////
		//Free Prop. Experiment
		////////////////////////

		double dx =		.025;
		double dy =		.005;
		double x_max =	1;
		double y_max =	1.5;
		double dt =		12.5e-6;
		double n =		20000;
		vector<double> norms;
		double normloss;
	
		domain* mydom = new domain(int(y_max/dy),int(x_max/dx),dx,dy);
		mydom->allocate_nodes();
		mydom->init_psi(0,20,.1,x_max/2,y_max/3);
		
		//exp(1x) * 10
		mydom->make_wall_barrier(5.,30);		
		mydom->copy_to_last();
		mydom->cn_init(dt);

		for(int i = 0; i < n; i++)
		{			
			mydom->timestep_cranknicholson(dt);
			if(i % 100 == 0)
			{
				cout << "\n++++++++++++\nStep: " << i << endl;
				norms.push_back( mydom->psi_sum());
				normloss = ( (norms[norms.size()-1] - norms[0]) / norms[0])*100;
				cout << "Norm Loss is: " << normloss << "%" << endl;

				mydom->dump_separate( string( "output_step" + to_string( (long long) i) +"_") );
			}
		}
	}

	void run_cn3()
	{
	
		/////////////////////////
		//Potential Barrier/Tunneling Experiment
		////////////////////////

		double dx =		.01;
		double dy =		.01;
		double x_max =	1.5;
		double y_max =	1;
		double dt =		5e-6;
		double n =		20000;
		vector<double> norms;
		double normloss;
		
		domain* mydom = new domain(int(y_max/dy),int(x_max/dx),dx,dy);
		mydom->allocate_nodes();
		mydom->init_psi(20,0,.1,x_max/3,y_max/2);
		
		//exp(1x) * 10
		mydom->make_wall_barrier(5.,30);
		
		//Make potential barrier and init.
		//Crank-Nicholson Matrix, A
		mydom->gauss_barrier_x(x_max/2,0,y_max,1e3,x_max/2,.01);
		mydom->copy_to_last();
		mydom->cn_init(dt);

		//do timestepping
		for(int i = 0; i < n; i++)
		{			
			mydom->timestep_cranknicholson(dt);
		
			if(i % 100 == 0)
			{
				cout << "\n++++++++++++\nStep: " << i << endl;
				norms.push_back( mydom->psi_sum());
				normloss = ( (norms[norms.size()-1] - norms[0]) / norms[0])*100;
				cout << "Norm Loss is: " << normloss << "%" << endl;

				mydom->dump_separate( string( "output_step" + to_string( (long long) i) +"_") );
			}
		}
	}

	void run_cn4()
	{

		/////////////////////////
		//Double Slit Experiment
		////////////////////////
			
		double dx =		.01;
		double dy =		.01;
		double x_max =	2.5;
		double y_max =	4;
		double dt =		10e-6;
		double n =		20000;
		double k =		60;
		double aperature_size = (2*3.14159) / k;
		vector<double> norms;
		double normloss;
		
		domain* mydom = new domain(int(y_max/dy),int(x_max/dx),dx,dy);
		mydom->allocate_nodes();
		mydom->init_psi(k,0,.2,x_max/3,y_max/2);
		
		//exp(1x) * 10
		mydom->make_wall_barrier(5.,30);
		
		//Make the double slit:
		mydom->gauss_barrier_x(x_max/2,0,y_max/2-1.5*aperature_size,1e3,x_max/2,.015);
		mydom->gauss_barrier_x(x_max/2,y_max/2+1.5*aperature_size,y_max,1e3,x_max/2,.015);
		mydom->gauss_barrier_x(x_max/2,y_max/2-0.5*aperature_size,y_max/2+0.5*aperature_size,1e3,x_max/2,.015);
				
		mydom->copy_to_last();
		mydom->cn_init(dt);

		for(int i = 0; i < n; i++)
		{			
			mydom->timestep_cranknicholson(dt);
			
			if(i % 100 == 0)
			{
				cout << "\n++++++++++++\nStep: " << i << endl;
				norms.push_back( mydom->psi_sum());
				normloss = ( (norms[norms.size()-1] - norms[0]) / norms[0])*100;
				cout << "Norm Loss is: " << normloss << "%" << endl;

				mydom->dump_separate( string( "output_step" + to_string( (long long) i) +"_") );
			}
		}
	}
}





#endif
