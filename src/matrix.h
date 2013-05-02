#ifndef cmatX
#define cmatX

#include <vector>
#include <complex>
#include <iostream>
#include <cstring>

using namespace std;

/////////////////////////
//Sparse Data Member
////////////////////////

template<class T>
class dat
{
public:

	dat(T x, int col):
	  x_(x),col_(col){}

	T	x_;
	int	col_;	
};

/////////////////////////
//Sparse Matrix Class
////////////////////////


template<class T>
class cmat
{
	public:
		
		//Sparse matrix data vector and 
		//row index vector
		vector< dat<T>  > 	mat_data_;
		vector< int > 		row_idx_;			
	
		//auxilliary vectors
		T*	diag_table_;
		T*	x_extra;

		cmat();
		~cmat(){/*delete diag_table_*/}
		void print();
		void insert( dat<T> k, int i, int j);
				
		//vector operations
		static complex<double> 	vect_dot( complex<double>* a, complex<double>* b, int n);
		static void vect_add( T* a, T* b, T* res, T alpha, int n);
		static void vect_sub( T* a, T* b, T* res, T alpha, int n);
		static void vect_copy( T* a, T* b, int n);
		
		//Math functions: conjugate gradient, jacobi,
		//and matrix vector product
		void cg(T* x, T* b, int iters);
		void jacobi(T* x, T* b, int iters);
		void mv( T* x, T* b);
		
		void new_row();	//add a new row to the sparse matrix
		inline int m(); //get size of (square) sparse matrix
		
		//make vector of diagonal elements of matrix
		void make_diag_table();
		void get_range(int &start, int &end, int i);
		
		//debug, etc.
		void print_full(ostream &outStr);
		double print_norm(T* a, T* b, int n, int k);
		static void print_dense_vector(T* x, int n, string desc = "");
		
	
};

template<class T>
cmat<T>::cmat():
	x_extra(NULL)
{}
	
template<class T>	
void cmat<T>::print()
{

	int start_idx,end_idx;
	for(unsigned int i = 0 ; i < row_idx_.size() ; i ++ )
	{							
		get_range(start_idx,end_idx,i);

		cout << "row: " << i << endl;
		for(int j = start_idx ; j < end_idx+1; j++)
			cout << "(" << mat_data_[j].col_ <<"-->" << mat_data_[j].x_ << ")" << endl;
					
		cout << endl;
	}
}

template<class T>		
void cmat<T>::insert( dat<T> k, int i, int j)
{
			///////////////
			//insert data k at A(i,j)
			//////////////
}

template<class T>
void cmat<T>::print_dense_vector(T* x, int n, string desc)
{
	cout << "Printing vector: " << desc << endl;
	for(int i =0; i < n; i++)
		cout << x[i] << endl;

	cout << endl;
}

template<class T>
void cmat<T>::make_diag_table()
{
	//it is up to the user to make sure
	//a diag exists
	//for each row	
	diag_table_ = new T[row_idx_.size()];
	int start_idx,end_idx;

	for(int i = 0 ; i < int(row_idx_.size()) ; i ++ )
	{							
		get_range(start_idx,end_idx,i);

		//it would actually be more efficient to 
		//do the diag term and then subtract it later
		for(int j = start_idx ; j < end_idx+1; j++)
		{
			if(mat_data_[j].col_ == i)
				diag_table_[i] = mat_data_[j].x_; 
		}
	}
}


template<class T>
void cmat<T>::get_range(int &start, int &end, int i)
{
	start	= row_idx_[i];
		
	if(i != int(row_idx_.size()) - 1)	end =  row_idx_[i+1] - 1;
	else								end = mat_data_.size() - 1;
}

template<class T>		
complex<double> cmat<T>::vect_dot( complex<double>* a, complex<double>* b, int n)
{			
	complex<double> sum = 0;
	for(int i =0; i < n ;i++)
		sum += a[i]*conj(b[i]);
		
	return sum;		
}

template<class T>
void cmat<T>::vect_add( T* a, T* b, T* res, T alpha, int n)
{
	//a + b
	for(int i = 0; i < n ;i++)
		res[i] = a[i] + alpha*b[i];
}

template<class T>
void cmat<T>::vect_sub( T* a, T* b, T* res, T alpha, int n)
{
	//a - b
	for(int i = 0; i < n ;i++)
		res[i] = a[i] - alpha*b[i];
}		

template<class T>
void cmat<T>::vect_copy( T* a, T* b, int n)
{
	//dest,src,size in bytes
	memcpy(a,b,sizeof(T)*n);			
}

template<class T>
void cmat<T>::cg(T* x, T* b, int iters)
{
	
	double tol = 1e-6;
	
	T* r_cur = 		new T[m()];
	T* r_next = 	new T[m()];
	T* p = 			new T[m()];
	T* temp_Ap = 	new T[m()];
	T alpha,beta,rr_cur,rr_next;
	//last values for x and p not needed
	
	//init

	vect_copy(r_cur,b,m());
	mv(x,temp_Ap);
	vect_sub(r_cur,temp_Ap,r_cur,1.,m());
	vect_copy(p,r_cur,m());
	
	rr_cur = vect_dot(r_cur,r_cur,m());
	
	
	for(int i = 0; i < iters; i++)
	{							
		
		mv(p,temp_Ap); //Ap=A*p
		alpha = rr_cur / vect_dot(p,temp_Ap,m());
	
		vect_add(x,p,x,alpha,m());
		vect_sub(r_cur,temp_Ap,r_next,alpha,m());
		
		rr_next = vect_dot(r_next,r_next,m());
		
		beta = rr_next / rr_cur;
		
		if(iters % 100 == 0)
			//cout << "Res: " << rr_next << endl;
		if( sqrt(abs(rr_next)) < tol) 
			break;
	
		vect_add(r_next,p,p,beta,m());
	
		//next mv and update
		swap(r_cur,r_next);
		rr_cur = rr_next;
	} 
	
	delete[] r_cur;
	delete[] r_next;
	delete[] p;
	delete[] temp_Ap;
			
}

template<class T>
void cmat<T>::print_full(ostream &outstr)
{
	cout << setprecision(4);
	int start_idx,end_idx;
	int rows = m();
	
	for(int i = 0 ; i < rows; i ++ )
	{							
	
		get_range(start_idx,end_idx,i);

		for(int j = 0; j < rows; j++)
		{
			bool found  = false;
			for(int k = start_idx ; k < end_idx+1; k++)
			{
			
				if(mat_data_[k].col_ == j)
				{
					outstr << setw(15) << mat_data_[k].x_;
					found = true;
					break;
				}
			}
			
			if(!found)	
				outstr << setw(15) << T(0);			
		}
		
		outstr << endl;
	}
}

template<class T>
void cmat<T>::jacobi(T* x, T* b, int iters)
{
	int start_idx,end_idx;
	int rows = m();
	double convergence;
	
	if(!x_extra) 
		x_extra = new T[row_idx_.size()];

	for(int k = 0 ; k < iters ; k++ )
	{				
		for(int i = 0 ; i < rows; i ++ )
		{							
			T sigma = b[i];
			get_range(start_idx,end_idx,i);

			//it would actually be more efficient to 
			//do the diag term and then subtract it later
			for(int j = start_idx ; j < end_idx+1; j++)
			{
				if(mat_data_[j].col_ != i)
					sigma -= mat_data_[j].x_*x[mat_data_[j].col_];
			}
	
			x_extra[i] = sigma / diag_table_[i];
		
		}

		swap(x,x_extra);			
		if( k > 10 && k % 10 == 0)
		{
			if(print_norm(x,x_extra,row_idx_.size(),k) < 1e-12)
			{	
				cout << "Solution converged." << endl;
				break;
			}
		}
		//print_dense_vector(x,row_idx_.size() );
	
	}	
}

template<class T>
void cmat<T>::new_row()
{
	row_idx_.push_back( mat_data_.size() );
}

template<class T>
void cmat<T>::mv( T* x, T* b)
{
	//b is dest
	int start_idx,end_idx;
	
	for(int i=0;i<int(row_idx_.size());i++)
	{
		get_range(start_idx,end_idx,i);
		
		b[i] = 0;	
		for(int j = start_idx; j < end_idx+1; j++)
			b[i] += mat_data_[j].x_ * x[ mat_data_[j].col_ ];       									
	}			
}

template<class T>
double cmat<T>::print_norm(T* a, T* b, int n, int k)
{
	double sum =0;
	for(int i = 0; i < n; i++)
		sum += abs(a[i] - b[i])*abs(a[i] - b[i]);

	sum = sqrt(sum);
	cout << "Convergence at iteration " << k << " " << sum << endl;
	return sum;
}

template<class T>
inline int cmat<T>::m(){ return row_idx_.size(); }

#endif
