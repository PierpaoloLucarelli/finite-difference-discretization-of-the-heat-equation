#include <iostream>
#include<cmath>
#include <stdexcept>
#include <map> 
#include <array>
#include <type_traits>
#include <math.h>

template <typename T>
class Vector
{
	public:

		int length;
    	T* data;

		// Default
		Vector(): length(0),data(nullptr){
		}

		// With length
		Vector(int length_): length(length_),data(new T[length_]()){
		}	

		// With init list
		Vector(std::initializer_list<T> list): Vector((int)list.size())
    	{
        	std::uninitialized_copy(list.begin(), list.end(), data);
    	}

    	// Copy constructtor
    	Vector(const Vector& v2)
        : length(v2.length), data(new T[v2.length]) {
        	if(data){
            	for(int i = 0 ; i < length ; i++)
                	data[i] = v2.data[i];
        	}
    	}

    	// Copy operator
	    Vector& operator=(const Vector& v2){
        	if(this != &v2){
            	delete[] data;
            	length = v2.length;
            	data = new T[v2.length];
            	for (auto i=0; i<v2.length; i++){
                	data[i] = v2.data[i];
            	}
        	}
        	return *this;
		}

		// Move operator
		Vector& operator=(Vector&& v2){
        	if(this != &v2){
            	delete[] data;
            	length = v2.length;
            	data = v2.data;
            	v2.length = 0;
            	v2.data = nullptr;
        	}
        	return *this;
    	}
    	template<typename U>
    	auto operator+(const Vector<U>& v2) -> Vector<std::decay_t<decltype((*this).data[0] + v2.data[0])>>{
        	if(length != v2.length){
        		throw std::invalid_argument( "Vectors have different length" );
        	}
        	Vector<std::decay_t<decltype((*this).data[0] + v2.data[0])>> output(length);
        	for(int i = 0 ; i < length ; i++){
            	output.data[i] = data[i] + v2.data[i];
        	}
        	return output;
    	}

    	double norm(){
    		double sum = 0;
    		for(int i = 0 ; i < length ; i++){
            	sum = sum + pow(data[i],2);
        	}
        	return sqrt(sum);
    	}

    	template<typename U>
    	auto operator-(const Vector<U>& v2) -> Vector<std::decay_t<decltype((*this).data[0] + v2.data[0])>>{
        	if(length != v2.length)
        		throw std::invalid_argument( "Vectors have different length" );
        	Vector<std::decay_t<decltype((*this).data[0] + v2.data[0])>> output(length);
        	for(int i = 0 ; i < length ; i++){
            	output.data[i] = data[i] - v2.data[i];
        	}
        	return output;
    	}

    	template<typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value, U>::type>
    	auto operator*(const U scalar) -> Vector<std::decay_t<decltype((*this).data[0]*scalar)>>{
    		Vector<std::decay_t<decltype((*this).data[0]*scalar)>>output(length);
    		for(int i = 0 ; i < length ; i++){
            	output.data[i] = data[i] * scalar;
        	}
        	return output;	
		}

    	// Destructor
    	~Vector()
    	{
        	length=0;
        	delete[] data;
    	}

    	void printData(){
        	for(int i = 0 ; i < length ; i++){
            	std::cout << data[i] << " " << std::endl;
        	}
    	}
};

template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type,typename U>
auto operator*(T scalar, Vector<U> & v) {
	auto temp = v*scalar;
    return temp;
}

template<typename T>
T dot(const Vector<T>& lhs, const Vector<T>& rhs)
{	
	T result = 0;
	if(lhs.length != rhs.length)
        	throw std::invalid_argument( "Vectors have different length" );
    for(int i = 0 ; i < lhs.length ; i++){
    	result += lhs.data[i]*rhs.data[i];
    }
    return result;
}


template<typename T>
class Matrix{
	public:
		const int rows;
		const int columns;
		// using Pair = std::array<int, 2>; 
		std::map<std::array<int, 2>,T> data; 

		Matrix():rows(0),columns(0), data(nullptr){}

		Matrix(int rows_, int columns_): rows(rows_),columns(columns_){
			
		}

		// not sure if this is implemented correctly
		T& operator[] (std::array<int, 2> index){
			auto it = data.find(index);
			return data[index];
		}

		void printData(){
			for (auto const& x : data){
    			std::cout << "[" << x.first[0] << "," << x.first[1] << "]: "
              	<< x.second << std::endl ;
            }
		}
};

template<typename T>
Vector<T> operator*(const Matrix<T>& lhs, const Vector<T>& rhs)
{
   	Vector<T> output(lhs.rows);
	for (auto const& x : lhs.data)
	{
		std::array<int, 2> key = x.first;
		T value = x.second;
		output.data[key[0]] += value*rhs.data[key[1]];
	}
    return output;
}

template<typename T>
int cg(const Matrix<T> &A, Vector<T> &b, Vector<T> &x, T tol, int maxiter)
{
   Vector<T> r = b - (A*x);
   Vector<T> p = r;
   for(int k = 0 ; k < maxiter ; k++){
   		Vector<T> A_p = A*p;
   		auto rr = dot(r,r);
		auto alpha = rr/dot((A_p), p);
		x = x + alpha*p;
		r = r-alpha*A_p;

		if (dot(r, r) < tol*tol){
       		return k;
		}
       	auto beta  = dot(r, r) / rr;
       	p = r + beta*p;
   }
   return -1;
}


template <int n, typename T>
class Heat
{
	public:
		double alpha;
		int m;
		double dt;
		double dx;
		Matrix<T> M;
		Heat(double alpha_, int m_, double dt_): alpha(alpha_), m(m_), dt(dt_),M(pow(m,n),pow(m,n)){
			dx = double(1)/(m+1);	
			for(int j = 0 ; j < pow(m,n) ; j++){
				// iterate each row of M matrix
				for(int k = 0 ; k < n ; k++){
					// each dimention
					int right = j+pow(m,k);
					int left = j-pow(m,k);
					if(right >= 0 && right < pow(m,n))// 0 < x<m^n
						M[{j,right}]=0-(dt/pow(dx,2)*alpha);
					if(left >= 0 && left < pow(m,n))// 0 < x<m^n
						M[{j,left}]=0-(dt/pow(dx,2)*alpha);
				}
				M[{j,j}]=1+2*n*(dt/pow(dx,2))*alpha;
			}
		}

		Vector<T> exact(T t) const{
			Vector<T> result(pow(m,n));
			for(int j = 0 ; j < pow(m,n) ; j++){
				Vector<int> mapping = mapIndexToVector(j);
				Vector<T> x = mapping*dx;
				result.data[j] = getU(x,t);
			}
			return result;
		}

		void inspectMatrix(){
			M.printData();
		}


		Vector<T> solve(T t) const{
			// Vector<T> wl = exact(0);
			Vector<T> wl(pow(m,n));
			for(int j = 0 ; j < pow(m,n) ; j++){
				Vector<int> mapping = mapIndexToVector(j);
				Vector<T> x = mapping*dx;
				wl.data[j] = getU(x,t);
			}
			Vector<T> wl1(pow(m,n));
			for(double t_ = dt ; t_ <= t ; t_+=dt){
				cg<double>(M, wl, wl1, 0.1, 1000);
				wl = wl1;
			}
			return wl;
		}

		Vector<int> mapIndexToVector(int num) const{
			Vector<int> result(n);
			if(num < pow(m,n)){
				for(int i = n-1 ; i >=0 ; i--){
					result.data[i] = (floor(double(num)/pow(m,i)))+1;
					num = num-(pow(m,i)*floor(double(num)/pow(m,i)));
				}
				return result;
			} else{
				std::cout << "Out of bounds" << std::endl;
				return result;
			}
		}

		double getU(Vector<T>x, T t)const{
			if(vectorIsOnBoundary(x)){
				return 0;
			} else if(t==0){
				double result = 1;
				for(int i = 0 ; i < x.length ; i++)
					result = result*sin(M_PI*x.data[i]);
				return result;
			} else{
				return exp(-n*pow(M_PI,2)*alpha*t)*getU(x,0);
			}
		}

		bool vectorIsOnBoundary(Vector<T> v) const{
			for(int i = 0 ; i < v.length ; i++){
				if(v.data[i] == 0 || v.data[i]==1){
					return true;
				}
			}
			return false;
		}


};


int vector_and_matrix_tests(){

	std::cout << "Testing vectors and matrices" << std::endl;

	Matrix<double> M(3,3);
	M[{0,0}]=1;
	M[{0,1}]=2;
	M[{0,2}]=3;
	M[{1,0}]=4;
	M[{1,1}]=5;
	M[{1,2}]=6;
	M[{2,0}]=7;
	M[{2,1}]=8;
	M[{2,2}]=9;

	Vector<double> v1 = {1,2,3};
	Vector<double> v2 = {4,5,6};

	auto v3 = v1 + v2;
	auto v4 = v1 - v2;
	auto v5 = 0.5 * v2;
	auto v6 = 3 * v2;
	auto s1 = dot(v1, v2);
	auto v7 = M * v2;

	v3.printData();
	v4.printData();
	v5.printData();
	v6.printData();
	v7.printData();

	std::cout << "" << std::endl;
	
	return 1;
}

int cg_tests(){

	std::cout << "Testing conjugate gradient" << std::endl;

	Matrix<double> M(2,2);
	M[{0,0}]=4;
	M[{0,1}]=1;
	M[{1,0}]=1;
	M[{1,1}]=3;

	Vector<double> b = { 1,2 };
	Vector<double> x = { 2,1};
	
	cg<double>(M, b, x, 0.02, 1000);
	x.printData();

	std::cout << "" << std::endl;

	return 1;
}

int heat_tests(double alpha, int m, double dt, int t){
	
	std::cout << "Testing heat equation solution" << std::endl;
	Heat<2,double> h(alpha, m, dt);
	Vector<double> e = h.exact(t);
	std::cout << "Exact" << std::endl;
	Heat<2,double> h2(alpha, m, dt);
	Vector<double> s = h2.solve(t);
	e.printData();
	std::cout << "" << std::endl;
	std::cout << "Solve" << std::endl;
	s.printData();
	std::cout << "" << std::endl;

	std::cout << "difference betwen solutions" << std::endl;
	Vector<double> err = e-s;
	std::cout << err.norm() << std::endl;
	return 1;
}



int main(){

	// Vector<int> v = Vector<int>({1,2,3,4});



	// Vector<double> a(5);
	// Vector<double> b = { 1, 2, 3, 4 };
	// Vector<double> c = { 1, 2, 3, 4 ,5};

	// Vector<double>d(b);
	// a=b;
	// a.printData();

	// Vector<double> v5;
	// v5 = std::move(c);
	// v5.printData();
	// c.printData();

	// Vector<int> a1 = { 1, 2, 3, 4 };
	// Vector<double> a2 = { 1.1, 2.1, 3.1, 4.1 };
	// auto v6 = a1+a2;
	// v6.printData();
	// auto v7 = a1-a2;
	// v7.printData();
	// auto v8 = a1*2.1;
	// v8.printData();
	// auto v9 = 2.1*a1;
	// v9.printData();

	// Vector<int> v10 = {1,2,3};
	// Vector<int> v11 = {1,2,3};

	// std::cout << dot(v10,v11) << std::endl;


	// Heat<3,double> h(0.3125, 3, 0.1);
	// Vector<double> e = h.exact(10);
	// Heat<3,double> h2(0.3125, 3, 0.1);
	// Vector<double> s = h2.solve(10);
	// // Vector<double> diff = e-s;
	// e.printData();
	// s.printData();



	vector_and_matrix_tests();
	cg_tests();
	heat_tests(0.3125, 3, 0.1, 1);




	return 0;
}