#include <iostream>
#include <stdexcept>

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
		Vector(int length_): data(new T[length_]),length(length_){
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
        	std::cout <<  "Copy operator" << std::endl;
        	return *this;
		}

		// Move operator
		Vector& operator=(Vector&& v2){
        	std::cout << "move operator" << std::endl;
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
    		std::cout << "Plus operator" << std::endl;
        	if(length != v2.length){
        		throw std::invalid_argument( "Vectors have different length" );
        	}
        	Vector<std::decay_t<decltype((*this).data[0] + v2.data[0])>> output(length);
        	for(int i = 0 ; i < length ; i++){
            	output.data[i] = data[i] + v2.data[i];
        	}
        	return output;
    	}

    	template<typename U>
    	auto operator-(const Vector<U>& v2) -> Vector<std::decay_t<decltype((*this).data[0] + v2.data[0])>>{
    		std::cout << "Minus operator" << std::endl;
        	if(length != v2.length)
        		throw std::invalid_argument( "Vectors have different length" );
        	Vector<std::decay_t<decltype((*this).data[0] + v2.data[0])>> output(length);
        	for(int i = 0 ; i < length ; i++){
            	output.data[i] = data[i] - v2.data[i];
        	}
        	return output;
    	}

    	template<typename U>
    	auto operator*(const U scalar) -> Vector<std::decay_t<decltype((*this).data[0]*scalar)>>{
    		std::cout << "* operator" << std::endl;
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

template <typename T, typename U>
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



int main(){

	Vector<int> v = Vector<int>({1,2,3,4});



	Vector<int> a(5);
	Vector<int> b = { 1, 2, 3, 4 };
	Vector<int> c = { 1, 2, 3, 4 ,5};

	Vector<int>d(b);
	a=b;
	a.printData();

	Vector<int> v5;
	v5 = std::move(c);
	v5.printData();
	c.printData();

	Vector<int> a1 = { 1, 2, 3, 4 };
	// Vector<double> a2 = { 1.1, 2.1, 3.1, 4.1 };
	Vector<double> a2 = { 1.1, 2.1, 3.1, 4.1 };
	auto v6 = a1+a2;
	v6.printData();
	auto v7 = a1-a2;
	v7.printData();
	auto v8 = a1*2.1;
	v8.printData();
	auto v9 = 2.1*a1;
	v9.printData();

	Vector<int> v10 = {1,2,3};
	Vector<int> v11 = {1,2,3};
	std::cout << "Dot" << std::endl;
	std::cout << dot(v10,v11) << std::endl;
	




	return 0;
}