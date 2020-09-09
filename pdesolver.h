#pragma once
#include<vector>
#include<functional>


//start with the heat-like PDE: u_t = s u_xx
//Assume time step size a, space step size b
//Then  r below is just  (s*a)/(b*b)
 
using stepperfunctype = std::function<void(vector<double>&, 
                 size_t, size_t, double, double, double)>;

//explicit version has boundary terms just to
//be consistent with the implict_step's signature
void explicit_step(vector<double>& data, 
        size_t start, size_t n, double r, 
        double boundA = 0.0, double boundB = 0.0)
{
    //take an initial vector of size n 
    //from the given vector starting at the given index
    //assume start_index+n and start_index+2n-1 are known
    //boundary conditions
    //r is a parameter for the underlying PDE

    size_t end = start+n-1;
    //first and last are boundary and given   
    //start and end are for the initial vector 
    for(size_t i = start+1; i < end; ++i)
    {
        data[i+n] = r * data[i-1] + (1-2*r) *data[i] + r*data[i+1];
    }
}

void implicit_step(vector<double>& data, 
        size_t start, size_t n, double r, double boundA, double boundB)
{
    //should produce x_{i+1} in
    // A_imp x_{i+1} = x_i + boundary vector
    //If space has J nodes, A_imp is (J-2) x(J-2)
    // x are column vectors of size (J-2)
    // boundary vector is a column vec of size J
    //with all but last and first entries are zero
    //(first is boundA, last is boundB)

    matrix<double> A(n-2,n-2); 
    for (size_t i=0; i< n-2; ++i)
    {
        A[i][i] = 1+2*r;
        if (i< n-3) 
        {
            A[i][i+1] = -r;
            A[i+1][i] = -r;
        }        
    }
    //solve: A* nextV = y(prevV+ boundary terms)
    std::vector<double> prevV(data.begin()+start+1, data.begin()+start+n-1);
    prevV[0] += boundA;
    prevV[n-2] += boundB; 
    std::vector<double> nextV = solveEqnTridiagonal(A, prevV);

    size_t end = start+n-1;
    //first and last are boundary and the given   
    //start and end are for the initial vector 
    for(size_t i = 1; i < n-1; ++i)
    {
        data[i+start+n] = nextV[i-1];
    }
}

//a generalized pointer (Gottschling, 1st edition, p.205)
stepperfunctype explicitstepper = explicit_step; 

//the parameter r could be templatized to allow for
//simple constant values but probabaly it is better
// a non-template version where parameter is always a vector.
//Let it be time-dependent only, so a vector over the time axis.

class discretepde
{    
    //see my notes to the code
protected:    
    size_t I; //time: 0,1,   ..., I (number of time steps)
    size_t J; //space: 0,1,...,J
    float Idelta = 0.02; //time increment = 0.02 approx 1 week
    //assuming V_t = 0.5 V_xx is the eqn //space increment > sqrt(Idelta*2*r)
    //Idelta and Jdelta give a stable configuration
    //see my notes... (assume r is always at most 1)
    float Jdelta = 0.1; //space increment
    float Izero = 0.0; //initial time point
    float Jzero; //initial space point
    vector<double> grid; //initial vector and one iterate
    double r; //r is passed to the stepper, look there for the meaning
    stepperfunctype stepper;

public:
    //ctor
    discretepde(size_t _I, size_t _J, double _r,
        stepperfunctype _stepper, float _Jzero): 
          I(_I), J(_J), r(_r),
          stepper(_stepper), Jzero(_Jzero),grid(_I*_J){}
    virtual double upperboundary(double i) =0;
    virtual double lowerboundary(double i) =0;
    virtual void initialvalue() =0;
    virtual void do_step(vector<double>& data, 
        size_t start, size_t n, double r) = 0;

    //accessor
    double* getgrid(const size_t time)
    {
        return &grid[time*J];
    }    
};


class pydiscretepde: public discretepde{
    public:
        //inherit constructors
        using discretepde::discretepde;
    double upperboundary(double i) override {
        PYBIND11_OVERLOAD_PURE(
            double, //return type
            discretepde, //parent class
            upperboundary, //function name in C++
            i   //function argument
        );
    }   
    double lowerboundary(double i) override {
        PYBIND11_OVERLOAD_PURE(
            double, //return type
            discretepde, //parent class
            lowerboundary, //function name in C++
            i   //function argument
        );
    }  
    void initialvalue() override {
        PYBIND11_OVERLOAD_PURE(
            void, //return type
            discretepde, //parent class
            initialvalue //function name in C++
        );
    }   
    void do_step(std::vector<double>& data, 
        size_t start, size_t n, double r) override {
            PYBIND11_OVERLOAD_PURE(
                void,
                discretepde,
                do_step,
                data, start, n, r
            );
        }
 

};

//can templatize over products whose boundary values
//can be calculated
class lgmpde: public discretepde
{
    public:
    lgmpde(size_t _I, size_t _J, double _r, 
        stepperfunctype _stepper, float _Jzero, variancemap alphaSqr): 
        discretepde(_I,_J, _r, _stepper, _Jzero)
    {
        //1 week = 7/365, time step size, so T= 7/365 * I
        //space = 0.01, so total space interval is 0.01* J
        initialvalue();
        double r_time;
        for (size_t i=0; i<I; ++i)
        {
            //fill in the boundary values
            grid[i*J] = lowerboundary(i * Idelta);
            grid[i*J+J-1] = upperboundary(i* Idelta);
            //initial value (or final)
            r_time = 0.5*alphaSqr(i*Idelta); 
            //i-th time step while time increment is Idelta
            //see variancemap definition for alphaSqr
            do_step(grid, i*J, J, r_time);
        }
    }
    //you should bring the LGM swaption prices here 
    double upperboundary(double i){return i;};
    double lowerboundary(double i){return 0.0;};
    void initialvalue(){};
    void do_step(vector<double>& data, 
        size_t start, size_t n, double r) 
        {
            return explicitstepper(data, start, n, r,
                   lowerboundary(0 *Idelta), upperboundary(0*Idelta));
        };
};

class heatpde: public discretepde
{
    public:
    heatpde(size_t _I, size_t _J, double _r, 
        stepperfunctype _stepper, float _Jzero): 
        discretepde(_I,_J, _r, _stepper, _Jzero)
    {
        //1 week = 7/365, time step size, so T= 7/365 * I
        //space = 0.01, so total space interval is 0.01* J
        initialvalue();
        //r = 0.5 * Jdelta*Jdelta/Idelta;
        for (size_t i=0; i<I-1; ++i)
        {
            //fill in the boundary values
            grid[i*J] = lowerboundary(i *Idelta);
            grid[i*J+J-1] = upperboundary(i*Idelta);
            //initial value (or final)
            
            stepper(grid, i*J, J, r, lowerboundary(i *Idelta), upperboundary(i*Idelta));
        }
    }
    //you should bring the LGM swaption prices here 
    double upperboundary(double i){return 0.0;};
    double lowerboundary(double i){return 0.0;};
    void initialvalue()
    {
        for (size_t j=0; j<J; ++j)
        {
            double aux = j*Jdelta+Jzero;
            grid[j] = (aux >0.000001 || aux<-0.000001)?  0: 1;
        }         
    }
    void do_step(vector<double>& data, 
        size_t start, size_t n, double r){}
};


//you can do all this by passing the pde to the stpper function!
//this might be a good exercise, and
//maybe a better implementation

//there are two competing implementations already present:
// do_step calls explicit stpper
//or I use the stpper that I used in the initialization

