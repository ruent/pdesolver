#pragma once
#include<vector>
#include<functional>
//#include "lgmarkov/model.h"
#include<cassert>

//start with the heat-like PDE: u_t = s u_xx
// U_t = s U_xx is the convention
//Assume time step size a, space step size b
//Then r below is just  (s*a)/(b*b)

//data below is the whole grid
//by providing a start, part of the grid is treated as an initial value

//FIX GRID NOTATION: x_{t,k}: t is time, k is space coordinate
 
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
    //first and last are boundary and given (first is i =0, last is i = n-1)   
    //start and end are for the initial vector 
    for(size_t i = start+1; i < end; ++i)
    {
        //i+n is the ith space point on the next time iterate, 
        //i.e. data[i+n] is the grid point x_{initial time +1, i}
        data[i+n] = r * data[i-1] + (1-2*r) *data[i] + r*data[i+1];
    }
    //WARNING: boundary conditions are not implemented above because
    //It is ASSUMED THAT boundaries are implemented before this is called
    //so, the boundA, boundB are misleading. I think I first put them there for
    //some interface like appearance over different iteration schemes
}

void implicit_step(vector<double>& data, 
        size_t start, size_t n, double r, double boundA, double boundB)
{
    //J or n, are supposed to be a fixed time t size of partitioning of space dimension
    //should produce x_{i+1} in
    // A x_{i+1} = x_i + boundary vector
    //If space has J nodes, A is (J-2) x(J-2)
    // x are column vectors of size (J-2)
    // boundary vector is a column vec of size J
    //with all but last and first entries are zero
    //(first is boundA, last is boundB)

    //A is used to get all non-baoundary space point values
    // v_(prev,0) and v_(prev,n-1) are boundary points
    // - r v_0 + (1+ 2r) v_1 - r v_2 = v_(prev,1)
    // - r v_1 + (1+ 2r) v_2 - r v_3 = v_(prev,2)
    //...
    //- r v_n-3 + (1+ 2r) v_n-2 - r v_n-1 = v_(prev,n-2)
    //IN The above equations v_0 (i.e. v_{next,0}), and v_n-1 (i.e. v_{next, n-1})
    // are added to the solution later
    // (they are excluded from the matrix A)
    //see the line prevV[0] += ...
    //SO THE FOLLOWING IS THE eqn A* nextV = prevV
    // (1+ 2r) v_1 - r v_2           = v_(prev,1) + r v_0
    // - r v_1 + (1+ 2r) v_2 - r v_3 = v_(prev,2)
    //...
    //- r v_n-3 + (1+ 2r) v_n-2      = v_(prev,n-2) +r v_n-1
    matrix<double> A(n-2,n-2); 
    for (size_t i=0; i< n-2; ++i){
        A[i][i] = 1+2*r;
        if (i< n-3){
            A[i][i+1] = -r;
            A[i+1][i] = -r;
        }        
    }
    //solve: A* nextV = y (y = prevV+ boundary terms)
    //boundary terms vector is of size J-2 with last entry
    //prevV[n-3]. The first entry should be v_{next,0}, and the last one should be 
    //NOTE: boundary terms are from the next time point,
    // the remaining entries of prevV, from 1 through n-4 are from initial time slice
    std::vector<double> prevV(data.begin()+start+1, data.begin()+start+n-1);
    prevV[0] += r*boundA;
    prevV[n-3] += r*boundB; 
     
    //WARNING: I trust the compiler here, where I assume the rest of the entries of prevV
    //are zero. I should fix this at a good time. FIXED ABOVE!
    std::vector<double> nextV = solveEqnTridiagonal(A, prevV); //A * nextV = prevV
    //first and last are boundary and the given   
    //start and end are for the initial vector 
    for(size_t i = 1; i < n-1; ++i){
        //std::cout << "nextV for i= " << i << " is: "<< nextV[i-1] << "\n";
        data[i+start+n] = nextV[i-1];
    }
    //fill the boundary terms as well
    data[start+n] += boundA;
    data[n-1+start+n] += boundB; 
}

//a generalized pointer (Gottschling, 1st edition, p.205)
stepperfunctype explicitstepper = explicit_step; 

//the parameter r could be templatized to allow for
//simple constant values but probabaly it is better
// a non-template version where parameter is always a vector.
//Let it be time-dependent only, so a vector over the time axis.

class discretepde
{    
    //see my notes
protected:    
    size_t I; //time: 0,1,   ..., I-1 (number of time steps)
    size_t J; //space: 0,1,...,J-1
    float Idelta = 0.02; //time increment = 0.02 approx 1 week
    //assuming V_t = 0.5 V_xx is the eqn //space increment > sqrt(Idelta*2*r)
    //Idelta and Jdelta give a stable configuration
    //see my notes... (assume r is always at most 1)
    float Jdelta = 0.01; //space increment
    float Izero = 0.0; //initial time point
    float Jzero; //initial space point
    vector<double> grid; //initial vector and one iterate
    double r; //r is passed to the stepper, look there for the meaning
    stepperfunctype stepper;

public:
    //ctor
    discretepde(size_t _I, size_t _J, double _r,
        stepperfunctype _stepper, float _Jzero, double _Jdelta): 
          I(_I), J(_J), r(_r),
          stepper(_stepper), Jzero(_Jzero), grid(_I*_J), Jdelta(_Jdelta) {}
    virtual double upperboundary(double i) =0;
    virtual double lowerboundary(double i) =0;
    virtual double initialvalue(double x) =0;
    virtual void do_step(vector<double>& data, 
        size_t start, size_t n, double r) = 0;

    //accessor
    double* getgrid(const size_t time)
    {
        if (time * J >= I * J ){
            std::cout << "not a good idea: memory overflow will happen!" << "\n";
            return &grid[0];
        }        
        return &grid[time*J];
    } 

    //accessor
    double getgrid_new(const size_t i)
    {
        if (i >= I * J ){
            std::cout << "not a good idea: memory overflow will happen!" << "\n";
            return grid[0];
        }        
        return grid[i];
    }     

    double getgrid_timespace(const size_t time, const size_t space)
    {
        return getgrid(time)[space];
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
    double initialvalue(double x) override {
        PYBIND11_OVERLOAD_PURE(
            double, //return type
            discretepde, //parent class
            initialvalue, //function name in C++
            x // function argument
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
    variancemap& alpha;
    lgmswapprice& exercise_value;
    std::vector<double> continuation_value;
    unsigned int Idelta_days;
    unsigned int period_days;     
    lgmpde(size_t _I, size_t _J, stepperfunctype _stepper, 
              float _Jzero, variancemap& _alpha, lgmswapprice& _exercise_value,
              std::vector<double> _continuation_value, unsigned int _Idelta_days,
              unsigned int _period_days, double _Jdelta): 
        // the 0.0 in the initilizer is really a design problem, it is required in discrete pde
        // can be remmoved ones lgm_pde is verified to work      
        discretepde(_I, _J, 0.0, _stepper, _Jzero, _Jdelta), alpha(_alpha), 
                  exercise_value(_exercise_value), continuation_value(_continuation_value),
                    Idelta_days(_Idelta_days), period_days(_period_days) 
    {
        //1 week = 7/365, time step size, so T= 7/365 * I
        //space = 0.01, so total space interval is 0.01* J
        fillInitialValues(0);
        double r_;        
        //Idelta = (double)Idelta_days/(double)period_days; 
        Idelta = (double)Idelta_days/365.25; //above one is stupid!       
        //loop hould go over i = 0, upto I-2, because
        // only
        size_t pDate = exercise_value.pDate;
        DayCountCalculator dayc = exercise_value.dayc; //Do you have this class copyable? 
        unsigned int expiry =  exercise_value.s.expiry;
        for (size_t i=0; i<I-1; ++i){
            //fill in the boundary values
            double time_lookup = static_cast<double>(i)*Idelta; 
            //casting is for to mute any concern about C size_t being a trouble, 
            //but I don't think this is really necessary
            //year frac looking is wrong: you need excel dates
            //time_lookup = dayc.yFrac(pDate, expiry - i * Idelta_days);    
            time_lookup = expiry - i * Idelta_days;          
            grid[i*J] = lowerboundary(time_lookup);
            grid[i*J+J-1] = upperboundary(time_lookup);
            //initial value (or final)
            double vol = alpha(time_lookup);
            r_ = 0.5*vol *vol; 
            // A criteria for correct Idelta would be for (I-1)*Idelta
            //should be exercise_date - previous exercise date
            r = r_ * Idelta/(Jdelta * Jdelta); 
            //details in the variance map function 
            //minus sign is because of time change of variables
            //value satisfies the pde: V_t + r_ V_xx = 0 with a final boundary.
            //If t is changed to T-t, PDE becomes V_t - r_ V_xx = 0,
            // i.e. V_t = r_ V_xx , and confirms to our convention.
            //So we can pass this r_ based r (without sign change) into the solver   
            //std::cout <<time_lookup <<" "<< i << " "<< J <<" "<<   r << "\n";     
            stepper(grid, i*J, J, r, lowerboundary(time_lookup), upperboundary(time_lookup));
        }
    }
    //you should bring the LGM swaption prices here 
    double upperboundary(double i){
        double val = grid[J-1];
        return val;
    };
    //lower boundary assumes sufficient thickness and small 3month dicount difference
    //under such assumptions lower boundary should be the continuation value
    //for my current calibration, continuation value is always given by exercising at one future 
    //exercise date
    double lowerboundary(double i){
        //return exercise_value(Jzero);
        double val = grid[0];
        //std::cout << "lower boundary: "<< val <<"\n";
        return grid[0];
        };
    double initialvalue(double x){
        return exercise_value.rebased_swap_price(x);
    };
    void fillInitialValues(size_t start_time_index){
        for (size_t j = 0; j<J; ++j){
            double value = initialvalue(Jzero + j * Jdelta);
            if(continuation_value[j] > value) value = continuation_value[j];
            //if (j%100 == 0)std::cout << "initial values: "<< value << " (index: "<< start_time_index+j<< ")\n";
            grid[start_time_index+j] = value;
        }
    };
    void do_step(vector<double>& data, //YOU NEED TO HAVE THIS FOR INTERFACE PURPOSES
        size_t start, size_t n, double r){
            ;
        }
};

 

class heatpde: public discretepde
{
    public:
    heatpde(size_t _I, size_t _J, double _r, 
        stepperfunctype _stepper, float _Jzero): 
        discretepde(_I,_J, _r, _stepper, _Jzero, 0.01)
    {
        //1 week = 7/365, time step size, so T= 7/365 * I
        //space = 0.01, so total space interval is 0.01* J
        initialvalue(0.0); //fills in the initial value, 0.0 is just a place holder
        //r = 0.5 * Jdelta*Jdelta/Idelta;
        for (size_t i=0; i<I-1; ++i)
        {
            //fill in the boundary values
            grid[i*J] = lowerboundary(i *Idelta);
            grid[i*J+J-1] = upperboundary(i*Idelta);
            //initial value (or final)            
            stepper(grid, i*J, J, r, lowerboundary(i *Idelta), upperboundary(i*Idelta));
            //do_step(grid, i*J, J, r_time);
        }
    }
    //you should bring the LGM swaption prices here 
    double upperboundary(double i){return 0.0;};
    double lowerboundary(double i){return 0.0;};
    double initialvalue(double x) {
    //double x is an after addition, not necessary here, but necessary for lgm_pde
        for (size_t j=0; j<J; ++j)
        {
            double aux = j*Jdelta+Jzero;
            grid[j] = (aux >0.000001 || aux<-0.000001)?  0: 1;
        }   
        return 0.0;      
    }
    void do_step(vector<double>& data, 
        size_t start, size_t n, double r){
             return explicitstepper(data, start, n, r,
                   lowerboundary(0 *Idelta), upperboundary(0*Idelta));
        }
};

 


//you can do all this by passing the pde to the stpper function!
//this might be a good exercise, and
//maybe a better implementation

//there are two competing implementations already present:
// do_step calls explicit stpper
//or I use the stpper that I used in the initialization