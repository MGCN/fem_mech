#include "Eikonal.h"
registerMooseObject("MECHApp", Eikonal);

template <>
InputParameters
validParams<Eikonal>()
{
  InputParameters params = Kernel::validParams();
    
  //params.addRequiredParam<Real>("c0","c0");

  params.addRequiredParam<Real>("tau","tau");

  params.addRequiredParam<Real>("c0","c0");

  return params;

}

Eikonal::Eikonal(const InputParameters &parameters) :
Kernel(parameters),
//_c0(getParam<Real>("c0")),
_tau(getParam<Real>("tau")),
_c0(getParam<Real>("c0")),
_conductivity(getMaterialProperty<RealTensorValue>("_K"))
{}

Real
Eikonal::computeQpResidual(){
	
    Real c0 = _c0;

    //Real sqrsigma = sqrt(1400.0);

    _sigma = _conductivity[_qp];

    Real a = _grad_u[_qp]*(_sigma*_grad_test[_i][_qp]);
	
    Real b = c0*std::sqrt( _grad_u[_qp]*(_sigma*_grad_u[_qp]))*_test[_i][_qp];

    Real c = - _tau*_test[_i][_qp];

    //Real rv=a+(theta/sqrsigma)*b+c;
    
    Real rv=a+b+c;

    if (isnan(rv))
	
    {

      std::cout<<"nan in residual  a       "<<a<< std::endl;
	
      std::cout<<"nan in residual  b       "<<b<< std::endl;
	
      std::cout<<"nan in residual  c       "<<c<< std::endl;
	
      std::cout<<"nan in residual  c0   "<<c0<< std::endl;
	
      //std::cout<<"nan in residual  sqrsigma"<<sqrsigma<< std::endl;
	
      std::cout<<_sigma<<std::endl;
	
      exit(1);
	
    }

    if (isinf(rv))
	
    {
	
      std::cout<<"inf in residual  a       "<<a<< std::endl;
	
      std::cout<<"inf in residual  b       "<<b<< std::endl;
	
      std::cout<<"inf in residual  c       "<<c<< std::endl;
	
      std::cout<<"inf in residual  c0   "<<c0<< std::endl;
	
      //std::cout<<"inf in residual  sqrsigma"<<sqrsigma<< std::endl;
	
      std::cout<<_sigma<<std::endl;
	
      exit(1);
	
    }
	
    return rv;
	
}
	
Real
Eikonal::computeQpJacobian(){

    Real c0 = _c0;

    //Real sqrsigma = sqrt(1400.0);

    _sigma =_conductivity[_qp];

    Real a = _grad_phi[_j][_qp]*(_sigma*_grad_test[_i][_qp]);
    
    Real b = c0/std::sqrt( _grad_u[_qp]*(_sigma*_grad_u[_qp]))*(_grad_phi[_j][_qp]*(_sigma*_grad_u[_qp]))*_test[_i][_qp];

    //Real rv=a+(theta/sqrsigma)*b;
    
    Real rv=a+b;
    
    
    if (isnan(rv))
    {
	
            std::cout<<"nan in jacobian"<<std::endl;

            exit(1);

    }

    if (isinf(rv))
    {

      std::cout<<"inf in jacobian"<<std::endl;
	
      exit(1);
        
    }

    return rv;

}

