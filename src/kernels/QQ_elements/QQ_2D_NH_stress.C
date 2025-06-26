#include "QQ_2D_NH_stress.h"

#include "MooseMesh.h"
#include "Assembly.h"
#include "MooseVariable.h"
#include "libmesh/quadrature.h"
registerMooseObject("MECHApp", QQ_2D_NH_stress);
template <>
InputParameters validParams<QQ_2D_NH_stress>() {
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<Real>("mu", "mu");
    params.addRequiredParam<Real>("lambda", "lambda");
    params.addRequiredCoupledVar("disp_x", "");
    params.addRequiredCoupledVar("disp_y", "");
    //params.addCoupledVar        ("disp_z", "");
    
    return params;
}

QQ_2D_NH_stress::QQ_2D_NH_stress(const InputParameters & params) :
Kernel(params),
//_component(getParam<unsigned>("component")),
_disp_x_var(coupled("disp_x")),
_disp_y_var(coupled("disp_y")),
//_disp_z_var(_mesh.dimension() == 3 ? coupled("disp_z") : 100000),
_disp_x(coupledValue("disp_x")),
_disp_y(coupledValue("disp_y")),
//_disp_z(coupledValue("disp_z")),
_mu(getParam<Real>("mu")),_lambda(getParam<Real>("lambda"))
{
    
    if (_mesh.dimension() == 3)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    
    if (_mesh.dimension() == 2)
        _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0);
    
    simpson_to_tri6= new int [6];
    
    simpson_to_tri6[0]=0;
    simpson_to_tri6[1]=1;
    simpson_to_tri6[2]=2;
    simpson_to_tri6[3]=3;
    simpson_to_tri6[4]=5;
    simpson_to_tri6[5]=4;
    
    _u_x= new Real [6];
    _u_y= new Real [6];
    
    _local_to_global= new int *[4];
    for (int i=0; i<4;++i)
        _local_to_global[i]=new int [3];
    
    _local_to_global[0][0]=0;
    _local_to_global[0][1]=3;
    _local_to_global[0][2]=5;
    
    _local_to_global[1][0]=1;
    _local_to_global[1][1]=4;
    _local_to_global[1][2]=3;
    
    _local_to_global[2][0]=2;
    _local_to_global[2][1]=5;
    _local_to_global[2][2]=4;
    
    _local_to_global[3][0]=0;
    _local_to_global[3][1]=1;
    _local_to_global[3][2]=2;
    
    
    // here we will store the gradient of basis functions: 0,1,2 are the fine ones, 3 the coarse ones
    _qq_grad_test=new RealVectorValue * [4];
    for (int i=0 ; i < 4; ++i)
        _qq_grad_test[i]=new RealVectorValue [3];
    
    // here we allocate mechanical quantities per element
    U = new RealTensorValue [4];
    F = new RealTensorValue [4];
    C = new RealTensorValue [4];
    S = new RealTensorValue [4];
    
    // SQQ is per node coarse
    SQQ = new RealTensorValue [3];
    
    // from here on, whatever concerns the assembly
    _v_x= new Real [6];
    _v_y= new Real [6];
    
    _V = new RealTensorValue [4];
    _E_lin = new RealTensorValue [4];
    
    _E_lin_QQ = new RealTensorValue [3];
    
    if (_var.number()==_disp_x_var)
    {
        _component=0;
    }
    else if (_var.number()==_disp_y_var)
    {
        _component=1;
    }
    else
    {
        _console<<"error\n exiting\n";
        exit(1);
    }
    
    _h_x= new Real [6];
    _h_y= new Real [6];
    
    _H = new RealTensorValue [4];
    _HV = new RealTensorValue [4];
    _E_lin_H = new RealTensorValue [4];
    _Slin = new RealTensorValue [4];
    
    _HV_QQ= new RealTensorValue [3];
    _Slin_QQ=new RealTensorValue [3];
    
    _console<<"KernelQQ: constructor end";
}

void QQ_2D_NH_stress::computeResidual()
{
    //    std::cout<<"chiamato\n";
    if (_mesh.dimension()==2)
        computeResidual2D();
}

void QQ_2D_NH_stress::computeResidual2D()
{
    /*std::cout<<"numero di punti "<<_qrule->n_points()<<std::endl;
     for (_qp = 0; _qp < _qrule->n_points(); _qp++)
     std::cout<<_JxW[_qp]<<std::endl;
     exit(1);*/
    if (_qrule->n_points()!=6)
    {
        _console<<"You are using a wrong quadrature rule,\n Exiting...";
        exit(1);
    }
    
    for (_qp=0; _qp<6; ++_qp)
    {
        _u_x[_qp]=_disp_x[simpson_to_tri6[_qp]];
        _u_y[_qp]=_disp_y[simpson_to_tri6[_qp]];
    }
    
    for (int tri=0; tri<4; ++tri)
        computeGradient(_q_point[simpson_to_tri6[_local_to_global[tri][0]]],
                        _q_point[simpson_to_tri6[_local_to_global[tri][1]]],
                        _q_point[simpson_to_tri6[_local_to_global[tri][2]]],
                        _qq_grad_test[tri]);
    
    
    // Displacement gradient
    
    for (int tri=0; tri<4; ++tri) {
        
        RealVectorValue Temp1,Temp2,Temp3;
        
        Temp1 = _qq_grad_test[tri][0]*_u_x[_local_to_global[tri][0]]+
        _qq_grad_test[tri][1]*_u_x[_local_to_global[tri][1]]+
        _qq_grad_test[tri][2]*_u_x[_local_to_global[tri][2]];
        Temp2 = _qq_grad_test[tri][0]*_u_y[_local_to_global[tri][0]]+
        _qq_grad_test[tri][1]*_u_y[_local_to_global[tri][1]]+
        _qq_grad_test[tri][2]*_u_y[_local_to_global[tri][2]];
        Temp3 = RealVectorValue(0.0,0.0,0.0);
        
        U[tri]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));
    }
    
    for (int i=0; i<4; ++i)
    {
        F[i]=U[i]+_identity;
        C[i]=F[i].transpose()*F[i];
        
        C[i](2,2)=1.0;
        RealTensorValue invC=C[i].inverse();
        C[i](2,2)=0.0;
        for (int kk=0; kk<3; ++kk)
        {
            invC(2,kk)=0.0;
            invC(kk,2)=0.0;
        }
        Real J=F[i](0,0)*F[i](1,1)-F[i](0,1)*F[i](1,0);
        S[i]=_mu*(_identity-invC)+_lambda*std::log(J)*invC;
    }
    
    for (int nodo_coarse=0; nodo_coarse<3; ++nodo_coarse)
    {
        // we project stress
        SQQ[nodo_coarse]=2.0*S[nodo_coarse]-S[3];
    }
    
    // Here we compute the mass matrix in order to assemble
    Real area=0.0;
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        area+=_JxW[_qp];
    
    DenseMatrix<Number> mass_matrix;
    mass_matrix.resize(3,3);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            if (i==j)
                mass_matrix(i,j)=area/6.0;
            else
                mass_matrix(i,j)=area/12.0;
    
    DenseVector<Number> & re = _assembly.residualBlock(_var.number());
    
    if (re.size()!=6)
    {
        _console<<"somethign went wrong\n, exiting...\n";
        exit(1);
    }
    
    _local_re.resize(re.size());
    _local_re.zero();
    
    for (int _i=0; _i<6; ++_i)
    {
        // we just clean the increment vectors to be on the safe side
        for (int iiii=0; iiii<6; ++iiii)
        {
            _v_x[iiii]=0.0;
            _v_y[iiii]=0.0;
        }
        // we set the increment to one ot the write component
        if (_component==0)
            _v_x[_i]=1.0;
        if (_component==1)
            _v_y[_i]=1.0;
        
        for (int tri=0; tri<4; ++tri) {
            
            RealVectorValue Temp1,Temp2,Temp3;
            
            Temp1 = _qq_grad_test[tri][0]*_v_x[_local_to_global[tri][0]]+
            _qq_grad_test[tri][1]*_v_x[_local_to_global[tri][1]]+
            _qq_grad_test[tri][2]*_v_x[_local_to_global[tri][2]];
            Temp2 = _qq_grad_test[tri][0]*_v_y[_local_to_global[tri][0]]+
            _qq_grad_test[tri][1]*_v_y[_local_to_global[tri][1]]+
            _qq_grad_test[tri][2]*_v_y[_local_to_global[tri][2]];
            Temp3 = RealVectorValue(0.0,0.0,0.0);
            
            _V[tri]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));
        }
        
        for (int tri=0; tri<4; ++tri)
        {
            _E_lin[tri] = 0.5*(F[tri].transpose()*_V[tri]+_V[tri].transpose()*F[tri]);
        }
        
        for (int nodo=0; nodo<3; ++nodo)
            _E_lin_QQ[nodo]=2.0*_E_lin[nodo]-_E_lin[3];
        
        DenseVector<Number> stress(3);
        DenseVector<Number> strain_lin(3);
        DenseVector<Number> res(3);
        
        for (int i_local=0; i_local<2; ++i_local)
            for (int j_local=0; j_local<2; ++j_local)
            {
                for (int nodo_coarse = 0; nodo_coarse < 3; ++nodo_coarse)  // nuova quadratura
                {
                    stress(nodo_coarse)=SQQ[nodo_coarse](i_local,j_local);
                    strain_lin(nodo_coarse)=_E_lin_QQ[nodo_coarse](i_local,j_local);
                }
                mass_matrix.vector_mult(res,stress);
                _local_re(_i) += res.dot(strain_lin);
            }
    }
    
    re += _local_re;
}

void QQ_2D_NH_stress::computeGradient(RealVectorValue x0, RealVectorValue x1, RealVectorValue x2, RealVectorValue * Gradient)
{
    RealTensorValue coordinates;
    coordinates(0,0)=x0(0);
    coordinates(0,1)=x0(1);
    coordinates(0,2)=1.0;
    coordinates(1,0)=x1(0);
    coordinates(1,1)=x1(1);
    coordinates(1,2)=1.0;
    coordinates(2,0)=x2(0);
    coordinates(2,1)=x2(1);
    coordinates(2,2)=1.0;
    
    coordinates=coordinates.inverse();
    Gradient[0]=RealVectorValue(coordinates(0,0),coordinates(1,0),0.0);
    Gradient[1]=RealVectorValue(coordinates(0,1),coordinates(1,1),0.0);
    Gradient[2]=RealVectorValue(coordinates(0,2),coordinates(1,2),0.0);
}



void QQ_2D_NH_stress::computeOffDiagJacobian(unsigned int jvar)
{
    //_console<<"Ho chiamato il computeJacobian. Stai usando le FD?";
    if (_mesh.dimension()==2)
    {
        if (jvar==_disp_x_var || jvar==_disp_y_var)
            computeJacobian2D(jvar);
    }
};

void QQ_2D_NH_stress::computeJacobian2D(unsigned int jvar)
{
    int jcomponent=-999;
    if (jvar==_disp_x_var)
    {
        jcomponent=0;
    }
    else
    {
        if (jvar==_disp_y_var)
        {
            jcomponent=1;
        }
        else
        {
            _console<<"you should not be here, exiting...\n";
            exit(1);
        }
    }
    if (jcomponent<0)
    {
        _console<<"component has not been fixed correctly, exiting...\n";
        exit(1);
    }
    
    
    for (_qp=0; _qp<6; ++_qp)
    {
        _u_x[_qp]=_disp_x[simpson_to_tri6[_qp]];
        _u_y[_qp]=_disp_y[simpson_to_tri6[_qp]];
    }
    
    for (int tri=0; tri<4; ++tri)
        computeGradient(_q_point[simpson_to_tri6[_local_to_global[tri][0]]],
                        _q_point[simpson_to_tri6[_local_to_global[tri][1]]],
                        _q_point[simpson_to_tri6[_local_to_global[tri][2]]],
                        _qq_grad_test[tri]);
    
    
    // Displacement gradient
    
    for (int tri=0; tri<4; ++tri) {
        
        RealVectorValue Temp1,Temp2,Temp3;
        
        Temp1 = _qq_grad_test[tri][0]*_u_x[_local_to_global[tri][0]]+
        _qq_grad_test[tri][1]*_u_x[_local_to_global[tri][1]]+
        _qq_grad_test[tri][2]*_u_x[_local_to_global[tri][2]];
        Temp2 = _qq_grad_test[tri][0]*_u_y[_local_to_global[tri][0]]+
        _qq_grad_test[tri][1]*_u_y[_local_to_global[tri][1]]+
        _qq_grad_test[tri][2]*_u_y[_local_to_global[tri][2]];
        Temp3 = RealVectorValue(0.0,0.0,0.0);
        
        U[tri]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));
    }
    
    for (int i=0; i<4; ++i)
    {
        F[i]=U[i]+_identity;
        C[i]=F[i].transpose()*F[i];
        
        C[i](2,2)=1.0;
        RealTensorValue invC=C[i].inverse();
        C[i](2,2)=0.0;
        for (int kk=0; kk<3; ++kk)
        {
            invC(2,kk)=0.0;
            invC(kk,2)=0.0;
        }
        Real J=F[i](0,0)*F[i](1,1)-F[i](0,1)*F[i](1,0);
        S[i]=_mu*(_identity-invC)+_lambda*std::log(J)*invC;
    }
    
    for (int nodo_coarse=0; nodo_coarse<3; ++nodo_coarse)
    {
        // we project stress
        SQQ[nodo_coarse]=2.0*S[nodo_coarse]-S[3];
    }
    
    // Here we compute the mass matrix in order to assemble
    Real area=0.0;
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        area+=_JxW[_qp];
    
    DenseMatrix<Number> mass_matrix;
    mass_matrix.resize(3,3);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            if (i==j)
                mass_matrix(i,j)=area/6.0;
            else
                mass_matrix(i,j)=area/12.0;
    
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
    _local_ke.resize(6,6);
    _local_ke.zero();
    
    for (int _i=0; _i<6; ++_i)
    {
        // we just clean the increment vectors to be on the safe side
        for (int iiii=0; iiii<6; ++iiii)
        {
            _v_x[iiii]=0.0;
            _v_y[iiii]=0.0;
        }
        // we set the increment to one ot the write component
        if (_component==0)
            _v_x[_i]=1.0;
        if (_component==1)
            _v_y[_i]=1.0;
        
        for (int tri=0; tri<4; ++tri) {
            
            RealVectorValue Temp1,Temp2,Temp3;
            
            Temp1 = _qq_grad_test[tri][0]*_v_x[_local_to_global[tri][0]]+
            _qq_grad_test[tri][1]*_v_x[_local_to_global[tri][1]]+
            _qq_grad_test[tri][2]*_v_x[_local_to_global[tri][2]];
            Temp2 = _qq_grad_test[tri][0]*_v_y[_local_to_global[tri][0]]+
            _qq_grad_test[tri][1]*_v_y[_local_to_global[tri][1]]+
            _qq_grad_test[tri][2]*_v_y[_local_to_global[tri][2]];
            Temp3 = RealVectorValue(0.0,0.0,0.0);
            
            _V[tri]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));
        }
        
        for (int tri=0; tri<4; ++tri)
        {
            _E_lin[tri] = 0.5*(F[tri].transpose()*_V[tri]+_V[tri].transpose()*F[tri]);
        }
        
        for (int nodo=0; nodo<3; ++nodo)
            _E_lin_QQ[nodo]=2.0*_E_lin[nodo]-_E_lin[3];
        
        
        for (int _j=0; _j<6; ++_j)
        {
            // we just clean the increment vectors to be on the safe side
            for (int jjjj=0; jjjj<6; ++jjjj)
            {
                _h_x[jjjj]=0.0;
                _h_y[jjjj]=0.0;
            }
            // we set the increment to one ot the write component
            if (jcomponent==0)
                _h_x[_j]=1.0;
            if (jcomponent==1)
                _h_y[_j]=1.0;
            
            for (int tri=0; tri<4; ++tri) {
                
                RealVectorValue Temp1,Temp2,Temp3;
                
                Temp1 = _qq_grad_test[tri][0]*_h_x[_local_to_global[tri][0]]+
                _qq_grad_test[tri][1]*_h_x[_local_to_global[tri][1]]+
                _qq_grad_test[tri][2]*_h_x[_local_to_global[tri][2]];
                Temp2 = _qq_grad_test[tri][0]*_h_y[_local_to_global[tri][0]]+
                _qq_grad_test[tri][1]*_h_y[_local_to_global[tri][1]]+
                _qq_grad_test[tri][2]*_h_y[_local_to_global[tri][2]];
                Temp3 = RealVectorValue(0.0,0.0,0.0);
                
                _H[tri]=RealTensorValue(Temp1(0),Temp1(1),Temp1(2),Temp2(0),Temp2(1),Temp2(2),Temp3(0),Temp3(1),Temp3(2));
            }
            
            for (int tri=0; tri<4; ++tri)
            {
                _E_lin_H[tri] = (F[tri].transpose()*_H[tri]+_H[tri].transpose()*F[tri]); // I removed the 0.5*
                C[tri](2,2)=1.0;
                F[tri](2,2)=1.0;
                RealTensorValue invC=C[tri].inverse();
                RealTensorValue invFtr=F[tri].inverse().transpose();
                C[tri](2,2)=0.0;
                F[tri](2,2)=0.0;
                for (int kk=0; kk<3; ++kk)
                {
                    invC(2,kk)=0.0;
                    invC(kk,2)=0.0;
                    invFtr(2,kk)=0.0;
                    invFtr(kk,2)=0.0;
                }
                Real J=F[tri](0,0)*F[tri](1,1)-F[tri](0,1)*F[tri](1,0);
                _Slin[tri]= _mu*invC*_E_lin_H[tri]*invC -_lambda*std::log(J)*invC*_E_lin_H[tri]*invC+_lambda*invFtr.contract(_H[tri])*invC;
                
                _HV[tri]=0.5*(_V[tri].transpose()*_H[tri]+_H[tri].transpose()*_V[tri]);

            }
            
            for (int nodo=0; nodo<3; ++nodo)
            {
                _Slin_QQ[nodo]=2.0*_Slin[nodo]-_Slin[3];
                _HV_QQ[nodo]=2.0*_HV[nodo]-_HV[3];
            }
            
            
            
            DenseVector<Number> stress(3);
            DenseVector<Number> strain_lin(3);
            DenseVector<Number> stress_lin(3);
            DenseVector<Number> strain_lin_lin(3);
            DenseVector<Number> res(3);
            DenseVector<Number> res2(3);
            
            for (int i_local=0; i_local<2; ++i_local)
                for (int j_local=0; j_local<2; ++j_local)
                {
                    for (int nodo_coarse = 0; nodo_coarse < 3; ++nodo_coarse)  // nuova quadratura
                    {
                        stress(nodo_coarse)=SQQ[nodo_coarse](i_local,j_local);
                        stress_lin(nodo_coarse)=_Slin_QQ[nodo_coarse](i_local,j_local);
                        strain_lin(nodo_coarse)=_E_lin_QQ[nodo_coarse](i_local,j_local);
                        strain_lin_lin(nodo_coarse)=_HV_QQ[nodo_coarse](i_local,j_local);
                    }
                    mass_matrix.vector_mult(res,stress);
                    mass_matrix.vector_mult(res2,stress_lin);
                    _local_ke(_i,_j) += (res.dot(strain_lin_lin)+res2.dot(strain_lin));
                }
        }
    }
    ke += _local_ke;
    
    
}





















