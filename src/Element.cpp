#include "Element.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------WEIGHTED INTEGRATION POINTS-------------------------
//------------------------------------------------------------------------------
template <int DIM>
void Element<DIM>::setIntegPointWeightFunction_FEM()
{
    QuadShapFunction shapeQuad;

    // NORMAL QUADRATURE
    NormalQuad nQuad = NormalQuad();
    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {
        intPointWeightFunctionPrev_FEM[index] = intPointWeightFunction_FEM[index];
        intPointWeightFunction_FEM[index] = 0.;
        index++;
    };

    int LNN = 4 * DIM - 2;

    double xsi[DIM], phi_[LNN];
    index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        for (int i = 0; i < DIM; i++) xsi[i] = nQuad.PointListFem(index, i);
        
        shapeQuad.evaluateFem(xsi, phi_);

        for (int i = 0; i < LNN; i++)
        {
            intPointWeightFunction_FEM[index] += (*nodes_)[connect_[i]]->getWeightFunction() * phi_[i];
        };

        index++;
    };

    // SPECIAL QUADRATURE
    SpecialQuad sQuad = SpecialQuad();
    index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {
        intPointWeightFunctionSpecialPrev_FEM[index] = intPointWeightFunctionSpecial_FEM[index];
        intPointWeightFunctionSpecial_FEM[index] = 0.;
        index++;
    };

    index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        for (int i = 0; i < DIM; i++) xsi[i] = sQuad.PointListFem(index, i);
        
        shapeQuad.evaluateFem(xsi, phi_);

        for (int i = 0; i < LNN; i++)
        {
            intPointWeightFunctionSpecial_FEM[index] += (*nodes_)[connect_[i]]->getWeightFunction() * phi_[i];
        };

        index++;
    };
};

template <int DIM>
void Element<DIM>::setIntegPointWeightFunctionFINE_ISO(double &glueZoneThickness, double &arlequinEpsilon)
{

    QuadShapFunction shapeQuad;

    int LNN = 18 * DIM - 27;

    // NORMAL QUADRATURE
    NormalQuad nQuad = NormalQuad();
    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {
        intPointWeightFunctionPrev_ISO[index] = intPointWeightFunction_ISO[index];
        intPointWeightFunction_ISO[index] = 0.;
        index++;
    };

    // data for evaluation of NURBS functions
    double wpc[LNN], phi_[LNN], xsi[DIM];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListIso(index, i);
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        double distFuncIP = 0.;
        for (int i = 0; i < LNN; i++)
        {
            double r = (*nodes_)[connect_[i]]->getDistFunction();
            distFuncIP += r * phi_[i];
        };

        double wFuncValue;
        if (distFuncIP >= glueZoneThickness)
        {
            wFuncValue = 1. - arlequinEpsilon;
        }
        else
        {
            wFuncValue = (1. - arlequinEpsilon) / glueZoneThickness * distFuncIP;
            if (wFuncValue > (1. - arlequinEpsilon))
                wFuncValue = 1. - arlequinEpsilon;
        };

        intPointWeightFunction_ISO[index] = wFuncValue;

        index++;
    };

    // SPECIAL QUADRATURE
    SpecialQuad sQuad = SpecialQuad();
    index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {
        intPointWeightFunctionSpecialPrev_ISO[index] = intPointWeightFunctionSpecial_ISO[index];
        intPointWeightFunctionSpecial_ISO[index] = 0.;
        index++;
    };

    index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {

        for (int i = 0; i < DIM; i++)
            xsi[i] = sQuad.PointListIso(index, i);
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        double distFuncIP = 0.;
        for (int i = 0; i < LNN; i++)
        {
            double r = (*nodes_)[connect_[i]]->getDistFunction();
            distFuncIP += r * phi_[i];
        };

        double wFuncValue;
        if (distFuncIP >= glueZoneThickness)
        {
            wFuncValue = 1. - arlequinEpsilon;
        }
        else
        {
            wFuncValue = (1. - arlequinEpsilon) / glueZoneThickness * distFuncIP;
            if (wFuncValue > (1. - arlequinEpsilon))
                wFuncValue = 1. - arlequinEpsilon;
        };

        intPointWeightFunctionSpecial_ISO[index] = wFuncValue;

        index++;
    };
};

template <int DIM>
void Element<DIM>::setIntegPointWeightFunctionCOARSE_ISO(double &glueZoneThickness, double &arlequinEpsilon)
{

    QuadShapFunction shapeQuad;

    int LNN = 18 * DIM - 27;

    // NORMAL QUADRATURE
    NormalQuad nQuad = NormalQuad();
    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {
        intPointWeightFunctionPrev_ISO[index] = intPointWeightFunction_ISO[index];
        intPointWeightFunction_ISO[index] = 0.;
        index++;
    };

    // data for evaluation of NURBS functions
    double wpc[LNN], phi_[LNN], xsi[DIM];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListIso(index, i);
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        double distFuncIP = 0;
        for (int i = 0; i < LNN; i++)
        {
            double r = (*nodes_)[connect_[i]]->getDistFunction();
            distFuncIP += r * phi_[i];
        };

        double wFuncValue;
        if (distFuncIP < 0)
        {
            wFuncValue = 1.;
        }
        else
        {
            if (distFuncIP >= glueZoneThickness)
            {
                wFuncValue = arlequinEpsilon;
            }
            else
            {
                wFuncValue = 1. - (1. - arlequinEpsilon) / glueZoneThickness * distFuncIP;
                if (wFuncValue < arlequinEpsilon)
                    wFuncValue = arlequinEpsilon;
            };
        };

        intPointWeightFunction_ISO[index] = wFuncValue;

        index++;
    };

    // SPECIAL QUADRATURE
    SpecialQuad sQuad = SpecialQuad();
    index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {
        intPointWeightFunctionSpecialPrev_ISO[index] = intPointWeightFunctionSpecial_ISO[index];
        intPointWeightFunctionSpecial_ISO[index] = 0.;
        index++;
    };

    index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {

        for (int i = 0; i < DIM; i++)
            xsi[i] = sQuad.PointListIso(index, i);
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        double distFuncIP = 0;
        for (int i = 0; i < LNN; i++)
        {
            double r = (*nodes_)[connect_[i]]->getDistFunction();
            distFuncIP += r * phi_[i];
        };

        double wFuncValue;
        if (distFuncIP < 0)
        {
            wFuncValue = 1.;
        }
        else
        {
            if (distFuncIP >= glueZoneThickness)
            {
                wFuncValue = arlequinEpsilon;
            }
            else
            {
                wFuncValue = 1. - (1. - arlequinEpsilon) / glueZoneThickness * distFuncIP;
                if (wFuncValue < arlequinEpsilon)
                    wFuncValue = arlequinEpsilon;
            };
        };

        intPointWeightFunctionSpecial_ISO[index] = wFuncValue;

        index++;
    };
};

//------------------------------------------------------------------------------
//-----------------SPATIAL TRANSFORMATION - JACOBIAN MATRIXES-------------------
//------------------------------------------------------------------------------
template <int DIM>
void Element<DIM>::getJacobianMatrix_FEM(double &djac_, double *xsi, double **Jac, double **ainv_)
{
    int LNN = 4 * DIM - 2;
    QuadShapFunction shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi[i] = new double[LNN];

    shapeQuad.evaluateGradientFem(xsi, dphi);

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            Jac[i][j] = 0.0;

    for (int i = 0; i < LNN; i++)
    {
        double xna[DIM];
        for (int j = 0; j < DIM; j++)
        {
            xna[j] = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(j) +
                     (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(j);
            for (int k = 0; k < DIM; k++) Jac[j][k] += xna[j] * dphi[k][i];
        };
    };

    MatrixDouble ainvAux(DIM, DIM);

    // Defining the Jacobian matrix
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainvAux(i, j) = Jac[i][j];

    // Computing the jacobian determinant
    djac_ = ainvAux.determinant();

    // if (djac_ <=0) std::cout << "Houston we have a problem: Negative jacobian" << std::endl;

    djac_ = abs(djac_);

    // Computing Jacobian inverse (transposed)
    ainvAux = ainvAux.inverse().transpose();
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            ainv_[i][j] = ainvAux(i,j);
        };
    };

    for (int i = 0; i < DIM; ++i) delete[] dphi[i];
    delete[] dphi;

    return;
};

template <int DIM>
void Element<DIM>::getJacobianMatrix_COARSE_FEM(double *xsi, double **ainv_)
{
    int LNN = 4 * DIM - 2;
    QuadShapFunction shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi[i] = new double[LNN];

    shapeQuad.evaluateGradientFem(xsi, dphi);

    double Jac[DIM][DIM] = {};

    for (int i = 0; i < LNN; i++)
    {
        double xna[DIM];
        for (int j = 0; j < DIM; j++)
        {
            xna[j] = alpha_f * (*nodesC_)[connectC_[i]]->getCoordinateValue(j) +
                     (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousCoordinateValue(j);
            for (int k = 0; k < DIM; k++)
                Jac[j][k] += xna[j] * dphi[k][i];
        };
    };

    MatrixDouble ainvAux(DIM, DIM);

    // Defining the Jacobian matrix
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainvAux(i, j) = Jac[i][j];

    // Computing Jacobian inverse (transposed)
    ainvAux = ainvAux.inverse().transpose();
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainv_[i][j] = ainvAux(i, j);

    for (int i = 0; i < DIM; ++i)
        delete[] dphi[i];
    delete[] dphi;

    return;
};

template <int DIM>
void Element<DIM>::getJacobianMatrix_ISO(double &djac_, double *xsi, double **quadJacMat, double **ainv_)
{

    int LNN = 18 * DIM - 27;
    QuadShapFunction shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi[i] = new double[LNN];

    double wpc[LNN];
    for (int i = 0; i < LNN; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparameters), Npatch_);

    double dx_dxsi[DIM][DIM] = {};
    for (int i = 0; i < LNN; i++)
    {
        double xna[DIM];
        for (int j = 0; j < DIM; j++)
        {
            xna[j] = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(j) +
                     (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(j);
            for (int k = 0; k < DIM; k++) dx_dxsi[j][k] += xna[j] * dphi[k][i];
        };
    };

    MatrixDouble ainvAux(DIM, DIM);

    // Defining the Jacobian matrix
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainvAux(i, j) = dx_dxsi[i][j];

    // Computing Jacobian inverse (transposed)
    ainvAux = ainvAux.inverse().transpose();
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            ainv_[i][j] = ainvAux(i, j);
        };
    };

    // Computing the quadrature jacobian matrix - dx/dxqsi
    double dxsi_dqxsi[DIM][DIM] = {};
    for (int i = 0; i < DIM; i++)
    {
        int deg = (*iparameters)[Npatch_]->getDegree(i);
        int npc = (*iparameters)[Npatch_]->getNcp(i);
        int size = deg + npc + 1;
        double knot[size];
        (*iparameters)[Npatch_]->getKnot(i, size, knot);
        double uF = knot[inc_[i]];
        double uS = knot[inc_[i] + 1];
        dxsi_dqxsi[i][i] = 0.5 * (uS - uF);
    };

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            quadJacMat[i][j] = 0.0;
        };
    };

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k  < DIM; k++) quadJacMat[i][j] += dxsi_dqxsi[k][j] * dx_dxsi[i][k];
        };
    };

    MatrixDouble quadJacMatAux(DIM, DIM);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            quadJacMatAux(i, j) = quadJacMat[i][j];

    djac_ = quadJacMatAux.determinant();

    djac_ = abs(djac_);

    for (int i = 0; i < DIM; ++i)
        delete[] dphi[i];
    delete[] dphi;

    return;
};

template <int DIM>
void Element<DIM>::getJacobianMatrix_COARSE_ISO(double *xsi, double **ainv_)
{
    int LNN = 18 * DIM - 27;
    QuadShapFunction shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi[i] = new double[LNN];

    double wpc[LNN];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *inc_ = (*nodesC_)[connectC_[LNN - 1]]->getINC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparametersC), NpatchC_);

    double dx_dxsi[DIM][DIM] = {};
    for (int i = 0; i < LNN; i++)
    {
        double xna[DIM];
        for (int j = 0; j < DIM; j++)
        {
            xna[j] = alpha_f * (*nodesC_)[connectC_[i]]->getCoordinateValue(j) +
                     (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousCoordinateValue(j);
            for (int k = 0; k < DIM; k++) dx_dxsi[j][k] += xna[j] * dphi[k][i];
        };
    };

    MatrixDouble ainvAux(DIM, DIM);

    // Defining the Jacobian matrix
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainvAux(i, j) = dx_dxsi[i][j];

    // Computing Jacobian inverse (transposed)
    ainvAux = ainvAux.inverse().transpose();
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainv_[i][j] = ainvAux(i, j);
  

    for (int i = 0; i < DIM; ++i)
        delete[] dphi[i];
    delete[] dphi;

    return;
};

template <int DIM>
void Element<DIM>::getQuadJacobianMatrix_ISO(double *xsi, double **quadJacMatInv)
{

    int LNN = 18 * DIM - 27;
    QuadShapFunction shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[DIM];
    for (int i = 0; i < DIM; ++i) dphi[i] = new double[LNN];

    double wpc[LNN];
    for (int i = 0; i < LNN; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparameters), Npatch_);

    double dx_dxsi[DIM][DIM] = {};
    for (int i = 0; i < LNN; i++)
    {
        double xna[DIM];
        for (int j = 0; j < DIM; j++)
        {
            xna[j] = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(j) +
                     (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(j);
            for (int k = 0; k < DIM; k++) dx_dxsi[j][k] += xna[j] * dphi[k][i];
        };
    };

    // Computing the quadrature jacobian matrix - dx/dxqsi
    double dxsi_dqxsi[DIM][DIM] = {};
    for (int i = 0; i < DIM; i++)
    {
        int deg = (*iparameters)[Npatch_]->getDegree(i);
        int npc = (*iparameters)[Npatch_]->getNcp(i);
        int size = deg + npc + 1;
        double knot[size];
        (*iparameters)[Npatch_]->getKnot(i, size, knot);
        double uF = knot[inc_[i]];
        double uS = knot[inc_[i] + 1];
        dxsi_dqxsi[i][i] = 0.5 * (uS - uF);
    };

    double quadJacMat[DIM][DIM] = {};
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k  < DIM; k++) quadJacMat[i][j] += dxsi_dqxsi[k][j] * dx_dxsi[i][k];

    MatrixDouble quadJacMatInvAux(DIM, DIM);

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            quadJacMatInvAux(i, j) = quadJacMat[i][j];

    quadJacMatInvAux = quadJacMatInvAux.inverse();

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            quadJacMatInv[i][j] = quadJacMatInvAux(i, j);

    for (int i = 0; i < DIM; ++i)
        delete[] dphi[i];
    delete[] dphi;

    return;
}

//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template <int DIM>
void Element<DIM>::getSpatialDerivatives_FEM(double *xsi, double **ainv_, double **dphi_dx)
{

    QuadShapFunction shapeQuad;

    int LNN = 4 * DIM - 2;

    double **dphi;
    dphi = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi[i] = new double[LNN];

    shapeQuad.evaluateGradientFem(xsi, dphi);

    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < LNN; ++j)
            dphi_dx[i][j] = 0.;

    // Quadratic shape functions spatial first derivatives
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < LNN; k++)
                dphi_dx[i][k] += ainv_[i][j] * dphi[j][k];

    for (int i = 0; i < DIM; ++i)
        delete[] dphi[i];
    delete[] dphi;

    return;
};

template <int DIM>
void Element<DIM>::getSpatialDerivatives_ISO(double *xsi, double **ainv_, double **dphi_dx)
{

    QuadShapFunction shapeQuad;

    int LNN = 18 * DIM - 27;

    double **dphi;
    dphi = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi[i] = new double[LNN];

    double wpc[LNN] = {};
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparameters), Npatch_);

    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < LNN; ++j)
            dphi_dx[i][j] = 0.;

    // Quadratic shape functions spatial first derivatives
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < LNN; k++)
                dphi_dx[i][k] += ainv_[i][j] * dphi[j][k];

    for (int i = 0; i < DIM; ++i)
        delete[] dphi[i];
    delete[] dphi;
    return;
};

template <int DIM>
void Element<DIM>::getSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double **dphi_dx)
{

    QuadShapFunction shapeQuad;

    int LNN = 18 * DIM - 27;

    double **dphi;
    dphi = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi[i] = new double[LNN];

    double wpc[LNN];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *inc_ = (*nodesC_)[connectC_[LNN - 1]]->getINC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparametersC), NpatchC_);

    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < LNN; ++j)
            dphi_dx[i][j] = 0.;

    // Quadratic shape functions spatial first derivatives
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < LNN; k++)
                dphi_dx[i][k] += ainv_[i][j] * dphi[j][k];

    for (int i = 0; i < DIM; ++i)
        delete[] dphi[i];
    delete[] dphi;
    return;
};

template <int DIM>
void Element<DIM>::getSecondSpatialDerivatives_FEM(double **ainv_, double ***ddphi_dx)
{

    QuadShapFunction shapeQuad;

    int LNN = 4 * DIM - 2;

    double ainvT_[DIM][DIM] = {};
    double inter[DIM][DIM] = {};

    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            for (int k = 0; k < LNN; ++k)
                ddphi_dx[i][j][k] = 0.;

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainvT_[i][j] = ainv_[j][i];

    double ***ddphi;
    ddphi = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphi[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphi[i][j] = new double[LNN];
    };

    shapeQuad.evaluateHessianFem(ddphi);

    for (int nf = 0; nf < LNN; nf++)
    {

        // Quadratic shape functions spatial second derivatives
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                for (int k = 0; k < DIM; k++)
                    inter[i][j] += ainv_[i][k] * ddphi[k][j][nf];

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                for (int k = 0; k < DIM; k++)
                    ddphi_dx[i][j][nf] += inter[i][k] * ainvT_[k][j];

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                inter[i][j] = 0.0;
    };

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; j++)
        {
            delete[] ddphi[i][j];
        };
        delete[] ddphi[i];
    };
    delete[] ddphi;

    return;
};

template <int DIM>
void Element<DIM>::getSecondSpatialDerivatives_ISO(double *xsi, double **ainv_, double ***ddphi_dx)
{

    QuadShapFunction shapeQuad;
    int LNN = 18 * DIM - 27;

    double ainvT_[DIM][DIM] = {};
    double inter[DIM][DIM] = {};

    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            for (int k = 0; k < LNN; ++k)
                ddphi_dx[i][j][k] = 0.;

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainvT_[i][j] = ainv_[j][i];

    double ***ddphi;
    ddphi = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphi[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphi[i][j] = new double[LNN];
    };

    double wpc[LNN] = {};
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    shapeQuad.evaluateHessianIso(xsi, ddphi, wpc, inc_, (*iparameters), Npatch_);

    // Quadratic shape functions spatial second derivatives

    for (int nf = 0; nf < LNN; nf++)
    {

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                for (int k = 0; k < DIM; k++)
                    inter[i][j] += ainv_[i][k] * ddphi[k][j][nf];

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                for (int k = 0; k < DIM; k++)
                    for (int nf = 0; nf < LNN; nf++)
                        ddphi_dx[i][j][nf] += inter[i][k] * ainvT_[k][j];

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                inter[i][j] = 0.0;
    };

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; j++)
        {
            delete[] ddphi[i][j];
        };
        delete[] ddphi[i];
    };
    delete[] ddphi;

    return;
};

template <int DIM>
void Element<DIM>::getSecondSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double ***ddphi_dx)
{

    QuadShapFunction shapeQuad;
    int LNN = 18 * DIM - 27;

    double ainvT_[DIM][DIM] = {};
    double inter[DIM][DIM] = {};

    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            for (int k = 0; k < LNN; ++k)
                ddphi_dx[i][j][k] = 0.;

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ainvT_[i][j] = ainv_[j][i];

    double ***ddphi;
    ddphi = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphi[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphi[i][j] = new double[LNN];
    };

    double wpc[LNN] = {};
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *inc_ = (*nodesC_)[connectC_[LNN - 1]]->getINC();

    shapeQuad.evaluateHessianIso(xsi, ddphi, wpc, inc_, (*iparametersC), NpatchC_);

    // Quadratic shape functions spatial second derivatives

    for (int nf = 0; nf < LNN; nf++)
    {

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                for (int k = 0; k < DIM; k++)
                    inter[i][j] += ainv_[i][k] * ddphi[k][j][nf];

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                for (int k = 0; k < DIM; k++)
                    for (int nf = 0; nf < LNN; nf++)
                        ddphi_dx[i][j][nf] += inter[i][k] * ainvT_[k][j];

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                inter[i][j] = 0.0;
    };

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; j++)
        {
            delete[] ddphi[i][j];
        };
        delete[] ddphi[i];
    };
    delete[] ddphi;

    return;
};

//------------------------------------------------------------------------------
//--------------------------INTERPOLATES VARIABLES------------------------------
//------------------------------------------------------------------------------

template <int DIM>
void Element<DIM>::getInterpForce(int &LNN, double *phi_, double *f_)
{

    for (int i = 0; i < DIM; i++)
    {
        f_[i] = 0.;
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f_[j] += (*nodes_)[connect_[i]]->getNodalForce(j) * phi_[i];
        };
    };
}

template <int DIM>
void Element<DIM>::getInterpCoord(int &LNN, double *phi_, double *x_, double *xPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        x_[i] = 0.;
        xPrev_[i] = 0.;
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            x_[j] += (*nodes_)[connect_[i]]->getCoordinateValue(j) * phi_[i];
            xPrev_[j] += (*nodes_)[connect_[i]]->getPreviousCoordinateValue(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpCoord_ISO(int &LNN, double *phi_, double *x_, double *xPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        x_[i] = 0.;
        xPrev_[i] = 0.;
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            x_[j] += (*nodes_)[connect_[i]]->getCoordinateValue(j) * phi_[i];
            xPrev_[j] += (*nodes_)[connect_[i]]->getPreviousCoordinateValue(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpVel(int &LNN, double *phi_, double *u_, double *uPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        u_[i] = 0.;
        uPrev_[i] = 0.;
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            u_[j] += (*nodes_)[connect_[i]]->getVelocity(j) * phi_[i];
            uPrev_[j] += (*nodes_)[connect_[i]]->getPreviousVelocity(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpVelDer(int &LNN, double **dphi_dx, double **du_dx, double **duPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            du_dx[i][j] = 0.;
            duPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                du_dx[j][k] += (*nodes_)[connect_[i]]->getVelocity(j) * dphi_dx[k][i];
                duPrev_dx[j][k] += (*nodes_)[connect_[i]]->getPreviousVelocity(j) * dphi_dx[k][i];
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpSecondVelDer(int &LNN, double ***ddphi_dx, double ***ddu_dxdx, double ***dduPrev_dxdx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                ddu_dxdx[i][j][k] = 0.;
                dduPrev_dxdx[i][j][k] = 0.;
            };
        };
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int p = 0; p < DIM; p++)
                {
                    ddu_dxdx[j][k][p] += (*nodes_)[connect_[i]]->getVelocity(j) * ddphi_dx[k][p][i];
                    dduPrev_dxdx[j][k][p] += (*nodes_)[connect_[i]]->getPreviousVelocity(j) * ddphi_dx[k][p][i];
                };
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpVelCoarse(int &LNN, double *phi_, double *u_, double *uPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        u_[i] = 0.;
        uPrev_[i] = 0.;
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            u_[j] += (*nodesC_)[connectC_[i]]->getVelocity(j) * phi_[i];
            uPrev_[j] += (*nodesC_)[connectC_[i]]->getPreviousVelocity(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpVelDerCoarse(int &LNN, double **dphi_dx, double **du_dx, double **duPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            du_dx[i][j] = 0.;
            duPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                du_dx[j][k] += (*nodesC_)[connectC_[i]]->getVelocity(j) * dphi_dx[k][i];
                duPrev_dx[j][k] += (*nodesC_)[connectC_[i]]->getPreviousVelocity(j) * dphi_dx[k][i];
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpSecondVelDerCoarse(int &LNN, double ***ddphi_dx, double ***ddu_dxdx, double ***dduPrev_dxdx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                ddu_dxdx[i][j][k] = 0.;
                dduPrev_dxdx[i][j][k] = 0.;
            };
        };
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int p = 0; p < DIM; p++)
                {
                    ddu_dxdx[j][k][p] += (*nodesC_)[connectC_[i]]->getVelocity(j) * ddphi_dx[k][p][i];
                    dduPrev_dxdx[j][k][p] += (*nodesC_)[connectC_[i]]->getPreviousVelocity(j) * ddphi_dx[k][p][i];
                };
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpMeshVel(int &LNN, double *phi_, double *uMesh_, double *uMeshPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        uMesh_[i] = 0.;
        uMeshPrev_[i] = 0.;
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            uMesh_[j] += (*nodes_)[connect_[i]]->getMeshVelocity(j) * phi_[i];
            uMeshPrev_[j] += (*nodes_)[connect_[i]]->getPreviousMeshVelocity(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpMeshVelDer(int &LNN, double **dphi_dx, double **duMesh_dx, double **duMeshPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            duMesh_dx[i][j] = 0.;
            duMeshPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                duMesh_dx[j][k] += (*nodes_)[connect_[i]]->getMeshVelocity(j) * dphi_dx[k][i];
                duMeshPrev_dx[j][k] += (*nodes_)[connect_[i]]->getPreviousMeshVelocity(j) * dphi_dx[k][i];
            };
        };
    };

 };

template <int DIM>
void Element<DIM>::getInterpMeshVelCoarse(int &LNN, double *phi_, double *uMesh_, double *uMeshPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        uMesh_[i] = 0.;
        uMeshPrev_[i] = 0.;
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            uMesh_[j] += (*nodesC_)[connectC_[i]]->getMeshVelocity(j) * phi_[i];
            uMeshPrev_[j] += (*nodesC_)[connectC_[i]]->getPreviousMeshVelocity(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpMeshVelDerCoarse(int &LNN, double **dphi_dx, double **duMesh_dx, double **duMeshPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            duMesh_dx[i][j] = 0.;
            duMeshPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                duMesh_dx[j][k] += (*nodesC_)[connectC_[i]]->getMeshVelocity(j) * dphi_dx[k][i];
                duMeshPrev_dx[j][k] += (*nodesC_)[connectC_[i]]->getPreviousMeshVelocity(j) * dphi_dx[k][i];
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpAccel(int &LNN, double *phi_, double *accel_, double *accelPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        accel_[i] = 0.;
        accelPrev_[i] = 0.;
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            accel_[j] += (*nodes_)[connect_[i]]->getAcceleration(j) * phi_[i];
            accelPrev_[j] += (*nodes_)[connect_[i]]->getPreviousAcceleration(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpAccelDer(int &LNN, double **dphi_dx, double **daccel_dx, double **daccelPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            daccel_dx[i][j] = 0.;
            daccelPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                daccel_dx[j][k] += (*nodes_)[connect_[i]]->getAcceleration(j) * dphi_dx[k][i];
                daccelPrev_dx[j][k] += (*nodes_)[connect_[i]]->getPreviousAcceleration(j) * dphi_dx[k][i];
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpAccelDerCoarse(int &LNN, double **dphi_dx, double **daccel_dx, double **daccelPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            daccel_dx[i][j] = 0.;
            daccelPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                daccel_dx[j][k] += (*nodesC_)[connectC_[i]]->getAcceleration(j) * dphi_dx[k][i];
                daccelPrev_dx[j][k] += (*nodesC_)[connectC_[i]]->getPreviousAcceleration(j) * dphi_dx[k][i];
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpPress(int &LNN, double *phi_, double &press_)
{
    press_ = 0.;

    for (int i = 0; i < LNN; i++)
        press_ += (*nodes_)[connect_[i]]->getPressure() * phi_[i];
};

template <int DIM>
void Element<DIM>::getInterpPressDer(int &LNN, double **dphi_dx, double *dpress_dx)
{

    for (int i = 0; i < DIM; i++)
        dpress_dx[i] = 0.;

    for (int i = 0; i < LNN; i++)
        for (int j = 0; j < DIM; j++)
            dpress_dx[j] += (*nodes_)[connect_[i]]->getPressure() * dphi_dx[j][i];
};

template <int DIM>
void Element<DIM>::getInterpSecondPressDer(int &LNN, double ***ddphi_dx, double **ddpress_dxdx)
{

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ddpress_dxdx[i][j] = 0.;

    for (int i = 0; i < LNN; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                ddpress_dxdx[j][k] += (*nodes_)[connect_[i]]->getPressure() * ddphi_dx[j][k][i];
};

template <int DIM>
void Element<DIM>::getInterpSecondPressDerCoarse(int &LNN, double ***ddphi_dx, double **ddpress_dxdx)
{

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ddpress_dxdx[i][j] = 0.;

    for (int i = 0; i < LNN; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                ddpress_dxdx[j][k] += (*nodesC_)[connectC_[i]]->getPressure() * ddphi_dx[j][k][i];
};

template <int DIM>
void Element<DIM>::getInterpLambda(int &LNN, double *phi_, double *lambda_)
{

    for (int i = 0; i < DIM; i++)
        lambda_[i] = 0.;

    for (int i = 0; i < LNN; i++)
        for (int j = 0; j < DIM; j++)
            lambda_[j] += (*nodes_)[connect_[i]]->getLagrangeMultiplier(j) * phi_[i];
};

template <int DIM>
void Element<DIM>::getInterpLambdaDer(int &LNN, double **dphi_dx, double **dlambda_dx)
{

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            dlambda_dx[i][j] = 0.;

    for (int i = 0; i < LNN; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                dlambda_dx[j][k] += (*nodes_)[connect_[i]]->getLagrangeMultiplier(j) * dphi_dx[k][i];
};

//------------------------------------------------------------------------------
//--------------------COMPUTE BEZIER TRANSFORMATION MATRIX----------------------
//------------------------------------------------------------------------------
template <int DIM>
void Element<DIM>::getDirMatrixC(int &deg, int &dim_, int &ind, double *knot_, double **matrixC_)
{

    double matrixCN_[deg + 1][deg + 1] = {};
    double alphas_[deg + 1] = {};

    int a_ = deg;
    int b_ = a_ + 1;
    int nb = 0;               // number of Bézier Elements
    int i_, mult, r, s, save; // auxiliar
    double numer, alpha;

    for (int i = 0; i <= deg; i++)
    {
        matrixC_[i][i] = 1.;
        matrixCN_[i][i] = 1.;
    };
    while (b_ <= ind + 1)
    {
        for (int i = 0; i <= deg; i++)
        {
            for (int j = 0; j <= deg; j++)
            {
                matrixC_[i][j] = matrixCN_[i][j];
            };
        };
        for (int i = 0; i <= deg; i++)
        {
            matrixCN_[i][i] = 1.;
        };
        i_ = b_;
        while ((b_ < dim_ - 1) && (knot_[b_ + 1] == knot_[b_]))
            b_++;
        mult = b_ - i_ + 1;
        if (mult < deg)
        {
            numer = knot_[b_] - knot_[a_];
            for (int j = deg; j > mult; j--)
            {
                alphas_[j - mult - 1] = numer / (knot_[a_ + j] - knot_[a_]);
            };
            r = deg - mult; // Insert knot r times
            for (int j = 1; j <= r; j++)
            {
                save = r - j;
                s = mult + j;
                for (int k = deg; k >= s; k--)
                {
                    alpha = alphas_[k - s];
                    for (int i = 0; i <= deg; i++)
                    {
                        matrixC_[i][k] = alpha * matrixC_[i][k] + (1.0 - alpha) * matrixC_[i][k - 1];
                    };
                };
                int cc = deg - j;
                for (int i = save; i <= j + save; i++)
                {
                    matrixCN_[i][save] = matrixC_[cc][deg];
                    cc++;
                };
            };
            nb = nb + 1;
            if (b_ <= ind + 1)
            {
                a_ = b_;
                b_++;
            };
        }
        else
        {
            nb = nb + 1;
            if (b_ <= ind + 1)
            {
                a_ = b_;
            };
        };
    };
};

template <>
void Element<2>::getMatrixC(double **MatrixC)
{
    int DIM = 2;
    int LNN = 18 * DIM - 27;

    int deg[DIM], npc[DIM];
    for (int i = 0; i < DIM; i++)
    {
        deg[i] = (*iparameters)[Npatch_]->getDegree(i);
        npc[i] = (*iparameters)[Npatch_]->getNcp(i);
    };
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    // u direction
    double **matrixCu;
    matrixCu = new double *[deg[0] + 1];
    for (int i = 0; i < deg[0] + 1; ++i)
        matrixCu[i] = new double[deg[0] + 1];

    int dim = deg[0] + npc[0] + 1;
    double uknot_[dim];
    int dir = 0;
    (*iparameters)[Npatch_]->getKnot(dir, dim, uknot_);
    getDirMatrixC(deg[0], dim, inc_[0], uknot_, matrixCu);

    // v direction
    double **matrixCv;
    matrixCv = new double *[deg[1] + 1];
    for (int i = 0; i < deg[1] + 1; ++i)
        matrixCv[i] = new double[deg[1] + 1];

    dim = deg[1] + npc[1] + 1;
    double vknot_[dim];
    dir = 1;
    (*iparameters)[Npatch_]->getKnot(dir, dim, vknot_);
    getDirMatrixC(deg[1], dim, inc_[1], vknot_, matrixCv);

    for (int k = 0; k < deg[1] + 1; k++)
        for (int l = 0; l < deg[1] + 1; l++)
            for (int i = 0; i < deg[0] + 1; i++)
                for (int j = 0; j < deg[0] + 1; j++)
                    MatrixC[i + k * (deg[0] + 1)][j + l * (deg[0] + 1)] = matrixCv[k][l] * matrixCu[i][j];

    for (int i = 0; i < deg[0] + 1; ++i)
        delete[] matrixCu[i];
    delete[] matrixCu;

    for (int i = 0; i < deg[1] + 1; ++i)
        delete[] matrixCv[i];
    delete[] matrixCv;

    return;
}

template <>
void Element<3>::getMatrixC(double **MatrixC)
{
    int DIM = 3;
    int LNN = 18 * DIM - 27;

    int deg[DIM], npc[DIM];
    for (int i = 0; i < DIM; i++)
    {
        deg[i] = (*iparameters)[Npatch_]->getDegree(i);
        npc[i] = (*iparameters)[Npatch_]->getNcp(i);
    };
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    // u direction
    double **matrixCu;
    matrixCu = new double *[deg[0] + 1];
    for (int i = 0; i < deg[0] + 1; ++i)
        matrixCu[i] = new double[deg[0] + 1];

    int dim = deg[0] + npc[0] + 1;
    double uknot_[dim];
    int dir = 0;
    (*iparameters)[Npatch_]->getKnot(dir, dim, uknot_);
    getDirMatrixC(deg[0], dim, inc_[0], uknot_, matrixCu);

    // v direction
    double **matrixCv;
    matrixCv = new double *[deg[1] + 1];
    for (int i = 0; i < deg[1] + 1; ++i)
        matrixCv[i] = new double[deg[1] + 1];
    dim = deg[1] + npc[1] + 1;
    double vknot_[dim];
    dir = 1;
    (*iparameters)[Npatch_]->getKnot(dir, dim, vknot_);
    getDirMatrixC(deg[1], dim, inc_[1], vknot_, matrixCv);

    // t direction
    double **matrixCt;
    matrixCt = new double *[deg[2] + 1];
    for (int i = 0; i < deg[2] + 1; ++i)
        matrixCt[i] = new double[deg[2] + 1];
    dim = deg[2] + npc[2] + 1;
    double tknot_[dim];
    dir = 2;
    (*iparameters)[Npatch_]->getKnot(dir, dim, tknot_);
    getDirMatrixC(deg[2], dim, inc_[2], tknot_, matrixCt);

    double matrixCtv[(deg[2] + 1) * (deg[1] + 1)][(deg[2] + 1) * (deg[1] + 1)];

    for (int k = 0; k < deg[2] + 1; k++)
        for (int l = 0; l < deg[2] + 1; l++)
            for (int i = 0; i < deg[1] + 1; i++)
                for (int j = 0; j < deg[1] + 1; j++)
                    matrixCtv[k * (deg[1] + 1) + i][l * (deg[1] + 1) + j] = matrixCt[k][l] * matrixCv[i][j];

    for (int k = 0; k < (deg[2] + 1) * (deg[1] + 1); k++)
        for (int l = 0; l < (deg[2] + 1) * (deg[1] + 1); l++)
            for (int i = 0; i < deg[0] + 1; i++)
                for (int j = 0; j < deg[0] + 1; j++)
                    MatrixC[k * (deg[0] + 1) + i][l * (deg[0] + 1) + j] = matrixCtv[k][l] * matrixCu[i][j];

    for (int i = 0; i < deg[0] + 1; ++i)
        delete[] matrixCu[i];
    delete[] matrixCu;

    for (int i = 0; i < deg[1] + 1; ++i)
        delete[] matrixCv[i];
    delete[] matrixCv;

    for (int i = 0; i < deg[2] + 1; ++i)
        delete[] matrixCt[i];
    delete[] matrixCt;

    return;
}

template <int DIM>
void Element<DIM>::getInvMatrixC(int &dir, double **MatrixCInv)
{

    int LNN = 18 * DIM - 27;
    int deg = (*iparameters)[Npatch_]->getDegree(dir);
    int npc = (*iparameters)[Npatch_]->getNcp(dir);
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();
    int dim_ = deg + npc + 1;
    double knot[dim_];
    (*iparameters)[Npatch_]->getKnot(dir, dim_, knot);

    // Matrix C
    double **matrixC;
    matrixC = new double *[deg + 1];
    for (int i = 0; i < deg + 1; ++i)
        matrixC[i] = new double[deg + 1];

    getDirMatrixC(deg, dim_, inc_[dir], knot, matrixC);

    MatrixDouble matrixCInvAux(deg + 1, deg + 1);

    for (int i = 0; i < deg + 1; i++)
        for (int j = 0; j < deg + 1; j++)
            matrixCInvAux(i, j) = matrixC[i][j];

    matrixCInvAux = matrixCInvAux.inverse();

    for (int i = 0; i < deg + 1; i++)
        for (int j = 0; j < deg + 1; j++)
            MatrixCInv[i][j] = matrixCInvAux(i, j);

    for (int i = 0; i < deg + 1; ++i)
        delete[] matrixC[i];
    delete[] matrixC;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<int DIM>
void Element<DIM>::getParameterSUPG_FEM(double &tSUPG_, double &tPSPG_, double &tLSIC_, double *phi_, double **dphi_dx) {

	double   r[2] = {}; double s[2] = {};
	
	double tSUGN1_;
	double tSUGN2_;
	double tSUGN3_;
	double hRGN_;
	double &dTime_ = parameters->getTimeStep();
	double &visc_ = parameters->getViscosity();
	double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();


	tSUPG_ = 0.;
	tSUGN1_ = 0.;
	tSUGN2_ = 0.;
	tSUGN3_ = 0.;
	hRGN_ = 0.;
	double hUGN_ = 0.;

	double ua_[2] = {};
	
    for (int i = 0; i < 6; i++){
        double ua[2], uma[2];
        for (int j = 0; j < 2; j++){
            ua[j] = alpha_f * (*nodes_)[connect_[i]]->getVelocity(j) +
                    (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(j);
            uma[j] = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(j) +
                     (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(j);
            ua[j] -= uma[j];
            ua_[j] += ua[j] * phi_[i];
        }
	};

	double uNorm = sqrt(ua_[0] * ua_[0] + ua_[1] * ua_[1]);
	if(uNorm > 1.e-10){
		s[0] = ua_[0] / uNorm;
		s[1] = ua_[1] / uNorm;
	}else{
		s[0] = 1. / sqrt(2.);
		s[1] = 1. / sqrt(2.);
	};

	for (int i = 0; i < 6; i++){

        double ua, va, uma, vma;
        
        ua = alpha_f * (*nodes_)[connect_[i]]->getVelocity(0) +
                (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(0);
        va = alpha_f * (*nodes_)[connect_[i]]->getVelocity(1) +
                (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(1);

        uma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(0) +
                    (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(0);
        vma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(1) +
                    (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(1);

        ua -= uma;
        va -= vma;
    
        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];

	};

	double rNorm = sqrt(r[0]*r[0] + r[1]*r[1]);

	if (rNorm >= 1.e-10){
		r[0] /= rNorm;
		r[1] /= rNorm;
	}else{
		r[0] = 1. / sqrt(2.);
		r[1] = 1. / sqrt(2.);
	};

	for (int i = 0; i < 6; i++){
		hRGN_ += fabs(r[0] * dphi_dx[0][i] + r[1] * dphi_dx[1][i]);
		hUGN_ += fabs(s[0] * dphi_dx[0][i] + s[1] * dphi_dx[1][i]);        
	};

	if (hRGN_ >= 1.e-10){
		hRGN_ = 2. / hRGN_;
	}else{
		hRGN_ = 2. / 1.e-10;
	};

	if (hUGN_ >= 1.e-10){
		hUGN_ = 2. / hUGN_;
	}else{
		hUGN_ = 2. / 1.e-10;
	};    

	if (uNorm >= 1.e-10){
		tSUGN1_ = hUGN_ / (2. * uNorm);
	}else{
		tSUGN1_ = hUGN_ / 2.e-10;
	};
			  
	tSUGN2_ = dTime_ / 2.;

	tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

   
	if (fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
	if (fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;

    //Computing tSUPG parameter
	tSUPG_ = (1. / sqrt(1. / (tSUGN1_ * tSUGN1_) + 
					   1. / (tSUGN2_ * tSUGN2_) + 
					   1. / (tSUGN3_ * tSUGN3_)));
   
	tPSPG_ = tSUPG_;
	tLSIC_ = tSUPG_ * uNorm * uNorm;
 
	return;
};


template <int DIM>
void Element<DIM>::getNewParameterSUPG_FEM(double &tSUPG_, double &tPSPG_, double &tLSIC_, double **Jac, double *phi_, double **dphi_dx)
{

    int LNN = 4 * DIM - 2;

    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();

    MatrixDouble MatrixInvD(DIM, DIM), MatrixInvQh(DIM, DIM), MatrixD(DIM, DIM),
        MatrixQh(DIM, DIM), MatrixG(DIM, DIM);
    double rr[DIM][DIM] = {};
    double r[DIM] = {};
    double ua_[DIM] = {};

    // setting zer
    MatrixD.setZero();
    MatrixQh.setZero();
    MatrixG.setZero();

    double hrqd, tSUGN1_, tSUGN2_, tSUGN3_;

    // Matrix D = polynomial order of the base functions
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            if (i == j)
                MatrixD(i, j) = 2.;

    // Inverse Matrix D
    MatrixInvD = MatrixD.inverse();

    // Matrix Q "hat"
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                MatrixQh(i, j) += Jac[i][k] * MatrixInvD(k, j);

    // Matrix Inverse Q "hat"
    MatrixInvQh = MatrixQh.inverse();

    // Matrix G
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                MatrixG(i, j) += MatrixInvQh(k, i) * MatrixInvQh(k, j);

    // Calculation hrqd
    for (int i = 0; i < LNN; i++)
    {
        double ua[DIM];
        double uma[DIM];

        for (int j = 0; j < DIM; j++)
        {
            ua[j] = alpha_f * (*nodes_)[connect_[i]]->getVelocity(j) +
                    (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(j);
            uma[j] = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(j) +
                     (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(j);
            ua[j] -= uma[j];
            ua_[j] += ua[j] * phi_[i];
        };

        double velnorm = 0.;
        for (int j = 0; j < DIM; j++)
            velnorm += ua[j] * ua[j];
        velnorm = sqrt(velnorm);

        for (int j = 0; j < DIM; j++)
            r[j] += velnorm * dphi_dx[j][i];
    };

    double rNorm = 0;
    for (int i = 0; i < DIM; i++)
        rNorm += r[i] * r[i];
    rNorm = sqrt(rNorm);

     double referpar = 0.000000000000000001;

    r[0] /= (rNorm+referpar);
    r[1] /= (rNorm+referpar);

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            rr[i][j] = r[i] * r[j];

    double rrMG = 0.;
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            rrMG += rr[i][j] * MatrixG(i, j);

    rrMG = sqrt(rrMG);

    hrqd = 2. / (rrMG + referpar);

    Eigen::EigenSolver<MatrixDouble> es(MatrixG);
    std::vector<double> lambda(DIM);
    for (int i = 0; i < DIM; i++)
        lambda[i] = real(es.eigenvalues()[i]);

    auto it = std::minmax_element(lambda.begin(), lambda.end());

    double lambdaMin = *it.first;
    double lambdaMax = *it.second;

    double hmin = 2. / sqrt(lambdaMax);
    double hmax = 2. / sqrt(lambdaMin);

    if (hrqd < hmin)
        hrqd = hmin;
    if (hrqd > hmax)
        hrqd = hmax;

    // tSUGN1_ (-2)
    tSUGN1_ = 0.;
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            tSUGN1_ += ua_[i] * ua_[j] * MatrixG(i, j);

    tSUGN2_ = (dTime_ / 2.);

    // tSUGN3_ (-1)
    tSUGN3_ = (visc_ / dens_) * (4. / (hrqd * hrqd));

    // Computing tSUPG parameter
    tSUPG_ = 1. / (sqrt(tSUGN1_ + (1. / (tSUGN2_ * tSUGN2_)) + (tSUGN3_ * tSUGN3_)));

    tPSPG_ = tSUPG_;

    tLSIC_ = hrqd * hrqd / tSUPG_;

    return;
};

template <int DIM>
void Element<DIM>::getNewParameterSUPG_ISO(double &tSUPG_, double &tPSPG_, double &tLSIC_, double **quadJacMat, double *phi_, double **dphi_dx)
{
    int LNN = 18 * DIM - 27;

    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();

    MatrixDouble MatrixInvD(DIM, DIM), MatrixInvQh(DIM, DIM), MatrixD(DIM, DIM),
        MatrixQh(DIM, DIM), MatrixG(DIM, DIM);
    double rr[DIM][DIM] = {};
    double r[DIM] = {};
    double ua_[DIM] = {};

    // setting zer
    MatrixD.setZero();
    MatrixQh.setZero();
    MatrixG.setZero();

    double hrqd, tSUGN1_, tSUGN2_, tSUGN3_;

    // Bezier parametric element lenght
    double dimxsi = 2.;

    for (int i = 0; i < DIM; i++)
    {

        int deg = (*iparameters)[Npatch_]->getDegree(i);

        double **MatrixCInv;
        MatrixCInv = new double *[deg + 1];
        for (int j = 0; j < deg + 1; ++j)
            MatrixCInv[j] = new double[deg + 1];

        getInvMatrixC(i, MatrixCInv);

        std::vector<double> Dxsi(deg);
        for (int j = 0; j < deg; j++) Dxsi[j] = 0.0;
        
        for (int j = 0; j < deg + 1; j++)
        {
            Dxsi[0] += j * (MatrixCInv[j][1] - MatrixCInv[j][0]);
            Dxsi[1] += j * (MatrixCInv[j][2] - MatrixCInv[j][1]);
        };
    
        for (int j = 0; j < deg; j++) Dxsi[j] *= (dimxsi / deg);

        auto it = std::max_element(std::begin(Dxsi), std::end(Dxsi));
 
        // RQD - MAX
        MatrixD(i, i) = *it;

        for (int j = 0; j < deg + 1; ++j)
            delete[] MatrixCInv[j];
        delete[] MatrixCInv;
    };


    // Inverse Matrix D
    MatrixInvD = MatrixD.inverse();

    // Matrix Q "hat"
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                MatrixQh(i, j) += quadJacMat[i][k] * MatrixInvD(k, j);


    // Matrix Inverse Q "hat"
    MatrixInvQh = MatrixQh.inverse();

    // Matrix G
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                MatrixG(i, j) += MatrixInvQh(k, i) * MatrixInvQh(k, j);

    // Calculation hrqd
    for (int i = 0; i < LNN; i++)
    {
        double ua[DIM];
        double uma[DIM];

        for (int j = 0; j < DIM; j++)
        {
            ua[j] = alpha_f * (*nodes_)[connect_[i]]->getVelocity(j) +
                    (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(j);
            uma[j] = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(j) +
                     (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(j);
            ua[j] -= uma[j];
            ua_[j] += ua[j] * phi_[i];
        };

        double velnorm = 0.;
        for (int j = 0; j < DIM; j++)
            velnorm += ua[j] * ua[j];
        velnorm = sqrt(velnorm);

        for (int j = 0; j < DIM; j++)
            r[j] += velnorm * dphi_dx[j][i];
    };

    double rNorm = 0;
    for (int i = 0; i < DIM; i++)
        rNorm += r[i] * r[i];
    rNorm = sqrt(rNorm);

    double referpar = 0.000000000000000001;

    // physical unitary direction from gradient velocity
    for (int i = 0; i < DIM; i++)
        r[i] /= (rNorm + referpar);

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            rr[i][j] = r[i] * r[j];

    double rrMG = 0.;
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            rrMG += rr[i][j] * MatrixG(i, j);
    

    rrMG = sqrt(rrMG);

    hrqd = 2. / (rrMG + referpar);

    Eigen::EigenSolver<MatrixDouble> es(MatrixG);
    std::vector<double> lambda(DIM);
    for (int i = 0; i < DIM; i++)
        lambda[i] = real(es.eigenvalues()[i]);

    auto it = std::minmax_element(lambda.begin(), lambda.end());

    double lambdaMin = *it.first;
    double lambdaMax = *it.second;

    double hmin = 2. / sqrt(lambdaMax);
    double hmax = 2. / sqrt(lambdaMin);

    if (hrqd < hmin)
        hrqd = hmin;
    if (hrqd > hmax)
        hrqd = hmax;

    // tSUGN1_ (-2)
    tSUGN1_ = 0.;
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            tSUGN1_ += ua_[i] * ua_[j] * MatrixG(i, j);

    tSUGN2_ = dTime_ / 2.;

    // tSUGN3_ (-1)
    tSUGN3_ = (visc_ / dens_) * (4. / (hrqd * hrqd));

    // Computing tSUPG parameter
    tSUPG_ = 1. / (sqrt(tSUGN1_ + (1. / (tSUGN2_ * tSUGN2_)) + (tSUGN3_ * tSUGN3_)));

    tPSPG_ = tSUPG_;

    tLSIC_ = hrqd * hrqd / tSUPG_;

    return;
};

template <int DIM>
void Element<DIM>::getParameterArlequinMN(int &LNN, int &LNNC, double &wna_, double &djac_, double &weight_,
                                          double &tARLQ_, double *phi_, double *phiC_,
                                          double **dphi_dx, double ***ddphi_dx)
{

    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh Velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVel(LNN, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

    // Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duMesh_dx[i] = new double[DIM];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duMeshPrev_dx[i] = new double[DIM];

    double dumeshna_dx[DIM][DIM];
    getInterpMeshVelDer(LNN, dphi_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];

    double lambda1[DIM * LNN][DIM * LNN] = {};
    double lambda0[DIM * LNN][DIM * LNNC] = {};
    double convec[DIM * LNN][DIM * LNN] = {};
    double iner[DIM * LNN][DIM * LNN] = {};
    double visc[DIM * LNN][DIM * LNN] = {};
    double press[DIM * LNN][DIM * LNN] = {};
    double gradlambda[DIM * LNN][DIM * LNN] = {};

    double AGDT = alpha_f * gamma * dTime_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {

            double inertia = 0.;
            double convection = 0.;
            double gradientLambda = 0.;
            double aux0 = 0.;

            for (int k = 0; k < DIM; k++)
            {

                inertia += dphi_dx[k][i] * dphi_dx[k][j];
                gradientLambda += dphi_dx[k][i] * dphi_dx[k][j];
                aux0 += ddphi_dx[k][k][i];

                for (int l = 0; l < DIM; l++)
                    convection += dphi_dx[k][i] * ((una_[l] - umeshna_[l]) * ddphi_dx[l][k][j] +
                                                   (duna_dx[l][k] - dumeshna_dx[l][k]) * dphi_dx[l][j]);
            };

            for (int k = 0; k < DIM; k++)
            {

                double pression = 0.;
                for (int l = 0; l < DIM; l++)
                {

                    pression += dphi_dx[l][i] * ddphi_dx[l][k][j];

                    double viscosity = ddphi_dx[k][l][j];
                    if (k == l)
                        for (int m = 0; m < DIM; m++)
                            viscosity += ddphi_dx[m][m][j];

                    visc[DIM * i + k][DIM * j + l] += aux0 * viscosity * visc_ * AGDT * wna_;
                };

                iner[DIM * i + k][DIM * j + k] += inertia * alpha_m * wna_;
                convec[DIM * i + k][DIM * j + k] += convection * AGDT * wna_;
                press[DIM * i + k][j] += pression / dens_ * wna_;
                gradlambda[DIM * i + k][DIM * j + k] += gradientLambda * k1 / dens_;
                lambda1[DIM * i + k][DIM * j + k] += phi_[i] * phi_[j] * k1;
            };
        };

        for (int j = 0; j < LNNC; j++)
            for (int k = 0; k < DIM; k++)
                lambda0[DIM * i + k][DIM * j + k] += phi_[i] * phiC_[j] * k1;
    };

    double normLambda0 = 0.;
    double normLambda1 = 0.;
    double normConvec = 0.;
    double normIner = 0.;
    double normVisc = 0.;
    double normPress = 0.;
    double normGradLambda = 0.;

    for (int i = 0; i < DIM * LNN; i++)
    {
        for (int j = 0; j < DIM * LNN; j++)
        {
            normLambda1 += lambda1[i][j] * lambda1[i][j];
            normConvec += convec[i][j] * convec[i][j];
            normIner += iner[i][j] * iner[i][j];
            normVisc += visc[i][j] * visc[i][j];
            normPress += press[i][j] * press[i][j];
            normGradLambda += gradlambda[i][j] * gradlambda[i][j];
        }
        for (int j = 0; j < DIM * LNNC; j++)
        {
            normLambda0 += lambda0[i][j] * lambda0[i][j];
        };
    };

    normLambda0 = sqrt(normLambda0);
    normLambda1 = sqrt(normLambda1);
    normConvec = sqrt(normConvec);
    normIner = sqrt(normIner);
    normVisc = sqrt(normVisc);
    normPress = sqrt(normPress);
    normGradLambda = sqrt(normGradLambda);

    if (fabs(normIner) <= 1.e-10)
        normIner = 1.e-10;
    if (fabs(normVisc) <= 1.e-10)
        normVisc = 1.e-10;
    if (fabs(normPress) <= 1.e-10)
        normPress = 1.e-10;
    if (fabs(normConvec) <= 1.e-10)
        normConvec = 1.e-10;
    if (fabs(normGradLambda) <= 1.e-10)
        normGradLambda = 1.e-10;

    double tA1 = std::fabs(normLambda1) / std::fabs(normConvec);
    double tA2 = std::fabs(normLambda0) / std::fabs(normConvec);
    double tB1 = std::fabs(normLambda1) / std::fabs(normIner);
    double tB2 = std::fabs(normLambda0) / std::fabs(normIner);
    double tC1 = std::fabs(normLambda1) / std::fabs(normVisc);
    double tC2 = std::fabs(normLambda0) / std::fabs(normVisc);
    double tD1 = std::fabs(normLambda1) / std::fabs(normPress);
    double tD2 = std::fabs(normLambda0) / std::fabs(normPress);
    double tE1 = std::fabs(normLambda1) / std::fabs(normGradLambda);
    double tE2 = std::fabs(normLambda0) / std::fabs(normGradLambda);
    if (tE1 > 1.)
        tE1 = 0.0000000001;
    if (tE2 > 1.)
        tE2 = 0.0000000001;

    double tA = 1. / sqrt(1. / (tA1 * tA1) + 1. / (tA2 * tA2));

    double tB = 1. / sqrt(1. / (tB1 * tB1) + 1. / (tB2 * tB2));

    double tC = 1. / sqrt(1. / (tC1 * tC1) + 1. / (tC2 * tC2));

    double tD = 1. / sqrt(1. / (tD1 * tD1) + 1. / (tD2 * tD2));

    double tE = 1. / sqrt(1. / (tE1 * tE1) + 1. / (tE2 * tE2));

    tARLQ_ = -1. / sqrt(1. / (tA * tA) +
                        1. / (tB * tB) +
                        1. / (tC * tC) +
                        1. / (tE * tE) +
                        1. / (tD * tD));

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;
};

template <int DIM>
void Element<DIM>::getParameterArlequinLaplaceMN(int &LNN, int &LNNC, double &wna_, double &djac_, double &weight_, double &tARLQ0_, double &tARLQ1_,
                                                 double *phi_, double *phiC_, double **dphi_dx, double **dphiC_dx, double ***ddphi_dx, double ***ddphiC_dx)
{

    double &k1 = parameters->getArlequinK1();

    // //Interpolated variables

    // // Fine velocity
    // double u_[DIM], uPrev_[DIM];
    // getInterpVel(LNN, phi_, u_, uPrev_);

    // // Coarse velocity
    // double uC_[DIM], uPrevC_[DIM];
    // getInterpVelCoarse(LNNC, phiC_, uC_, uPrevC_);

    // Fine velocity derivatives
    // double **du_dx;
    // du_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     du_dx[i] = new double[DIM];

    // double **duPrev_dx;
    // duPrev_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     duPrev_dx[i] = new double[DIM];

    // getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);

    // Coarse velocity derivatives
    // double **duC_dx;
    // duC_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     duC_dx[i] = new double[DIM];

    // double **duCPrev_dx;
    // duCPrev_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     duCPrev_dx[i] = new double[DIM];

    // getInterpVelDerCoarse(LNNC, dphiC_dx,duC_dx,duCPrev_dx);

    // Fine second velocity derivatives
    // double ***ddu_dxdx;
    // ddu_dxdx = new double **[DIM];
    // for (int i = 0; i < DIM; ++i)
    // {
    //     ddu_dxdx[i] = new double *[DIM];
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         ddu_dxdx[i][j] = new double[DIM];
    //     };
    // };

    // double ***dduPrev_dxdx;
    // dduPrev_dxdx = new double **[DIM];
    // for (int i = 0; i < DIM; ++i)
    // {
    //     dduPrev_dxdx[i] = new double *[DIM];
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         dduPrev_dxdx[i][j] = new double[DIM];
    //     };
    // };

    // getInterpSecondVelDer(LNN, ddphi_dx, ddu_dxdx,dduPrev_dxdx);

    // Coarse second velocity derivatives
    // double ***dduC_dxdx;
    // dduC_dxdx = new double **[DIM];
    // for (int i = 0; i < DIM; ++i)
    // {
    //     dduC_dxdx[i] = new double *[DIM];
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         dduC_dxdx[i][j] = new double[DIM];
    //     };
    // };

    // double ***dduCPrev_dxdx;
    // dduCPrev_dxdx = new double **[DIM];
    // for (int i = 0; i < DIM; ++i)
    // {
    //     dduCPrev_dxdx[i] = new double *[DIM];
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         dduCPrev_dxdx[i][j] = new double[DIM];
    //     };
    // };
    // getInterpSecondVelDerCoarse(LNNC, ddphiC_dx, dduC_dxdx, dduCPrev_dxdx);

    // // Lambda Derivative
    // double **dlambda_dx;
    // dlambda_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i) dlambda_dx[i] = new double[DIM];
    // getInterpLambdaDer(LNN,dphi_dx,dlambda_dx);

    double FC1[LNN][LNN] = {};
    double FE[LNN][LNN] = {};
    double FC0[LNN][LNNC] = {};

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            double gradientLambda = 0.;
            for (int k = 0; k < DIM; k++)
                gradientLambda += dphi_dx[k][i] * dphi_dx[k][j];

            FE[i][j] += gradientLambda * k1;
            FC1[i][j] += phi_[i] * phi_[j] * k1;
        };
        for (int j = 0; j < LNNC; j++)
            FC0[i][j] += phi_[i] * phiC_[j] * k1;
    };

    double normFC1 = 0.;
    double normFE = 0.;
    double normFC0 = 0.;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            normFC1 += FC1[i][j] * FC1[i][j];
            normFE += FE[i][j] * FE[i][j];
        }
        for (int j = 0; j < LNNC; j++)
            normFC0 += FC0[i][j] * FC0[i][j];
    };

    normFC1 = sqrt(normFC1);
    normFE = sqrt(normFE);
    normFC0 = sqrt(normFC0);

    if (fabs(normFE) <= 1.e-10)
        normFE = 1.e-10;
    if (fabs(normFC1) <= 1.e-10)
        normFC1 = 1.e-10;
    if (fabs(normFC0) <= 1.e-10)
        normFC0 = 1.e-10;

    // Option 1
    tARLQ1_ = std::fabs(normFC1) / std::fabs(normFE);
    tARLQ0_ = std::fabs(normFC0) / std::fabs(normFE);

    // Option 2
    //  tARLQ0_ = std::min({std::fabs(normFC0)/std::fabs(normFE),std::fabs(normFC1)/std::fabs(normFE)});
    //  tARLQ1_ = tARLQ0_;

    // Option 3
    //  double t1 = std::fabs(normFC1)/std::fabs(normFE);
    //  double fc = std::min({std::fabs(normFC0),std::fabs(normFC1)});
    //  double t0 = std::fabs(normFC0)/std::fabs(normFE);

    // tARLQ1_ = t1*fc/std::fabs(normFC1);
    // tARLQ0_ = t0*fc/std::fabs(normFC0);

    // Option 4
    //  tARLQ1_ = std::fabs(normFC1)/std::fabs(normFE);
    //  tARLQ0_ = 0.;

    // Option 5
    //  double fc = std::min({std::fabs(normFC0),std::fabs(normFC1)});
    //  tARLQ0_ = 0.;
    //  tARLQ1_ = fc/std::fabs(normFE);

    // // deallocating memory

    // for (int i = 0; i < DIM; ++i)
    //     delete[] du_dx[i];
    // delete[] du_dx;

    // for (int i = 0; i < DIM; ++i)
    //     delete[] duC_dx[i];
    // delete[] duC_dx;

    // for (int i = 0; i < DIM; ++i)
    //     delete[] duPrev_dx[i];
    // delete[] duPrev_dx;

    // for (int i = 0; i < DIM; ++i)
    //     delete[] duCPrev_dx[i];
    // delete[] duCPrev_dx;

    // for (int i = 0; i < DIM; ++i)
    // {
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         delete[] ddu_dxdx[i][j];
    //     };
    //     delete[] ddu_dxdx[i];
    // };
    // delete[] ddu_dxdx;

    // for (int i = 0; i < DIM; ++i)
    // {
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         delete[] dduPrev_dxdx[i][j];
    //     };
    //     delete[] dduPrev_dxdx[i];
    // };
    // delete[] dduPrev_dxdx;

    // for (int i = 0; i < DIM; ++i)
    // {
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         delete[] dduC_dxdx[i][j];
    //     };
    //     delete[] dduC_dxdx[i];
    // };
    // delete[] dduC_dxdx;

    // for (int i = 0; i < DIM; ++i)
    // {
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         delete[] dduCPrev_dxdx[i][j];
    //     };
    //     delete[] dduCPrev_dxdx[i];
    // };
    // delete[] dduCPrev_dxdx;

    // for (int i = 0; i < DIM; ++i)
    //     delete[] dlambda_dx[i];
    // delete[] dlambda_dx;
};

template <int DIM>
void Element<DIM>::getParameterArlequinMN_COARSE(int &LNN, int &LNNC, double &wna_, double &djac_, double &weight_, double &tARLQ_,
                                                 double *phi_, double *phiC_, double **dphi_dx, double **dphiC_dx, double ***ddphi_dx,
                                                 double ***ddphiC_dx)
{

    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVelCoarse(LNNC, phiC_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh Velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVelCoarse(LNNC, phiC_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDerCoarse(LNNC, dphiC_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    // Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duMesh_dx[i] = new double[DIM];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duMeshPrev_dx[i] = new double[DIM];

    double dumeshna_dx[DIM][DIM];
    getInterpMeshVelDerCoarse(LNNC, dphiC_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];
        };
    };

    double lambda1[DIM * LNN][DIM * LNN] = {};
    double lambda0[DIM * LNN][DIM * LNNC] = {};
    double convec[DIM * LNN][DIM * LNNC] = {};
    double iner[DIM * LNN][DIM * LNNC] = {};
    double visc[DIM * LNN][DIM * LNNC] = {};
    double press[DIM * LNN][DIM * LNNC] = {};
    double gradlambda[DIM * LNN][DIM * LNN] = {};

    double AGDT = alpha_f * gamma * dTime_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNNC; j++)
        {

            double inertia = 0.;
            double convection = 0.;
            double aux0 = 0.;

            for (int k = 0; k < DIM; k++)
            {

                inertia += dphi_dx[k][i] * dphiC_dx[k][j];
                aux0 += ddphi_dx[k][k][i];

                for (int l = 0; l < DIM; l++)
                    convection += dphi_dx[k][i] * ((una_[l] - umeshna_[l]) * ddphiC_dx[l][k][j] +
                                                   (duna_dx[l][k] - dumeshna_dx[l][k]) * dphiC_dx[l][j]);
            };

            for (int k = 0; k < DIM; k++)
            {

                double pression = 0.;
                for (int l = 0; l < DIM; l++)
                {

                    pression += dphi_dx[l][i] * ddphiC_dx[l][k][j];

                    double viscosity = ddphiC_dx[k][l][j];
                    if (k == l)
                        for (int m = 0; m < DIM; m++)
                            viscosity += ddphiC_dx[m][m][j];

                    visc[DIM * i + k][DIM * j + l] += aux0 * viscosity * visc_ * AGDT * wna_;
                };

                iner[DIM * i + k][DIM * j + k] += inertia * alpha_m * wna_;
                convec[DIM * i + k][DIM * j + k] += convection * AGDT * wna_;
                press[DIM * i + k][j] += pression / dens_ * wna_;
                lambda0[DIM * i + k][DIM * j + k] += phi_[i] * phiC_[j] * k1;
            };
        };

        for (int j = 0; j < LNN; j++)
        {

            double gradientLambda = 0.;
            for (int k = 0; k < DIM; k++)
                gradientLambda += dphi_dx[k][i] * dphi_dx[k][j];

            for (int k = 0; k < DIM; k++)
            {
                lambda1[DIM * i + k][DIM * j + k] += phi_[i] * phi_[j] * k1;
                gradlambda[DIM * i + k][DIM * j + k] += gradientLambda * k1 / dens_;
            };
        };
    };

    double normLambda0 = 0.;
    double normLambda1 = 0.;
    double normConvec = 0.;
    double normIner = 0.;
    double normVisc = 0.;
    double normPress = 0.;
    double normGradLambda = 0.;

    for (int i = 0; i < DIM * LNN; i++)
    {
        for (int j = 0; j < DIM * LNNC; j++)
        {
            normLambda0 += lambda0[i][j] * lambda0[i][j];
            normConvec += convec[i][j] * convec[i][j];
            normIner += iner[i][j] * iner[i][j];
            normVisc += visc[i][j] * visc[i][j];
            normPress += press[i][j] * press[i][j];
        }
        for (int j = 0; j < DIM * LNN; j++)
        {
            normLambda1 += lambda1[i][j] * lambda1[i][j];
            normGradLambda += gradlambda[i][j] * gradlambda[i][j];
        };
    };

    normLambda0 = sqrt(normLambda0);
    normLambda1 = sqrt(normLambda1);
    normConvec = sqrt(normConvec);
    normIner = sqrt(normIner);
    normVisc = sqrt(normVisc);
    normPress = sqrt(normPress);
    normGradLambda = sqrt(normGradLambda);

    if (fabs(normIner) <= 1.e-10)
        normIner = 1.e-10;
    if (fabs(normVisc) <= 1.e-10)
        normVisc = 1.e-10;
    if (fabs(normPress) <= 1.e-10)
        normPress = 1.e-10;
    if (fabs(normConvec) <= 1.e-10)
        normConvec = 1.e-10;
    if (fabs(normGradLambda) <= 1.e-10)
        normGradLambda = 1.e-10;

    double tA1 = std::fabs(normLambda1) / std::fabs(normConvec);
    double tA2 = std::fabs(normLambda0) / std::fabs(normConvec);
    double tB1 = std::fabs(normLambda1) / std::fabs(normIner);
    double tB2 = std::fabs(normLambda0) / std::fabs(normIner);
    double tC1 = std::fabs(normLambda1) / std::fabs(normVisc);
    double tC2 = std::fabs(normLambda0) / std::fabs(normVisc);
    double tD1 = std::fabs(normLambda1) / std::fabs(normPress);
    double tD2 = std::fabs(normLambda0) / std::fabs(normPress);
    double tE1 = std::fabs(normLambda1) / std::fabs(normGradLambda);
    double tE2 = std::fabs(normLambda0) / std::fabs(normGradLambda);
    if (tE1 > 1.)
        tE1 = 0.0000000001;
    if (tE2 > 1.)
        tE2 = 0.0000000001;

    double tA = 1. / sqrt(1. / (tA1 * tA1) + 1. / (tA2 * tA2));

    double tB = 1. / sqrt(1. / (tB1 * tB1) + 1. / (tB2 * tB2));

    double tC = 1. / sqrt(1. / (tC1 * tC1) + 1. / (tC2 * tC2));

    double tD = 1. / sqrt(1. / (tD1 * tD1) + 1. / (tD2 * tD2));

    double tE = 1. / sqrt(1. / (tE1 * tE1) + 1. / (tE2 * tE2));

    tARLQ_ = -1. / sqrt(1. / (tA * tA) +
                        1. / (tB * tB) +
                        1. / (tC * tC) +
                        1. / (tE * tE) +
                        1. / (tD * tD));

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;
};

//------------------------------------------------------------------------------
//-----------------------------ELEMENT LOCAL MATRIX-----------------------------
//------------------------------------------------------------------------------

template <int DIM>
void Element<DIM>::getElemMatrix(int &LNN, double &wna_, double &djac_, double &weight_, double &tSUPG_, double &tPSPG_, double &tLSIC_,
                                 double *phi_, double **dphi_dx, double **jacobianNRMatrix)
{
    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh Velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVel(LNN, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

    double WJ = weight_ * djac_ * wna_;
    double AGDT = alpha_f * gamma * dTime_;
    double VAGDT = visc_ * AGDT;
    double DAGDT = dens_ * AGDT;
    double DAM = dens_ * alpha_m;

    for (int i = 0; i < LNN; i++)
    {
        double wSUPGi = 0.;
        for (int k = 0; k < DIM; k++)
            wSUPGi += (una_[k] - umeshna_[k]) * dphi_dx[k][i];

        for (int j = 0; j < LNN; j++)
        {

            double wSUPGj = 0.;
            for (int k = 0; k < DIM; k++)
                wSUPGj += (una_[k] - umeshna_[k]) * dphi_dx[k][j];

            // Mass matrix
            double M = (phi_[i] * phi_[j] + wSUPGi * phi_[j] * tSUPG_) * DAM;

            // Convection matrix
            double C = (wSUPGj * phi_[i] + wSUPGi * wSUPGj * tSUPG_) * DAGDT;

            for (int k = 0; k < DIM; k++)
            {

                jacobianNRMatrix[DIM * i + k][DIM * j + k] += (M + C) * WJ;

                for (int l = 0; l < DIM; l++)
                {

                    // Convection derivatives
                    double Cuu = phi_[i] * duna_dx[k][l] * phi_[j] * DAGDT;
                    // LSIC
                    double KLS = dphi_dx[k][i] * dphi_dx[l][j] * tLSIC_ * DAGDT;
                    // Difusion matrix
                    double kk = dphi_dx[l][i] * dphi_dx[k][j] * VAGDT;
                    if (k == l)
                        for (int m = 0; m < DIM; m++)
                            kk += dphi_dx[m][i] * dphi_dx[m][j] * VAGDT;

                    jacobianNRMatrix[DIM * i + k][DIM * j + l] += (Cuu + KLS + kk) * WJ;
                };

                // Gradient operator
                double QSUPG = -(dphi_dx[k][i] * phi_[j]) + wSUPGi * dphi_dx[k][j] * tSUPG_;

                // divergent operator
                double QQ = dphi_dx[k][j] * phi_[i] * AGDT;

                jacobianNRMatrix[DIM * i + k][DIM * LNN + j] += QSUPG * WJ;
                jacobianNRMatrix[DIM * LNN + i][DIM * j + k] += QQ * WJ;

                // PSPG stabilization matrixes
                double H = dphi_dx[k][i] * phi_[j] * tPSPG_ * alpha_m;
                double G = dphi_dx[k][i] * wSUPGj * tPSPG_ * AGDT;

                jacobianNRMatrix[DIM * LNN + i][DIM * j + k] += (H + G) * WJ;
            };

            double Q = 0.;
            for (int m = 0; m < DIM; m++)
                Q += dphi_dx[m][i] * dphi_dx[m][j] * tPSPG_ / (dens_);

            jacobianNRMatrix[DIM * LNN + i][DIM * LNN + j] += Q * WJ;
        };
    };

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    return;
};

template <int DIM>
void Element<DIM>::getElemMatrixLaplace(int &LNN, double &wna_, double &djac_, double &weight_,
                                        double **dphi_dx, double **jacobianNRMatrix)
{

    double WJ = djac_ * weight_ * wna_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            double Lap = 0.;
            for (int k = 0; k < DIM; k++)
                Lap += dphi_dx[k][i] * dphi_dx[k][j];
            jacobianNRMatrix[i][j] += Lap * WJ;
        };
    };

    return;
};

template <int DIM>
void Element<DIM>::getElemMatrixElas(int &LNN, double &djac_, double &weight_, double &meshMovPar,
                                     double **dphi_dx, double **jacobianNRMatrix)
{

    double WJ = djac_ * weight_ * meshMovPar;

    // Elasticity version
    double E_ = 1.;
    double v_ = 0.3;
    double const1 = E_ * v_ / ((1. + v_) * (1. - 2. * v_));
    double const2 = E_ / (1. + v_);

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    double kk = dphi_dx[k][i] * dphi_dx[l][j] * const1 + 0.5 * dphi_dx[l][i] * dphi_dx[k][j] * const2;
                   
                    if (k==l)
                    for (int m = 0; m < DIM; m++) kk += 0.5 * dphi_dx[m][i] * dphi_dx[m][j] * const2;

                    jacobianNRMatrix[DIM * i + k][DIM * j + l] += kk * WJ;
                };
            };
        };
    };

    // for (int i = 0; i < 6; i++)
    // {
    //     for (int j = 0; j < 6; j++)
    //     {

    //         // First term
    //         jacobianNRMatrix[2 * i][2 * j] += dphi_dx[0][i] * dphi_dx[0][j] * weight_ * djac_ * meshMovingParameter * const1;
    //         jacobianNRMatrix[2 * i][2 * j + 1] += dphi_dx[0][i] * dphi_dx[1][j] * weight_ * djac_ * meshMovingParameter * const1;
    //         jacobianNRMatrix[2 * i + 1][2 * j] += dphi_dx[1][i] * dphi_dx[0][j] * weight_ * djac_ * meshMovingParameter * const1;
    //         jacobianNRMatrix[2 * i + 1][2 * j + 1] += dphi_dx[1][i] * dphi_dx[1][j] * weight_ * djac_ * meshMovingParameter * const1;

    //         // Second term
    //         jacobianNRMatrix[2 * i][2 * j] += (dphi_dx[0][i] * dphi_dx[0][j] + 0.5 * dphi_dx[1][i] * dphi_dx[1][j]) * weight_ * djac_ * meshMovingParameter * const2;
    //         jacobianNRMatrix[2 * i][2 * j + 1] += 0.5 * dphi_dx[1][i] * dphi_dx[0][j] * weight_ * djac_ * meshMovingParameter * const2;
    //         jacobianNRMatrix[2 * i + 1][2 * j] += 0.5 * dphi_dx[0][i] * dphi_dx[1][j] * weight_ * djac_ * meshMovingParameter * const2;
    //         jacobianNRMatrix[2 * i + 1][2 * j + 1] += (0.5 * dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * weight_ * djac_ * meshMovingParameter * const2;
    //     };
    // };

    return;
};

template <int DIM>
void Element<DIM>::getElemMatrixTNS(int &LNN, double &wna_, double &djac_, double &weight_,
                                    double &tSUPG_, double &tPSPG_, double &tLSIC_,
                                    double *phi_, double **dphi_dx, double **jacobianNRMatrix)
{

    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh Velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVel(LNN, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

    double WJ = weight_ * djac_ * wna_;
    double AGDT = alpha_f * gamma * dTime_;
    double VAGDT = visc_ * AGDT;
    double DAGDT = dens_ * AGDT;
    double DAM = dens_ * alpha_m;

    // Navier Stokes Problem
    for (int i = 0; i < LNN; i++)
    {
        double wSUPGi = 0.;
        for (int k = 0; k < DIM; k++)
            wSUPGi += (una_[k] - umeshna_[k]) * dphi_dx[k][i];

        for (int j = 0; j < LNN; j++)
        {
 
            double wSUPGj = 0.;
            for (int k = 0; k < DIM; k++)
                wSUPGj += (una_[k] - umeshna_[k]) * dphi_dx[k][j];

            // Mass matrix
            double M = (phi_[i] * phi_[j] + wSUPGi * phi_[j] * tSUPG_) * DAM;

            // Convection matrix
            double C = (wSUPGj * phi_[i] + wSUPGi * wSUPGj * tSUPG_) * DAGDT;

            for (int k = 0; k < DIM; k++)
            {
 
                jacobianNRMatrix[DIM * i + k][DIM * j + k] += (M + C) * WJ;

                for (int l = 0; l < DIM; l++)
                {

                    // Convection derivatives
                    double Cuu = phi_[i] * duna_dx[k][l] * phi_[j] * DAGDT;
                    // LSIC
                    double KLS = dphi_dx[k][i] * dphi_dx[l][j] * tLSIC_ * DAGDT;
                    // Difusion matrix
                    double kk = dphi_dx[l][i] * dphi_dx[k][j] * VAGDT;
                    if (k == l)
                        for (int m = 0; m < DIM; m++)
                            kk += dphi_dx[m][i] * dphi_dx[m][j] * VAGDT;

                    jacobianNRMatrix[DIM * i + k][DIM * j + l] += (Cuu + KLS + kk) * WJ;
                };

                // Gradient operator
                double QSUPG = -(dphi_dx[k][i] * phi_[j]) + wSUPGi * dphi_dx[k][j] * tSUPG_;

                // divergent operator
                double QQ = dphi_dx[k][j] * phi_[i] * AGDT;

                jacobianNRMatrix[DIM * i + k][DIM * LNN + j] += QSUPG * WJ;
                jacobianNRMatrix[DIM * LNN + i][DIM * j + k] += QQ * WJ;

                // PSPG stabilization matrixes
                double H = dphi_dx[k][i] * phi_[j] * tPSPG_ * alpha_m;
                double G = dphi_dx[k][i] * wSUPGj * tPSPG_ * AGDT;

                jacobianNRMatrix[DIM * LNN + i][DIM * j + k] += (H + G) * WJ;
            };

            double Q = 0.;
            for (int m = 0; m < DIM; m++)
                Q += dphi_dx[m][i] * dphi_dx[m][j] * tPSPG_ / (dens_);

            jacobianNRMatrix[DIM * LNN + i][DIM * LNN + j] += Q * WJ;
        };
    };

    // //Stokes Problem
    // for (int i = 0; i < LNN; i++)
    // {
    //     for (int j = 0; j < LNN; j++)
    //     {
    //         double M = phi_[i] * phi_[j] * DAM;

    //         for (int k = 0; k < DIM; k++)
    //         {
    //             jacobianNRMatrix[DIM * i + k][DIM * j + k] += M * WJ;

    //             for (int l = 0; l < DIM; l++)
    //             {
    //                 // Difusion matrix
    //                 double kk = dphi_dx[l][i] * dphi_dx[k][j] * VAGDT;
    //                 if (k == l)
    //                     for (int m = 0; m < DIM; m++)
    //                         kk += dphi_dx[m][i] * dphi_dx[m][j] * VAGDT;

    //                 jacobianNRMatrix[DIM * i + k][DIM * j + l] += kk * WJ;
    //             };

    //             // Gradient operator
    //             double QSUPG = -dphi_dx[k][i] * phi_[j];

    //             // divergent operator
    //             double QQ = dphi_dx[k][j] * phi_[i] * AGDT;

    //             jacobianNRMatrix[DIM * i + k][DIM * LNN + j] += QSUPG * WJ;
    //             jacobianNRMatrix[DIM * LNN + i][DIM * j + k] += QQ * WJ;

    //             // PSPG stabilization matrixes
    //             double H = dphi_dx[k][i] * phi_[j] * tPSPG_ * alpha_m;

    //             jacobianNRMatrix[DIM * LNN + i][DIM * j + k] += H * WJ;

    //         };

    //         double Q = 0.;
    //         for (int m = 0; m < DIM; m++) Q += dphi_dx[m][i] * dphi_dx[m][j] * tPSPG_ / (dens_);

    //         jacobianNRMatrix[DIM * LNN + i][DIM * LNN + j] += Q * WJ;
    //     };
    // };

    // Laplace Problem -> A atualizacao da velocidade fica diferente para esse problema.
    //  for (int i = 0; i < LNN; i++){
    //      for (int j = 0; j < LNN; j++){
    //          double Lap = 0.;
    //          for (int k = 0; k < DIM; k++) Lap += dphi_dx[k][i] * dphi_dx[k][j];
    //          for (int k = 0; k < DIM; k++) jacobianNRMatrix[DIM*i+k][DIM*j+k] += Lap * WJ;
    //      };
    //  };

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    return;
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsSameMesh(int &LNN, double &djac_, double &weight_, double *phi_, double **dphi_dx,
                                               double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    // // Fluid Data
    // double &k1 = parameters->getArlequinK1();
    // double &alpha_f = parameters->getAlphaF();

    // // Interpolated variables

    // // Velocity
    // double u_[DIM], uPrev_[DIM], una_[DIM];
    // getInterpVel(LNN, phi_, u_, uPrev_);
    // for (int i = 0; i < DIM; i++)
    // una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // // Lagrange Multiplieres
    // double lambda_[DIM];
    // getInterpLambda(LNN, phi_, lambda_);

    // double WJ = weight_ * djac_;

    // for (int i = 0; i < LNN; i++)
    // {
    //     for (int j = 0; j < LNN; j++)
    //     {
    //         // L2 operator
    //         double L2 = -phi_[i] * phi_[j] * k1;
    //         for (int k = 0; k < DIM; k++)
    //         lagrMultMatrix[DIM*i+k][DIM*j+k] += L2 * WJ;

    //     };

    //     for (int k = 0; k < DIM; k++){

    //         double L2L = phi_[i] * lambda_[k] * k1;
    //         double L2U = phi_[i] * una_[k] * k1;

    //         rhsVector1[DIM*i+k] += (L2L) * WJ;
    //         rhsVector2[DIM*i+k] += (L2U) * WJ;
    //     };

    // };

    // Fluid Data
    double &alpha_f = parameters->getAlphaF();
    double &dens_ = parameters->getDensity();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    // Interpolated variables

    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

    // Lagrange Multiplieres
    double lambda_[DIM];
    getInterpLambda(LNN, phi_, lambda_);

    // Lagrange Multiplieres derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dlambda_dx[i] = new double[DIM];
    getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

    double WJ = weight_ * djac_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            // L2 operator
            double L2 = -phi_[i] * phi_[j] * k1;

            for (int k = 0; k < DIM; k++)
            {

                lagrMultMatrix[DIM * i + k][DIM * j + k] += L2 * WJ;

                for (int l = 0; l < DIM; l++)
                {
                    double H1 = -0.5 * dphi_dx[l][i] * dphi_dx[k][j] * k2;
                    if (l == k)
                        for (int m = 0; m < DIM; m++)
                            H1 += -0.5 * dphi_dx[m][i] * dphi_dx[m][j] * k2;
                    lagrMultMatrix[DIM * i + k][DIM * j + l] += H1 * WJ;
                };
            };
        };

        for (int k = 0; k < DIM; k++)
        {

            double L2L = phi_[i] * lambda_[k] * k1;
            double L2U = phi_[i] * una_[k] * k1;

            double H1L = 0.;
            double H1U = 0.;
            for (int l = 0; l < DIM; l++)
            {
                H1L += 0.5 * (dphi_dx[l][i] * dlambda_dx[l][k] + dphi_dx[l][i] * dlambda_dx[k][l]) * k2;
                H1U += 0.5 * (dphi_dx[l][i] * duna_dx[l][k] + dphi_dx[l][i] * duna_dx[k][l]) * k2;
            };

            rhsVector1[DIM * i + k] += (L2L + H1L) * WJ;
            rhsVector2[DIM * i + k] += (L2U + H1U) * WJ;
        };
    };

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] dlambda_dx[i];
    delete[] dlambda_dx;
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsSameMeshLaplace(int &LNN, double &djac_, double &weight_, double *phi_,
                                                      double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    // Fluid Data
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables

    // Velocity
    double u_[DIM], uPrev_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);

    // Lagrange Multiplieres
    double lambda_[DIM];
    getInterpLambda(LNN, phi_, lambda_);

    double WJ = weight_ * djac_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            // L2 operator
            double L2 = -phi_[i] * phi_[j] * k1;
            lagrMultMatrix[i][j] += L2 * WJ;
        };

        double L2L = phi_[i] * lambda_[0] * k1;
        double L2U = phi_[i] * u_[0] * k1;

        rhsVector1[i] += (L2L)*WJ;
        rhsVector2[i] += (L2U)*WJ;
    };
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsSameMesh_tSUPG_tPSPG(int &LNN, double &djac_, double &weight_, double &tSUPG_, double &tPSPG_,
                                                           double *phi_, double **dphi_dx, double **jacobianNRMatrix, double *rhsVector)
{

    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh Velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVel(LNN, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Lagrange Multiplieres
    double lambda_[DIM];
    getInterpLambda(LNN, phi_, lambda_);

    double WJ = weight_ * djac_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            double msupg = 0.;
            for (int k = 0; k < DIM; k++)
                msupg += -(una_[k] - umeshna_[k]) * dphi_dx[k][i] * phi_[j] * tSUPG_ * k1;

            for (int k = 0; k < DIM; k++)
            {
                double mpspg = -(dphi_dx[k][i] * phi_[j]) * tPSPG_ * k1 / dens_;
                jacobianNRMatrix[DIM * i + k][DIM * j + k] += msupg * WJ;
                jacobianNRMatrix[DIM * LNN + i][DIM * j + k] += mpspg * WJ;
            };
        };

        double vpspg = 0.;
        for (int k = 0; k < DIM; k++)
        {
            double vsupg = 0.;
            for (int l = 0; l < DIM; l++)
                vsupg += (una_[l] - umeshna_[l]) * dphi_dx[l][i];
            vsupg *= lambda_[k] * tSUPG_ * k1;

            rhsVector[DIM * i + k] += vsupg * WJ;

            vpspg += dphi_dx[k][i] * lambda_[k] * tPSPG_ * k1 / dens_;
        };

        rhsVector[DIM * LNN + i] += vpspg * WJ;
    };
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsSameMeshArlqStabLaplace(int &LNN, double &wna_, double &djac_, double &weight_,
                                                              double &tARLQ0_, double &tARLQ1_, double *phi_,
                                                              double **dphi_dx, double ***ddphi_dx,
                                                              double **arlequinStabD, double *arlequinStabVectorD,
                                                              double **arlequinStab1, double *arlequinStabVector1)
{

    double &k1 = parameters->getArlequinK1();

    // Lagrange Multiplieres derivatives
    // double **dlambda_dx;
    // dlambda_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     dlambda_dx[i] = new double[DIM];
    // getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

    // Get force derivative
    double dForce_[DIM];
    double x_[DIM], xPrev_[DIM];
    getInterpCoord(LNN, phi_, x_, xPrev_);
    dForce_[0] = -(std::exp(x_[0]) * x_[0] * (x_[0] * (x_[1] * x_[1] - x_[1] + 2.) + 3 * x_[1] * x_[1] - 3 * x_[1] - 2.));
    dForce_[1] = 0.;

    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            double LL = 0.;
            for (int k = 0; k < DIM; k++)
                LL += -dphi_dx[k][i] * dphi_dx[k][j] * (tARLQ0_ + tARLQ1_) * k1;
            arlequinStabD[i][j] += LL * WJ / wna_;
        };

        // rhs vector
        // double LL = 0.;
        // for (int l = 0; l < DIM; l++) LL += dphi_dx[l][i] * dlambda_dx[k][l] * (tARLQ0_ + tARLQ1_) * k1;
        // arlequinStabVectorD[i] += LL * WJ / wna_;

        double gradf = 0.0;
        for (int l = 0; l < DIM; l++)
            gradf += -dphi_dx[l][i] * dForce_[l] * tARLQ0_;
        arlequinStabVector1[i] += gradf * WJ;
    };

    // for (int i = 0; i < DIM; ++i)
    //     delete[] dlambda_dx[i];
    // delete[] dlambda_dx;
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsSameMeshArlqStab(int &LNN, double &wna_, double &djac_, double &weight_, double &tARLQ_,
                                                       double *phi_, double **dphi_dx, double ***ddphi_dx,
                                                       double **arlequinStabD, double *arlequinStabVectorD,
                                                       double **arlequinStab1, double *arlequinStabVector1)
{
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();
    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh Velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVel(LNN, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

    // Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duMesh_dx[i] = new double[DIM];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duMeshPrev_dx[i] = new double[DIM];

    double dumeshna_dx[DIM][DIM];
    getInterpMeshVelDer(LNN, dphi_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];

    // Acceleration derivatives
    double **daccel_dx;
    daccel_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        daccel_dx[i] = new double[DIM];

    double **daccelPrev_dx;
    daccelPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        daccelPrev_dx[i] = new double[DIM];

    double daccelm_dx[DIM][DIM];
    getInterpAccelDer(LNN, dphi_dx, daccel_dx, daccelPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            daccelm_dx[i][j] = alpha_m * daccel_dx[i][j] + (1. - alpha_m) * daccelPrev_dx[i][j];

    // Lagrange Multiplieres derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dlambda_dx[i] = new double[DIM];
    getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

    // Second pressure derivatives
    double **ddpress_dxdx;
    ddpress_dxdx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ddpress_dxdx[i] = new double[DIM];
    getInterpSecondPressDer(LNN, ddphi_dx, ddpress_dxdx);

    // Second velocity derivatives
    double ***ddu_dxdx;
    ddu_dxdx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddu_dxdx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddu_dxdx[i][j] = new double[DIM];
    };

    double ***dduPrev_dxdx;
    dduPrev_dxdx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        dduPrev_dxdx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            dduPrev_dxdx[i][j] = new double[DIM];
    };

    double dduna_dxdx[DIM][DIM][DIM];
    getInterpSecondVelDer(LNN, ddphi_dx, ddu_dxdx, dduPrev_dxdx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                dduna_dxdx[i][j][k] = alpha_f * ddu_dxdx[i][j][k] + (1. - alpha_f) * dduPrev_dxdx[i][j][k];

    double AGDT = alpha_f * gamma * dTime_;
    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {

            double mass = 0.;
            double convec = 0.;
            double LL = 0.;
            for (int k = 0; k < DIM; k++)
            {
                mass += -dphi_dx[k][i] * dphi_dx[k][j] * tARLQ_;
                LL += dphi_dx[k][i] * dphi_dx[k][j] * tARLQ_ * k1 / dens_;
                for (int l = 0; l < DIM; l++)
                    convec += -dphi_dx[k][i] * ((una_[l] - umeshna_[l]) * ddphi_dx[l][k][j] + (duna_dx[l][k] - dumeshna_dx[l][k]) * dphi_dx[l][j]) * tARLQ_;
            };

            for (int k = 0; k < DIM; k++)
            {
                double press = 0.;
                for (int l = 0; l < DIM; l++)
                    press += -dphi_dx[l][i] * ddphi_dx[l][k][j] * tARLQ_ / dens_;

                arlequinStab1[DIM * i + k][DIM * LNN + j] += press * WJ;
                arlequinStab1[DIM * i + k][DIM * j + k] += (mass * alpha_m + convec * AGDT) * WJ;
                arlequinStabD[DIM * i + k][DIM * j + k] += LL * WJ / wna_;
            };
        };

        // rhs vector
        for (int k = 0; k < DIM; k++)
        {

            double mass = 0.;
            double convec = 0.;
            double press = 0.;
            double LL = 0.;
            for (int l = 0; l < DIM; l++)
            {

                mass += dphi_dx[l][i] * daccelm_dx[k][l] * tARLQ_;
                press += dphi_dx[l][i] * ddpress_dxdx[l][k] * tARLQ_ / dens_;
                LL += -dphi_dx[l][i] * dlambda_dx[k][l] * tARLQ_ * k1 / dens_;

                for (int m = 0; m < DIM; m++)
                    convec += dphi_dx[l][i] * ((una_[m] - umeshna_[m]) * dduna_dxdx[k][m][l] + (duna_dx[m][l] - dumeshna_dx[m][l]) * duna_dx[k][m]) * tARLQ_;
            };

            arlequinStabVector1[DIM * i + k] += (mass + convec + press) * WJ;
            arlequinStabVectorD[DIM * i + k] += LL * WJ / wna_;
        };
    };

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] daccel_dx[i];
    delete[] daccel_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] daccelPrev_dx[i];
    delete[] daccelPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] dlambda_dx[i];
    delete[] dlambda_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ddpress_dxdx[i];
    delete[] ddpress_dxdx;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; j++)
            delete[] ddu_dxdx[i][j];
        delete[] ddu_dxdx[i];
    };
    delete[] ddu_dxdx;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; j++)
            delete[] dduPrev_dxdx[i][j];
        delete[] dduPrev_dxdx[i];
    };
    delete[] dduPrev_dxdx;
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsDifferentMesh(double &djac_, double &weight_, int &LNN, double *phi_,
                                                    double **dphi_dx, int &LNNC, double *phiC_, double **dphiC_dx,
                                                    double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)

{

    // double &k1 = parameters->getArlequinK1();
    // double &alpha_f = parameters->getAlphaF();

    // // Interpolated variables
    // // Velocity
    // double u_[DIM], uPrev_[DIM], una_[DIM];
    // getInterpVelCoarse(LNNC, phiC_, u_, uPrev_);
    // for (int i = 0; i < DIM; i++)
    // una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // double lambda_[DIM];
    // getInterpLambda(LNN, phi_, lambda_);

    // double WJ = weight_ * djac_;

    // for (int i = 0; i < LNN; i++)
    // {
    //     for (int j = 0; j < LNNC; j++)
    //     {
    //         // L2 operator
    //         double L2 = phi_[i] * phiC_[j] * k1;
    //         for (int k = 0; k < DIM; k++) lagrMultMatrix[DIM*i+k][DIM*j+k] += L2 * WJ;
    //     };

    //     // L2 operator - Velocity
    //     for (int k = 0; k < DIM; k++){
    //         double L2U = -phi_[i] * una_[k] * k1;
    //         rhsVector2[DIM*i+k] += (L2U) * WJ;
    //     };

    // };

    // for (int i = 0; i < LNNC; i++)
    // {
    //     for (int k = 0; k < DIM; k++){
    //         double L2L = -phiC_[i] * lambda_[k] * k1;
    //         rhsVector1[DIM*i+k] += (L2L) * WJ;
    //     };

    // };

    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVelCoarse(LNNC, phiC_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDerCoarse(LNNC, dphiC_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

    double lambda_[DIM];
    getInterpLambda(LNN, phi_, lambda_);

    // Lagrange Multipliers derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dlambda_dx[i] = new double[DIM];
    getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

    double WJ = weight_ * djac_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNNC; j++)
        {

            // L2 operator
            double L2 = phi_[i] * phiC_[j] * k1;

            for (int k = 0; k < DIM; k++)
            {

                lagrMultMatrix[DIM * i + k][DIM * j + k] += L2 * WJ;

                for (int l = 0; l < DIM; l++)
                {
                    double H1 = 0.5 * dphi_dx[l][i] * dphiC_dx[k][j] * k2;
                    if (l == k)
                        for (int m = 0; m < DIM; m++)
                            H1 += 0.5 * dphi_dx[m][i] * dphiC_dx[m][j] * k2;
                    lagrMultMatrix[DIM * i + k][DIM * j + l] += H1 * WJ;
                };
            };
        };

        // L2 operator - Velocity
        for (int k = 0; k < DIM; k++)
        {
            double L2U = -phi_[i] * una_[k] * k1;

            double H1U = 0.;
            for (int l = 0; l < DIM; l++)
                H1U += -0.5 * (dphi_dx[l][i] * duna_dx[l][k] + dphi_dx[l][i] * duna_dx[k][l]) * k2;

            rhsVector2[DIM * i + k] += (L2U + H1U) * WJ;
        };
    };

    for (int i = 0; i < LNNC; i++)
    {

        for (int k = 0; k < DIM; k++)
        {

            double L2L = -phiC_[i] * lambda_[k] * k1;

            double H1L = 0.;
            for (int l = 0; l < DIM; l++)
                H1L += -0.5 * (dphiC_dx[l][i] * dlambda_dx[l][k] + dphiC_dx[l][i] * dlambda_dx[k][l]) * k2;

            rhsVector1[DIM * i + k] += (L2L + H1L) * WJ;
        };
    };

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] dlambda_dx[i];
    delete[] dlambda_dx;
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsDifferentMeshLaplace(double &djac_, double &weight_, int &LNN, double *phi_,
                                                           int &LNNC, double *phiC_,
                                                           double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)

{

    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM];
    getInterpVelCoarse(LNNC, phiC_, u_, uPrev_);

    double lambda_[DIM];
    getInterpLambda(LNN, phi_, lambda_);

    double WJ = weight_ * djac_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNNC; j++)
        {

            // L2 operator
            // phiC_[j]
            double L2 = phi_[i] * phiC_[j] * k1;
            lagrMultMatrix[i][j] += L2 * WJ;
        };

        // L2 operator - Velocity
        double L2U = -phi_[i] * u_[0] * k1;
        rhsVector2[i] += (L2U)*WJ;
    };

    for (int i = 0; i < LNNC; i++)
    {
        double L2L = -phiC_[i] * lambda_[0] * k1;
        rhsVector1[i] += (L2L)*WJ;
    };
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG(double &djac_, double &weight_, double &tSUPG_, double &tPSPG_,
                                                                int &LNN, double *phi_, int &LNNC, double *phiC_, double **dphiC_dx,
                                                                double **jacobianNRMatrix, double *rhsVector)
{

    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();

    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVelCoarse(LNNC, phiC_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh Velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVelCoarse(LNNC, phiC_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    double lambda_[DIM];
    getInterpLambda(LNN, phi_, lambda_);

    double WJ = weight_ * djac_;

    for (int i = 0; i < LNNC; i++)
    {
        for (int j = 0; j < LNN; j++)
        {

            double msupg = 0.;
            for (int k = 0; k < DIM; k++)
                msupg += (una_[k] - umeshna_[k]) * dphiC_dx[k][i] * phi_[j] * tSUPG_ * k1;

            for (int k = 0; k < DIM; k++)
            {

                jacobianNRMatrix[DIM * i + k][DIM * j + k] += msupg * WJ;

                double mpspg = (dphiC_dx[k][i] * phi_[j]) * tPSPG_ * k1 / dens_;
                jacobianNRMatrix[DIM * LNNC + i][DIM * j + k] += mpspg * WJ;
            };
        };

        double vpspg = 0.;
        for (int k = 0; k < DIM; k++)
        {

            double vsupg = 0.;
            for (int m = 0; m < DIM; m++)
            {
                vsupg += -(una_[m] - umeshna_[m]) * dphiC_dx[m][i];
            };
            vsupg *= lambda_[k] * tSUPG_ * k1;
            rhsVector[DIM * i + k] += vsupg * WJ;

            vpspg += -dphiC_dx[k][i] * lambda_[k] * tPSPG_ * k1 / dens_;
        };

        rhsVector[DIM * LNNC + i] += vpspg * WJ;
    };
};

template <int DIM>
void Element<DIM>::getMatrixAndVectorsDifferentMeshArlqStabLaplace(int &LNN, int &LNNC, double &wna_, double &djac_,
                                                                   double &weight_, double &tARLQ0_, double &tARLQ1_,
                                                                   double *phi_, double **dphi_dx, double **dphiC_dx, double ***ddphiC_dx,
                                                                   double **arlequinStabD, double *arlequinStabVectorD,
                                                                   double **arlequinStab0, double *arlequinStabVector0)
{

    double &k1 = parameters->getArlequinK1();

    // Lagrange Multiplieres derivatives
    // double **dlambda_dx;
    // dlambda_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     dlambda_dx[i] = new double[DIM];
    // getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

    // Coarse velocity derivatives
    // double **duC_dx;
    // duC_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     duC_dx[i] = new double[DIM];

    // double **duCPrev_dx;
    // duCPrev_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     duCPrev_dx[i] = new double[DIM];

    // getInterpVelDerCoarse(LNNC, dphiC_dx,duC_dx,duCPrev_dx);

    // Coarse second velocity derivatives
    // double ***dduC_dxdx;
    // dduC_dxdx = new double **[DIM];
    // for (int i = 0; i < DIM; ++i)
    // {
    //     dduC_dxdx[i] = new double *[DIM];
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         dduC_dxdx[i][j] = new double[DIM];
    //     };
    // };

    // double ***dduCPrev_dxdx;
    // dduCPrev_dxdx = new double **[DIM];
    // for (int i = 0; i < DIM; ++i)
    // {
    //     dduCPrev_dxdx[i] = new double *[DIM];
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         dduCPrev_dxdx[i][j] = new double[DIM];
    //     };
    // };
    // getInterpSecondVelDerCoarse(LNNC, ddphiC_dx, dduC_dxdx, dduCPrev_dxdx);

    // Get force derivative
    double dForce_[DIM];
    double x_[DIM], xPrev_[DIM];
    getInterpCoord(LNN, phi_, x_, xPrev_);
    dForce_[0] = -(std::exp(x_[0]) * x_[0] * (x_[0] * (x_[1] * x_[1] - x_[1] + 2.) + 3 * x_[1] * x_[1] - 3 * x_[1] - 2.));
    dForce_[1] = 0.;

    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNN; j++)
        {
            double LL = 0.;
            for (int k = 0; k < DIM; k++)
                LL += -dphi_dx[k][i] * dphi_dx[k][j] * (tARLQ0_ + tARLQ1_) * k1;

            arlequinStabD[i][j] += LL * WJ / wna_;
        };

        // rhs vector
        // double LL = 0.;
        // for (int l = 0; l < DIM; l++) LL += dphi_dx[l][i] * dlambda_dx[k][l] * (tARLQ0_ + tARLQ1_) * k1;
        // arlequinStabVectorD[i] += LL * WJ / wna_;

        double gradf = 0.0;
        for (int l = 0; l < DIM; l++)
            gradf += dphi_dx[l][i] * dForce_[l] * tARLQ1_;
        arlequinStabVector0[i] += gradf * WJ;
    };

    // for (int i = 0; i < DIM; ++i) delete[] dlambda_dx[i];
    // delete[] dlambda_dx;

    // for (int i = 0; i < DIM; ++i)
    //     delete[] duC_dx[i];
    // delete[] duC_dx;

    // for (int i = 0; i < DIM; ++i)
    //     delete[] duCPrev_dx[i];
    // delete[] duCPrev_dx;

    // for (int i = 0; i < DIM; ++i)
    // {
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         delete[] dduC_dxdx[i][j];
    //     };
    //     delete[] dduC_dxdx[i];
    // };
    // delete[] dduC_dxdx;

    // for (int i = 0; i < DIM; ++i)
    // {
    //     for (int j = 0; j < DIM; j++)
    //     {
    //         delete[] dduCPrev_dxdx[i][j];
    //     };
    //     delete[] dduCPrev_dxdx[i];
    // };
    // delete[] dduCPrev_dxdx;
}

template <int DIM>
void Element<DIM>::getMatrixAndVectorsDifferentMeshArlqStab(int &LNN, int &LNNC, double &wna_, double &djac_, double &weight_, double &tARLQ_,
                                                            double **dphi_dx, double *phiC_, double **dphiC_dx, double ***ddphiC_dx,
                                                            double **arlequinStabD, double *arlequinStabVectorD,
                                                            double **arlequinStab0, double *arlequinStabVector0)
{
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();
    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVelCoarse(LNNC, phiC_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh Velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVelCoarse(LNNC, phiC_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDerCoarse(LNNC, dphiC_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

    // Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duMesh_dx[i] = new double[DIM];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duMeshPrev_dx[i] = new double[DIM];

    double dumeshna_dx[DIM][DIM];
    getInterpMeshVelDerCoarse(LNNC, dphiC_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] +
                                (1. - alpha_f) * duMeshPrev_dx[i][j];

    // Acceleration derivatives
    double **daccel_dx;
    daccel_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        daccel_dx[i] = new double[DIM];

    double **daccelPrev_dx;
    daccelPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        daccelPrev_dx[i] = new double[DIM];

    double daccelm_dx[DIM][DIM];
    getInterpAccelDerCoarse(LNNC, dphiC_dx, daccel_dx, daccelPrev_dx);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            daccelm_dx[i][j] = alpha_m * daccel_dx[i][j] +
                               (1. - alpha_m) * daccelPrev_dx[i][j];

    // Second pressure derivatives
    double **ddpress_dxdx;
    ddpress_dxdx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ddpress_dxdx[i] = new double[DIM];
    getInterpSecondPressDerCoarse(LNNC, ddphiC_dx, ddpress_dxdx);

    // Second velocity derivatives
    double ***ddu_dxdx;
    ddu_dxdx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddu_dxdx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
        {
            ddu_dxdx[i][j] = new double[DIM];
        };
    };

    double ***dduPrev_dxdx;
    dduPrev_dxdx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        dduPrev_dxdx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
        {
            dduPrev_dxdx[i][j] = new double[DIM];
        };
    };

    double dduna_dxdx[DIM][DIM][DIM];
    getInterpSecondVelDerCoarse(LNNC, ddphiC_dx, ddu_dxdx, dduPrev_dxdx);

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                dduna_dxdx[i][j][k] = alpha_f * ddu_dxdx[i][j][k] +
                                      (1. - alpha_f) * dduPrev_dxdx[i][j][k];

    // Lagrange Multiplieres derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dlambda_dx[i] = new double[DIM];
    getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

    double AGDT = alpha_f * gamma * dTime_;
    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < LNN; i++)
    {
        for (int j = 0; j < LNNC; j++)
        {

            double mass = 0.;
            double convec = 0.;
            for (int k = 0; k < DIM; k++)
            {
                mass += dphi_dx[k][i] * dphiC_dx[k][j] * tARLQ_;
                for (int l = 0; l < DIM; l++)
                    convec += dphi_dx[k][i] * ((una_[l] - umeshna_[l]) * ddphiC_dx[l][k][j] + (duna_dx[l][k] - dumeshna_dx[l][k]) * dphiC_dx[l][j]) * tARLQ_;
            };

            for (int k = 0; k < DIM; k++)
            {

                double press = 0.;
                for (int l = 0; l < DIM; l++)
                    press += dphi_dx[l][i] * ddphiC_dx[l][k][j] * tARLQ_ / dens_;

                arlequinStab0[DIM * i + k][DIM * LNNC + j] += press * WJ;
                arlequinStab0[DIM * i + k][DIM * j + k] += (0.*mass * alpha_m + convec * AGDT) * WJ;
            };
        };

        // rhs vector
        for (int k = 0; k < DIM; k++)
        {

            double mass = 0.;
            double convec = 0.;
            double press = 0.;
            double LL = 0.;
            for (int l = 0; l < DIM; l++)
            {

                mass += -dphi_dx[l][i] * daccelm_dx[k][l] * tARLQ_;
                press += -dphi_dx[l][i] * ddpress_dxdx[l][k] * tARLQ_ / dens_;
                LL += -dphi_dx[l][i] * dlambda_dx[k][l] * tARLQ_ * k1 / dens_;

                for (int m = 0; m < DIM; m++)
                    convec += -dphi_dx[l][i] * ((una_[m] - umeshna_[m]) * dduna_dxdx[k][m][l] + (duna_dx[m][l] - dumeshna_dx[m][l]) * duna_dx[k][m]) * tARLQ_;
            };

            arlequinStabVector0[DIM * i + k] += (0.*mass + convec + press) * WJ;
            arlequinStabVectorD[DIM * i + k] += LL * WJ / wna_;
        };

        for (int j = 0; j < LNN; j++)
        {
            double LL = 0.;
            for (int k = 0; k < DIM; k++)
                LL += dphi_dx[k][i] * dphi_dx[k][j] * tARLQ_ * k1 / dens_;

            for (int k = 0; k < DIM; k++)
                arlequinStabD[DIM * i + k][DIM * j + k] += LL * WJ / wna_;
        };
    };

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] daccel_dx[i];
    delete[] daccel_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] daccelPrev_dx[i];
    delete[] daccelPrev_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] dlambda_dx[i];
    delete[] dlambda_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ddpress_dxdx[i];
    delete[] ddpress_dxdx;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; j++)
        {
            delete[] ddu_dxdx[i][j];
        };
        delete[] ddu_dxdx[i];
    };
    delete[] ddu_dxdx;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; j++)
        {
            delete[] dduPrev_dxdx[i][j];
        };
        delete[] dduPrev_dxdx[i];
    };
    delete[] dduPrev_dxdx;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template <int DIM>
void Element<DIM>::getResidualVector(int &LNN, double &wna_, double &djac_, double &weight_,
                                     double &tSUPG_, double &tPSPG_, double &tLSIC_,
                                     double *phi_, double **dphi_dx, double *rhsVector)
{

    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();

    double ff[DIM];
    for (int i = 0; i < DIM; i++)
        ff[i] = parameters->getFieldForce(i);

    // Interpolated variables
    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVel(LNN, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Acceleration
    double accel_[DIM], accelPrev_[DIM], accelm_[DIM];
    getInterpAccel(LNN, phi_, accel_, accelPrev_);
    for (int i = 0; i < DIM; i++)
        accelm_[i] = alpha_m * accel_[i] + (1. - alpha_m) * accelPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    // Pressure
    double press_;
    getInterpPress(LNN, phi_, press_);

    // Pressure derivatives
    double dpress_dx[DIM];
    getInterpPressDer(LNN, dphi_dx, dpress_dx);

    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < LNN; i++)
    {
        double aux0 = 0.;
        for (int k = 0; k < DIM; k++)
            aux0 += duna_dx[k][k];

        double wSUPGi = 0.;
        for (int k = 0; k < DIM; k++)
            wSUPGi += (una_[k] - umeshna_[k]) * dphi_dx[k][i];

        for (int j = 0; j < DIM; j++)
        {

            double mm = phi_[i] * accelm_[j] * dens_ + wSUPGi * accelm_[j] * tSUPG_ * dens_;

            double kls = dphi_dx[j][i] * aux0 * tLSIC_ * dens_;

            double kk = 0.;
            for (int k = 0; k < DIM; k++) kk += dphi_dx[k][i] * (duna_dx[j][k] + duna_dx[k][j]) * visc_;


            double aux1 = 0.;
            for (int k = 0; k < DIM; k++)
                aux1 += duna_dx[j][k] * (una_[k] - umeshna_[k]);
            double cc = aux1 * phi_[i] * dens_ +
                        wSUPGi * aux1 * tSUPG_ * dens_;

            double pp = -(dphi_dx[j][i] * press_) +
                        wSUPGi * dpress_dx[j] * tSUPG_;

            double ffv = phi_[i] * dens_ * ff[j] +
                         tSUPG_ * dens_ * wSUPGi * ff[j];

            rhsVector[DIM * i + j] += (-mm + ffv - kk - pp - cc - kls) * WJ;
        };

        double aux2 = 0.;
        double aux3 = 0.;
        double aux4 = 0.;
        double aux5 = 0.;
        for (int k = 0; k < DIM; k++)
        {
            aux2 += dphi_dx[k][i] * dpress_dx[k];
            aux3 += dphi_dx[k][i] * accelm_[k];
            for (int l = 0; l < DIM; l++)
                aux4 += dphi_dx[k][i] * ((una_[l] - umeshna_[l]) * duna_dx[k][l]);
            aux5 += dphi_dx[k][i] * ff[k];
        };

        double qq = aux0 * phi_[i] +
                    aux2 * tPSPG_ / dens_ +
                    aux3 * tPSPG_ +
                    aux4 * tPSPG_;

        double ffp = aux5 * tPSPG_;

        rhsVector[DIM * LNN + i] += (ffp - qq) * WJ;
    };

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    return;
};

template <int DIM>
void Element<DIM>::getResidualVectorTNS(int &iTime, int &LNN, double &wna_, double &djac_, double &weight_,
                                        double &tSUPG_, double &tPSPG_, double &tLSIC_,
                                        double *phi_, double **dphi_dx, double *rhsVector)
{

    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &dTime = parameters->getTimeStep();

    double ff[DIM] = {};
    // getInterpForce(LNN, phi_,ff);
    // for (int i = 0; i < DIM; i++)
    //     ff[i] = parameters->getFieldForce(i);

    // Interpolated variables
    // coord
    double coord_[DIM], coordPrev_[DIM];
    getInterpCoord(LNN, phi_, coord_, coordPrev_);

    // Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    getInterpVel(LNN, phi_, u_, uPrev_);
    for (int i = 0; i < DIM; i++)
        una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    // Mesh velocity
    double uMesh_[DIM], uMeshPrev_[DIM], umeshna_[DIM];
    getInterpMeshVel(LNN, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < DIM; i++)
        umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    // Acceleration
    double accel_[DIM], accelPrev_[DIM], accelm_[DIM];
    getInterpAccel(LNN, phi_, accel_, accelPrev_);
    for (int i = 0; i < DIM; i++)
        accelm_[i] = alpha_m * accel_[i] + (1. - alpha_m) * accelPrev_[i];

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    double duna_dx[DIM][DIM];
    getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

    // Pressure
    double press_;
    getInterpPress(LNN, phi_, press_);

    // Pressure derivatives
    double dpress_dx[DIM];
    getInterpPressDer(LNN, dphi_dx, dpress_dx);

    double WJ = weight_ * djac_ * wna_;

    // Navier Stokes
    for (int i = 0; i < LNN; i++)
    {
        double aux0 = 0.;
        for (int k = 0; k < DIM; k++)
            aux0 += duna_dx[k][k];

        double wSUPGi = 0.;
        for (int k = 0; k < DIM; k++)
            wSUPGi += (una_[k] - umeshna_[k]) * dphi_dx[k][i];

        for (int j = 0; j < DIM; j++)
        {

            double mm = phi_[i] * accelm_[j] * dens_ + wSUPGi * accelm_[j] * tSUPG_ * dens_;

            double kls = dphi_dx[j][i] * aux0 * tLSIC_ * dens_;

            double kk = 0.;
            for (int k = 0; k < DIM; k++)
                kk += dphi_dx[k][i] * (duna_dx[j][k] + duna_dx[k][j]) * visc_;


            double aux1 = 0.;
            for (int k = 0; k < DIM; k++)
                aux1 += duna_dx[j][k] * (una_[k] - umeshna_[k]);
            double cc = aux1 * phi_[i] * dens_ +
                        wSUPGi * aux1 * tSUPG_ * dens_;

            double pp = -(dphi_dx[j][i] * press_) +
                        wSUPGi * dpress_dx[j] * tSUPG_;

            double ffv = phi_[i] * dens_ * ff[j] +
                         tSUPG_ * dens_ * wSUPGi * ff[j];

            rhsVector[DIM * i + j] += (-mm + ffv - kk - pp - cc - kls) * WJ;
        };

        double aux2 = 0.;
        double aux3 = 0.;
        double aux4 = 0.;
        double aux5 = 0.;
        for (int k = 0; k < DIM; k++)
        {
            aux2 += dphi_dx[k][i] * dpress_dx[k];
            aux3 += dphi_dx[k][i] * accelm_[k];
            for (int l = 0; l < DIM; l++)
                aux4 += dphi_dx[k][i] * ((una_[l] - umeshna_[l]) * duna_dx[k][l]);
            aux5 += dphi_dx[k][i] * ff[k];
        };

        double qq = aux0 * phi_[i] +
                    aux2 * tPSPG_ / dens_ +
                    aux3 * tPSPG_ +
                    aux4 * tPSPG_;

        double ffp = aux5 * tPSPG_;
        
        rhsVector[DIM * LNN + i] += (ffp - qq) * WJ;
    };

    // Stokes
    // ff[0] = (1.-2.*pi) * (8.*pi*pi*visc_/dens_-1.) * cos(2*pi*coord_[0]) * sin(2*pi*coord_[1]) * exp(-iTime*dTime);
    // ff[1] = (1.+2.*pi) * (1.-8.*pi*pi*visc_/dens_) * cos(2*pi*coord_[1]) * sin(2*pi*coord_[0]) * exp(-iTime*dTime);

    //Stokes Problem
    // double aux0 = 0.;
    // for (int k = 0; k < DIM; k++) aux0 += duna_dx[k][k];

    // for (int i = 0; i < LNN; i++)
    // {
    //     for (int j = 0; j < DIM; j++)
    //     {

    //         double mm = phi_[i] * accelm_[j] * dens_;

    //         double kk = 0.;
    //         for (int k = 0; k < DIM; k++) kk += dphi_dx[k][i] * (duna_dx[j][k] + duna_dx[k][j]) * visc_;

    //         double pp = -(dphi_dx[j][i] * press_);

    //         double ffv = phi_[i] * dens_ * ff[j];

    //         rhsVector[DIM * i + j] += (-mm - kk - pp + ffv) * WJ;
    //     };

    //     double aux3 = 0.;
    //     double aux2 = 0.;
    //     double aux5 = 0.;
    //     for (int k = 0; k < DIM; k++){
    //         aux2 += dphi_dx[k][i] * dpress_dx[k];
    //         aux3 += dphi_dx[k][i] * accelm_[k];
    //         aux5 += dphi_dx[k][i] * ff[k];
    //     };

    //     double qq = aux0 * phi_[i]
    //                 + aux2 * tPSPG_ / dens_
    //                 + aux3 * tPSPG_;

    //     double ffp = aux5 * tPSPG_;

    //     rhsVector[DIM * LNN + i] += (ffp - qq) * WJ;

    // };

    // Laplace Problem

    // // Velocity Derivatives
    // double **du_dx;
    // du_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     du_dx[i] = new double[DIM];

    // double **duPrev_dx;
    // duPrev_dx = new double *[DIM];
    // for (int i = 0; i < DIM; ++i)
    //     duPrev_dx[i] = new double[DIM];

    // getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);

    // double WJ = wna_ * djac_ * weight_;

    // for (int i = 0; i < LNN; i++){
    //    for (int j=0; j < DIM; j++){
    //         double Lap = 0.;
    //         for (int k = 0; k < DIM; k++) Lap += dphi_dx[k][i] * du_dx[j][k];
    //         rhsVector[DIM*i+j] -= Lap * WJ;
    //    };
    // };

    // deallocating memory
    for (int i = 0; i < DIM; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;
};

template <int DIM>
void Element<DIM>::getResidualVectorLaplace(int &LNN, double &wna_, double &djac_,
                                            double &weight_, double *phi_,
                                            double **dphi_dx, double *rhsVector)
{

    double ff = parameters->getFieldForce(0);

    // Field force
    double x_[DIM], xPrev_[DIM];
    getInterpCoord(LNN, phi_, x_, xPrev_);
    double force_ = -(std::exp(x_[0]) * x_[0] * (-2. - 3. * x_[1] + 3. * x_[1] * x_[1] + x_[0] * (2. - x_[1] + x_[1] * x_[1])));

    // Velocity Derivatives
    double **du_dx;
    du_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        du_dx[i] = new double[DIM];

    double **duPrev_dx;
    duPrev_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        duPrev_dx[i] = new double[DIM];

    getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);

    double WJ = wna_ * djac_ * weight_;

    for (int i = 0; i < LNN; i++)
    {

        double ffv = phi_[i] * (ff + force_);

        double K = 0;
        for (int k = 0; k < DIM; k++)
            K += dphi_dx[k][i] * du_dx[0][k];

        rhsVector[i] += (ffv - K) * WJ;
    };
};

//------------------------------------------------------------------------------
//----------------------ARLEQUIN ELEMENT LOCAL MATRIX/VECTOR--------------------
//------------------------------------------------------------------------------
template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMesh_FEM(int &ip, double **lagrMultMatrix,
                                                      double *rhsVector1, double *rhsVector2)

{

    int LNN = 4 * DIM - 2;
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    double djac_, weight_;

    // data for computation of FEM basis functions
    double xsi[DIM], phi_[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    //int index = ip;

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

    // Defines the integration points adimentional coordinates
    for (int i = 0; i < DIM; i++)
        xsi[i] = nQuad.PointListFem(index, i);

    // Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    // Computes the velocity shape functions
    shapeQuad.evaluateFem(xsi, phi_);

    // Computes the jacobian matrix
    getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

    // Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

    // Computes matrixes and vectors
    getMatrixAndVectorsSameMesh(LNN, djac_, weight_, phi_, dphi_dx, lagrMultMatrix,
                                rhsVector1, rhsVector2);

    index++;
    }

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMesh_ISO(int &ip, double **lagrMultMatrix,
                                                      double *rhsVector1, double *rhsVector2)

{

    int LNN = 18 * DIM - 27;
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    double djac_, weight_;

    // data for IGA functions
    double wpc[LNN], phi_[LNN], xsi[DIM];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    double **quadJacMat;
    quadJacMat = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        quadJacMat[i] = new double[DIM];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListIso(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        // Computes the jacobian matrix
        getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

        // Computes matrixes and vectors
        getMatrixAndVectorsSameMesh(LNN, djac_, weight_, phi_, dphi_dx, lagrMultMatrix,
                                    rhsVector1, rhsVector2);
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] quadJacMat[i];
    delete[] quadJacMat;

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMeshLaplace_FEM(int &index, double **lagrMultMatrix,
                                                             double *rhsVector1, double *rhsVector2)

{
    int LNN = 4 * DIM - 2;
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    double djac_, weight_;

    // data for computation of FEM basis functions
    double xsi[DIM], phi_[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    // Defines the integration points adimentional coordinates
    for (int i = 0; i < DIM; i++)
        xsi[i] = nQuad.PointListFem(index, i);

    // Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    // Computes the velocity shape functions
    shapeQuad.evaluateFem(xsi, phi_);

    // Computes the jacobian matrix
    getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

    // Computes matrixes and vectors
    getMatrixAndVectorsSameMeshLaplace(LNN, djac_, weight_, phi_, lagrMultMatrix,
                                       rhsVector1, rhsVector2);

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMeshLaplace_ISO(double **lagrMultMatrix,
                                                             double *rhsVector1, double *rhsVector2)

{
    int LNN = 18 * DIM - 27;
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    double djac_, weight_;

    // data for IGA functions
    double wpc[LNN], phi_[LNN], xsi[DIM];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    double **quadJacMat;
    quadJacMat = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        quadJacMat[i] = new double[DIM];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListIso(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        // Computes the jacobian matrix
        getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

        // Computes matrixes and vectors
        getMatrixAndVectorsSameMeshLaplace(LNN, djac_, weight_, phi_, lagrMultMatrix,
                                           rhsVector1, rhsVector2);
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] quadJacMat[i];
    delete[] quadJacMat;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMesh_FEM_TESTE(double **lagrMultMatrix,
                                                            double *rhsVector1, double *rhsVector2)

{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    double djac_, weight_;

    // data for computation of FEM basis functions
    double xsi[3], phi_[10];

    // Fluid Data
    double &alpha_f = parameters->getAlphaF();
    double &dens_ = parameters->getDensity();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        double **dphi_dx;
        dphi_dx = new double *[3];
        for (int i = 0; i < 3; ++i)
            dphi_dx[i] = new double[10];

        double **ainv_;
        ainv_ = new double *[3];
        for (int i = 0; i < 3; ++i)
            ainv_[i] = new double[3];

        double **Jac;
        Jac = new double *[3];
        for (int i = 0; i < 3; ++i)
            Jac[i] = new double[3];

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < 3; i++)
            xsi[i] = nQuad.PointListFem(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListFem(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi, phi_);

        // Computes the jacobian matrix
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

        // Interpolated variables
        // Velocity
        int LNN = 10;
        double u_[3], uPrev_[3], una_[3];
        getInterpVel(LNN, phi_, u_, uPrev_);
        for (int i = 0; i < 3; i++)
            una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

        // Velocity Derivatives
        double **du_dx;
        du_dx = new double *[3];
        for (int i = 0; i < 3; ++i)
            du_dx[i] = new double[3];

        double **duPrev_dx;
        duPrev_dx = new double *[3];
        for (int i = 0; i < 3; ++i)
            duPrev_dx[i] = new double[3];

        double duna_dx[3][3];
        getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

        // Lagrange Multiplieres
        double lambda_[3];
        getInterpLambda(LNN, phi_, lambda_);

        // Lagrange Multiplieres derivatives
        double **dlambda_dx;
        dlambda_dx = new double *[3];
        for (int i = 0; i < 3; ++i)
            dlambda_dx[i] = new double[3];
        getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

        double WJ = weight_ * djac_;

        for (int i = 0; i < 10; i++)
        {
            for (int j = 0; j < 10; j++)
            {
                // L2 operator
                double L2 = -phi_[i] * phi_[j] * k1;

                lagrMultMatrix[3 * i][3 * j] += L2 * WJ;
                lagrMultMatrix[3 * i + 1][3 * j + 1] += L2 * WJ;
                lagrMultMatrix[3 * i + 2][3 * j + 2] += L2 * WJ;

                // H1 operator
                double H1xx = -(dphi_dx[0][i] * dphi_dx[0][j] +
                                0.5 * dphi_dx[1][i] * dphi_dx[1][j] +
                                0.5 * dphi_dx[2][i] * dphi_dx[2][j]) *
                              k2;
                double H1xy = -0.5 * dphi_dx[1][i] * dphi_dx[0][j] * k2;
                double H1xz = -0.5 * dphi_dx[2][i] * dphi_dx[0][j] * k2;
                double H1yx = -0.5 * dphi_dx[0][i] * dphi_dx[1][j] * k2;
                double H1yy = -(dphi_dx[1][i] * dphi_dx[1][j] +
                                0.5 * dphi_dx[0][i] * dphi_dx[0][j] +
                                0.5 * dphi_dx[2][i] * dphi_dx[2][j]) *
                              k2;
                double H1yz = -0.5 * dphi_dx[2][i] * dphi_dx[1][j] * k2;
                double H1zx = -0.5 * dphi_dx[0][i] * dphi_dx[2][j] * k2;
                double H1zy = -0.5 * dphi_dx[1][i] * dphi_dx[2][j] * k2;
                double H1zz = -(dphi_dx[2][i] * dphi_dx[2][j] +
                                0.5 * dphi_dx[0][i] * dphi_dx[0][j] +
                                0.5 * dphi_dx[1][i] * dphi_dx[1][j]) *
                              k2;

                lagrMultMatrix[3 * i][3 * j] += H1xx * WJ;
                lagrMultMatrix[3 * i][3 * j + 1] += H1xy * WJ;
                lagrMultMatrix[3 * i][3 * j + 2] += H1xz * WJ;
                lagrMultMatrix[3 * i + 1][3 * j] += H1yx * WJ;
                lagrMultMatrix[3 * i + 1][3 * j + 1] += H1yy * WJ;
                lagrMultMatrix[3 * i + 1][3 * j + 2] += H1yz * WJ;
                lagrMultMatrix[3 * i + 2][3 * j] += H1zx * WJ;
                lagrMultMatrix[3 * i + 2][3 * j + 1] += H1zy * WJ;
                lagrMultMatrix[3 * i + 2][3 * j + 2] += H1zy * WJ;
            };

            double H1Lx = dphi_dx[0][i] * dlambda_dx[0][0] +
                          0.5 * dphi_dx[1][i] * dlambda_dx[0][1] +
                          0.5 * dphi_dx[2][i] * dlambda_dx[0][2] +
                          0.5 * dphi_dx[1][i] * dlambda_dx[1][0] +
                          0.5 * dphi_dx[2][i] * dlambda_dx[2][0];
            double H1Ly = 0.5 * dphi_dx[0][i] * dlambda_dx[0][1] +
                          dphi_dx[1][i] * dlambda_dx[1][1] +
                          0.5 * dphi_dx[0][i] * dlambda_dx[1][0] +
                          0.5 * dphi_dx[2][i] * dlambda_dx[1][2] +
                          0.5 * dphi_dx[2][i] * dlambda_dx[2][1];
            double H1Lz = 0.5 * dphi_dx[0][i] * dlambda_dx[0][2] +
                          0.5 * dphi_dx[1][i] * dlambda_dx[1][2] +
                          0.5 * dphi_dx[0][i] * dlambda_dx[2][0] +
                          0.5 * dphi_dx[1][i] * dlambda_dx[2][1] +
                          dphi_dx[2][i] * dlambda_dx[2][2];

            rhsVector1[3 * i] += (phi_[i] * lambda_[0] * k1 + H1Lx * k2) * WJ;
            rhsVector1[3 * i + 1] += (phi_[i] * lambda_[1] * k1 + H1Ly * k2) * WJ;
            rhsVector1[3 * i + 2] += (phi_[i] * lambda_[2] * k1 + H1Lz * k2) * WJ;

            double H1Ux = dphi_dx[0][i] * duna_dx[0][0] +
                          0.5 * dphi_dx[1][i] * duna_dx[0][1] +
                          0.5 * dphi_dx[2][i] * duna_dx[0][2] +
                          0.5 * dphi_dx[1][i] * duna_dx[1][0] +
                          0.5 * dphi_dx[2][i] * duna_dx[2][0];
            double H1Uy = 0.5 * dphi_dx[0][i] * duna_dx[0][1] +
                          dphi_dx[1][i] * duna_dx[1][1] +
                          0.5 * dphi_dx[0][i] * duna_dx[1][0] +
                          0.5 * dphi_dx[2][i] * duna_dx[1][2] +
                          0.5 * dphi_dx[2][i] * duna_dx[2][1];
            double H1Uz = 0.5 * dphi_dx[0][i] * duna_dx[0][2] +
                          0.5 * dphi_dx[1][i] * duna_dx[1][2] +
                          0.5 * dphi_dx[0][i] * duna_dx[2][0] +
                          0.5 * dphi_dx[1][i] * duna_dx[2][1] +
                          dphi_dx[2][i] * duna_dx[2][2];

            rhsVector2[3 * i] += (phi_[i] * una_[0] * k1 + H1Ux * k2) * WJ;
            rhsVector2[3 * i + 1] += (phi_[i] * una_[1] * k1 + H1Uy * k2) * WJ;
            rhsVector2[3 * i + 2] += (phi_[i] * una_[2] * k1 + H1Uz * k2) * WJ;
        };

        // deallocating memory
        for (int i = 0; i < 3; ++i)
            delete[] du_dx[i];
        delete[] du_dx;

        for (int i = 0; i < 3; ++i)
            delete[] duPrev_dx[i];
        delete[] duPrev_dx;

        for (int i = 0; i < 3; ++i)
            delete[] dlambda_dx[i];
        delete[] dlambda_dx;

        for (int i = 0; i < 3; ++i)
            delete[] dphi_dx[i];
        delete[] dphi_dx;

        for (int i = 0; i < 3; ++i)
            delete[] ainv_[i];
        delete[] ainv_;

        for (int i = 0; i < 3; ++i)
            delete[] Jac[i];
        delete[] Jac;

        index++;
    };

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMesh_tSUPG_tPSPG_FEM(int &index, double **jacobianNRMatrix, double *rhsVector)
{
    int LNN = 4 * DIM - 2;

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // data for computation of IGA basis functions
    double phi_[LNN], xsi[DIM];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    double tSUPG_, tPSPG_, tLSIC_;
    double djac_, weight_;

    // Defines the integration points adimentional coordinates
    for (int i = 0; i < DIM; i++)
        xsi[i] = nQuad.PointListFem(index, i);

    // Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    // Computes the velocity shape functions
    shapeQuad.evaluateFem(xsi, phi_);

    // Computes the jacobian matrix
    getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

    // Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

    // Computes stabilization parameters
    getNewParameterSUPG_FEM(tSUPG_, tPSPG_, tLSIC_, Jac, phi_, dphi_dx);

    // Computes matrixes and vectors
    getMatrixAndVectorsSameMesh_tSUPG_tPSPG(LNN, djac_, weight_, tSUPG_, tPSPG_, phi_, dphi_dx,
                                            jacobianNRMatrix, rhsVector);

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    return;
}

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMeshArlqStab_FEM(double **arlequinStabD, double *arlequinStabVectorD)
{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // data for computation of IGA basis functions
    double xsi[2], phi_[6];
    int LNN = 4 * DIM - 2;

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i)
        dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i)
        ainv_[i] = new double[2];

    double **Jac;
    Jac = new double *[2];
    for (int i = 0; i < 2; ++i)
        Jac[i] = new double[2];

    double tSUPG_, tPSPG_, tLSIC_;

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        // Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListFem(index, 0);
        xsi[1] = nQuad.PointListFem(index, 1);
        // Returns the quadrature integration weight
        double weight_ = nQuad.WeightListFem(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi, phi_);

        // Computes the jacobian matrix
        double djac_ = 0.0;
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

        // Computes stabilization parameters
        getNewParameterSUPG_FEM(tSUPG_, tPSPG_, tLSIC_, Jac, phi_, dphi_dx);
        double tARLQ_ = tSUPG_ * 1.e-2;

        // Lagrange Multiplieres derivatives
        double **dlambda_dx;
        dlambda_dx = new double *[2];
        for (int i = 0; i < 2; ++i)
            dlambda_dx[i] = new double[2];
        getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

        double WJ = weight_ * djac_;
        double k1 = 1.;
        double dens_ = 1.;

        for (int i = 0; i < LNN; i++)
        {
            for (int j = 0; j < LNN; j++)
            {

                // Diagonal matrix
                double LL = -(dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ * k1 / dens_;

                arlequinStabD[2 * i][2 * j] += LL * WJ;
                arlequinStabD[2 * i + 1][2 * j + 1] += LL * WJ;
            };

            // Diagonal matrix
            double LLx = (dphi_dx[0][i] * dlambda_dx[0][0] + dphi_dx[1][i] * dlambda_dx[0][1]) * tARLQ_ * k1 / dens_;
            double LLy = (dphi_dx[0][i] * dlambda_dx[1][0] + dphi_dx[1][i] * dlambda_dx[1][1]) * tARLQ_ * k1 / dens_;

            arlequinStabVectorD[2 * i] += LLx * WJ;
            arlequinStabVectorD[2 * i + 1] += LLy * WJ;
        };

        // deallocating memory
        for (int i = 0; i < 2; ++i)
            delete[] dlambda_dx[i];
        delete[] dlambda_dx;

        index++;
    };

    for (int i = 0; i < 2; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;
    for (int i = 0; i < 2; ++i)
        delete[] ainv_[i];
    delete[] ainv_;
    for (int i = 0; i < 2; ++i)
        delete[] Jac[i];
    delete[] Jac;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMeshArlqStab_FEM_ISO(int &index, int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                  std::vector<IParameters *> &iparamC, double **arlequinStabD,
                                                                  double *arlequinStabVectorD, double **arlequinStab1, double *arlequinStabVector1)
{
    int LNN = 4 * DIM - 2;
    int LNNC = 18 * DIM - 27;

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // data for IGA coarse mesh computation
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[LNNC], phiC_[LNNC], xsiC[DIM];
    for (int i = 0; i < LNNC; i++)
        wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[LNNC - 1]]->getINC();

    // data for computation of FEM fine mesh computations
    double xsi[DIM], phi_[LNN];
    double tARLQ_, djac_, weight_;
    double &alpha_f = parameters->getAlphaF();

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double ***ddphi_dx;
    ddphi_dx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphi_dx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphi_dx[i][j] = new double[LNN];
    };

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    // Defines the integration points adimentional coordinates
    for (int i = 0; i < DIM; i++)
        xsi[i] = nQuad.PointListFem(index, i);

    // Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    // Computes the shape functions
    shapeQuad.evaluateFem(xsi, phi_);

    // Computes the jacobian matrix
    getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

    // Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

    // Computes spatial second derivatives
    getSecondSpatialDerivatives_FEM(ainv_, ddphi_dx);

    // Integration point in the coarse mesh element
    for (int k = 0; k < DIM; k++)
        xsiC[k] = intPointCorrespXsi_FEM[index][k];

    // Computes the coqrse shape functions
    shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

    // Computes tARLQ_  (matriz norm)
    double wna_ = alpha_f * intPointWeightFunctionSpecial_FEM[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev_FEM[index];
    getParameterArlequinMN(LNN, LNNC, wna_, djac_, weight_, tARLQ_, phi_, phiC_, dphi_dx, ddphi_dx);

    // Computes matrixes and vectors
    getMatrixAndVectorsSameMeshArlqStab(LNN, wna_, djac_, weight_, tARLQ_, phi_, dphi_dx, ddphi_dx,
                                        arlequinStabD, arlequinStabVectorD,
                                        arlequinStab1, arlequinStabVector1);

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; ++j)
            delete[] ddphi_dx[i][j];
        delete[] ddphi_dx[i];
    };
    delete[] ddphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMeshArlqStab_FEM_FEM(int &index, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                  double **arlequinStabD, double *arlequinStabVectorD,
                                                                  double **arlequinStab1, double *arlequinStabVector1)
{
    int LNN = 4 * DIM - 2;
    int LNNC = 4 * DIM - 2;
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // Data for FEM coarse mesh computation
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    double phiC_[LNNC], xsiC[DIM];

    // data for computation of FEM fine mesh computations
    double xsi[DIM], phi_[LNN];
    double tARLQ_, djac_, weight_;
    double &alpha_f = parameters->getAlphaF();

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double ***ddphi_dx;
    ddphi_dx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphi_dx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphi_dx[i][j] = new double[LNN];
    };

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    // Defines the integration points adimentional coordinates
    for (int i = 0; i < DIM; i++) xsi[i] = nQuad.PointListFem(index, i);

    // Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    // Computes the shape functions
    shapeQuad.evaluateFem(xsi, phi_);

    // Computes the jacobian matrix
    getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

    // Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

    // Computes spatial second derivatives
    getSecondSpatialDerivatives_FEM(ainv_, ddphi_dx);

    // Integration point in the coarse mesh element
    for (int k = 0; k < DIM; k++)
        xsiC[k] = intPointCorrespXsi_FEM[index][k];

    // Computes the coqrse shape functions
    shapeQuad.evaluateFem(xsiC, phiC_);

    // Computes tARLQ_  (matriz norm)
    double wna_ = alpha_f * intPointWeightFunctionSpecial_FEM[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev_FEM[index];
    getParameterArlequinMN(LNN, LNNC, wna_, djac_, weight_, tARLQ_, phi_, phiC_, dphi_dx, ddphi_dx);

    // Computes matrixes and vectors
    getMatrixAndVectorsSameMeshArlqStab(LNN, wna_, djac_, weight_, tARLQ_, phi_, dphi_dx, ddphi_dx,
                                        arlequinStabD, arlequinStabVectorD,
                                        arlequinStab1, arlequinStabVector1);

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; ++j)
            delete[] ddphi_dx[i][j];
        delete[] ddphi_dx[i];
    };
    delete[] ddphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersSameMeshArlqStabLaplace_FEM_FEM(int &index, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                         double **arlequinStabD, double *arlequinStabVectorD,
                                                                         double **arlequinStab1, double *arlequinStabVector1)
{
    int LNN = 4 * DIM - 2;
    int LNNC = 4 * DIM - 2;
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // data for FEM coarse mesh computation
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    double phiC_[LNNC], xsiC[DIM];

    double **dphiC_dx;
    dphiC_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphiC_dx[i] = new double[LNN];

    double **ainvC_;
    ainvC_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainvC_[i] = new double[DIM];

    double ***ddphiC_dx;
    ddphiC_dx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphiC_dx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphiC_dx[i][j] = new double[LNN];
    };

    // data for FEM fine mesh computations
    double xsi[DIM], phi_[LNN];
    double tARLQ0_, tARLQ1_, djac_, weight_;
    double &alpha_f = parameters->getAlphaF();

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    double ***ddphi_dx;
    ddphi_dx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphi_dx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphi_dx[i][j] = new double[LNN];
    };

    // Defines the integration points adimentional coordinates
    for (int i = 0; i < DIM; i++)
        xsi[i] = nQuad.PointListFem(index, i);

    // Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    // Computes the shape functions
    shapeQuad.evaluateFem(xsi, phi_);

    // Computes the jacobian matrix
    getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

    // Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

    // Computes spatial second derivatives
    getSecondSpatialDerivatives_FEM(ainv_, ddphi_dx);

    // Integration point in the coarse mesh element
    for (int k = 0; k < DIM; k++)
        xsiC[k] = intPointCorrespXsi_FEM[index][k];

    // Computes the coarse shape functions
    shapeQuad.evaluateFem(xsiC, phiC_);

    // Computes the coarse jacobian matrix
    getJacobianMatrix_COARSE_FEM(xsiC, ainvC_);

    // Computes coarse spatial derivatives
    getSpatialDerivatives_FEM(xsiC, ainvC_, dphiC_dx);

    // Computes coarse spatial second derivatives
    getSecondSpatialDerivatives_FEM(ainvC_, ddphiC_dx);

    // Computes tARLQ_  (matriz norm)
    double wna_ = intPointWeightFunctionSpecial_FEM[index];
    getParameterArlequinLaplaceMN(LNN, LNNC, wna_, djac_, weight_, tARLQ0_, tARLQ1_,
                                  phi_, phiC_, dphi_dx, dphiC_dx, ddphi_dx, ddphiC_dx);

    // Computes matrixes and vectors
    getMatrixAndVectorsSameMeshArlqStabLaplace(LNN, wna_, djac_, weight_, tARLQ0_, tARLQ1_,
                                               phi_, dphi_dx, ddphi_dx, arlequinStabD, arlequinStabVectorD,
                                               arlequinStab1, arlequinStabVector1);

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; ++j)
            delete[] ddphi_dx[i][j];
        delete[] ddphi_dx[i];
    };
    delete[] ddphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; ++j)
            delete[] ddphiC_dx[i][j];
        delete[] ddphiC_dx[i];
    };
    delete[] ddphiC_dx;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMesh_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                               std::vector<IParameters *> &iparamC, int &ielem,
                                                               double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int LNN = 4 * DIM - 2;
    int LNNC = 18 * DIM - 27;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[LNNC], phiC_[LNNC], xsiC[DIM];
    for (int i = 0; i < LNNC; i++) wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[LNNC - 1]]->getINC();

    double **dphiC_dx;
    dphiC_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphiC_dx[i] = new double[LNNC];

    double **ainvC_;
    ainvC_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainvC_[i] = new double[DIM];

    // Data for FEM fine mesh computation (Lagrange field)
    double phi_[LNN], xsi[DIM];
    double djac_, weight_;

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC, ainvC_, dphiC_dx);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMesh(djac_, weight_, LNN, phi_, dphi_dx, LNNC, phiC_, dphiC_dx,
                                             lagrMultMatrix, rhsVector1, rhsVector2);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < DIM; ++i)
        delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainvC_[i];
    delete[] ainvC_;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMeshLaplace_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                      std::vector<IParameters *> &iparamC, int &ielem,
                                                                      double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int LNN = 4 * DIM - 2;
    int LNNC = 18 * DIM - 27;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[LNNC], phiC_[LNNC], xsiC[DIM];
    for (int i = 0; i < LNNC; i++)
        wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[LNNC - 1]]->getINC();

    // Data for FEM fine mesh computation (Lagrange field)
    double phi_[LNN], xsi[DIM];
    double djac_, weight_;

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMeshLaplace(djac_, weight_, LNN, phi_, LNNC, phiC_,
                                                    lagrMultMatrix, rhsVector1, rhsVector2);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMesh_FEM_FEM(std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                               int &ielem, double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int LNN = 4 * DIM - 2;
    int LNNC = 4 * DIM - 2;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    double xsiC[DIM], phiC_[LNNC];

    double **dphiC_dx;
    dphiC_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphiC_dx[i] = new double[LNNC];

    double **ainvC_;
    ainvC_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainvC_[i] = new double[DIM];

    double **JacC;
    JacC = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        JacC[i] = new double[DIM];

    // Data for FEM fine mesh computation (Lagrange field)
    double phi_[LNN], xsi[DIM];
    double djac_, weight_;

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsiC, phiC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC, ainvC_, dphiC_dx);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMesh(djac_, weight_, LNN, phi_, dphi_dx, LNNC, phiC_, dphiC_dx,
                                             lagrMultMatrix, rhsVector1, rhsVector2);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < DIM; ++i)
        delete[] JacC[i];
    delete[] JacC;

    for (int i = 0; i < DIM; ++i)
        delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMesh_ISO_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                               std::vector<IParameters *> &iparamC, int &ielem,
                                                               double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int LNN = 18 * DIM - 27;
    int LNNC = 18 * DIM - 27;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[LNNC], phiC_[LNNC], xsiC[DIM];
    for (int i = 0; i < LNNC; i++)
        wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[LNNC - 1]]->getINC();

    double **ainvC_;
    ainvC_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainvC_[i] = new double[DIM];

    double **dphiC_dx;
    dphiC_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphiC_dx[i] = new double[LNNC];

    // Data for IGA fine mesh computation (Lagrange field)
    double wpc[LNN], phi_[LNN], xsi[DIM];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **quadJacMat;
    quadJacMat = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        quadJacMat[i] = new double[DIM];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {

        if ((intPointCorrespElem_ISO[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListIso(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListIso(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

            // Computes the jacobian matrix
            getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_ISO[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC, ainvC_, dphiC_dx);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMesh(djac_, weight_, LNN, phi_, dphi_dx, LNNC, phiC_, dphiC_dx,
                                             lagrMultMatrix, rhsVector1, rhsVector2);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < DIM; ++i)
        delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] quadJacMat[i];
    delete[] quadJacMat;

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMeshLaplace_FEM_FEM(std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                      int &ielem, double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int LNN = 4 * DIM - 2;
    int LNNC = 4 * DIM - 2;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    double xsiC[DIM], phiC_[LNNC];

    // Data for FEM fine mesh computation (Lagrange field)
    double phi_[LNN], xsi[DIM];
    double djac_, weight_;

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsiC, phiC_);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMeshLaplace(djac_, weight_, LNN, phi_, LNNC, phiC_,
                                                    lagrMultMatrix, rhsVector1, rhsVector2);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMeshLaplace_ISO_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                      std::vector<IParameters *> &iparamC, int &ielem,
                                                                      double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int LNN = 18 * DIM - 27;
    int LNNC = 18 * DIM - 27;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[LNNC], phiC_[LNNC], xsiC[DIM];
    for (int i = 0; i < LNNC; i++)
        wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[LNNC - 1]]->getINC();

    // Data for IGA fine mesh computation (Lagrange field)
    double wpc[LNN], phi_[LNN], xsi[DIM];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    double djac_, weight_;

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **quadJacMat;
    quadJacMat = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        quadJacMat[i] = new double[DIM];

    int index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {

        if ((intPointCorrespElem_ISO[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListIso(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListIso(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

            // Computes the jacobian matrix
            getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_ISO[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMeshLaplace(djac_, weight_, LNN, phi_, LNNC, phiC_,
                                                    lagrMultMatrix, rhsVector1, rhsVector2);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] quadJacMat[i];
    delete[] quadJacMat;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMesh_FEM_ISO_TESTE(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                     std::vector<IParameters *> &iparamC, int &ielem,
                                                                     double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int LNN = 10;
    int LNNC = 27;

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[27], phiC_[27], xsiC[3];
    for (int i = 0; i < 27; i++)
        wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[26]]->getINC();

    // Data for FEM fine mesh computation (Lagrange field)
    double phi_[10], xsi[3];
    double djac_, weight_;

    double &alpha_f = parameters->getAlphaF();
    double &dens_ = parameters->getDensity();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            double **dphi_dx;
            dphi_dx = new double *[3];
            for (int i = 0; i < 3; ++i)
                dphi_dx[i] = new double[10];

            double **ainv_;
            ainv_ = new double *[3];
            for (int i = 0; i < 3; ++i)
                ainv_[i] = new double[3];

            double **Jac;
            Jac = new double *[3];
            for (int i = 0; i < 3; ++i)
                Jac[i] = new double[3];

            double **dphiC_dx;
            dphiC_dx = new double *[3];
            for (int i = 0; i < 3; ++i)
                dphiC_dx[i] = new double[27];

            double **ainvC_;
            ainvC_ = new double *[3];
            for (int i = 0; i < 3; ++i)
                ainvC_[i] = new double[3];

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < 3; k++)
                xsi[k] = nQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = nQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC, ainvC_, dphiC_dx);

            // Interpolated variables
            // Velocity
            double u_[3], uPrev_[3], una_[3];
            getInterpVelCoarse(LNNC, phiC_, u_, uPrev_);
            for (int i = 0; i < 3; i++)
                una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

            // Velocity Derivatives
            double **du_dx;
            du_dx = new double *[3];
            for (int i = 0; i < 3; ++i)
                du_dx[i] = new double[3];

            double **duPrev_dx;
            duPrev_dx = new double *[3];
            for (int i = 0; i < 3; ++i)
                duPrev_dx[i] = new double[3];

            double duna_dx[3][3];
            getInterpVelDerCoarse(LNNC, dphiC_dx, du_dx, duPrev_dx);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];

            // Lagrange Multiplieres
            double lambda_[3];
            getInterpLambda(LNN, phi_, lambda_);

            // Lagrange Multiplieres derivatives
            double **dlambda_dx;
            dlambda_dx = new double *[3];
            for (int i = 0; i < 3; ++i)
                dlambda_dx[i] = new double[3];
            getInterpLambdaDer(LNN, dphi_dx, dlambda_dx);

            double WJ = weight_ * djac_;

            for (int i = 0; i < 10; i++)
            {
                for (int j = 0; j < 27; j++)
                {
                    // L2 operator
                    double L2 = phi_[i] * phiC_[j] * k1;

                    lagrMultMatrix[3 * i][3 * j] += L2 * WJ;
                    lagrMultMatrix[3 * i + 1][3 * j + 1] += L2 * WJ;
                    lagrMultMatrix[3 * i + 2][3 * j + 2] += L2 * WJ;

                    // H1 operator
                    double H1xx = (dphi_dx[0][i] * dphiC_dx[0][j] +
                                   0.5 * dphi_dx[1][i] * dphiC_dx[1][j] +
                                   0.5 * dphi_dx[2][i] * dphiC_dx[2][j]) *
                                  k2;
                    double H1xy = 0.5 * dphi_dx[1][i] * dphiC_dx[0][j] * k2;
                    double H1xz = 0.5 * dphi_dx[2][i] * dphiC_dx[0][j] * k2;
                    double H1yx = 0.5 * dphi_dx[0][i] * dphiC_dx[1][j] * k2;
                    double H1yy = (dphi_dx[1][i] * dphiC_dx[1][j] +
                                   0.5 * dphi_dx[0][i] * dphiC_dx[0][j] +
                                   0.5 * dphi_dx[2][i] * dphiC_dx[2][j]) *
                                  k2;
                    double H1yz = 0.5 * dphi_dx[2][i] * dphiC_dx[1][j] * k2;
                    double H1zx = 0.5 * dphi_dx[0][i] * dphiC_dx[2][j] * k2;
                    double H1zy = 0.5 * dphi_dx[1][i] * dphiC_dx[2][j] * k2;
                    double H1zz = (dphi_dx[2][i] * dphiC_dx[2][j] +
                                   0.5 * dphi_dx[0][i] * dphiC_dx[0][j] +
                                   0.5 * dphi_dx[1][i] * dphiC_dx[1][j]) *
                                  k2;

                    lagrMultMatrix[3 * i][3 * j] += H1xx * WJ;
                    lagrMultMatrix[3 * i][3 * j + 1] += H1xy * WJ;
                    lagrMultMatrix[3 * i][3 * j + 2] += H1xz * WJ;
                    lagrMultMatrix[3 * i + 1][3 * j] += H1yx * WJ;
                    lagrMultMatrix[3 * i + 1][3 * j + 1] += H1yy * WJ;
                    lagrMultMatrix[3 * i + 1][3 * j + 2] += H1yz * WJ;
                    lagrMultMatrix[3 * i + 2][3 * j] += H1zx * WJ;
                    lagrMultMatrix[3 * i + 2][3 * j + 1] += H1zy * WJ;
                    lagrMultMatrix[3 * i + 2][3 * j + 2] += H1zy * WJ;
                };

                double H1Ux = dphi_dx[0][i] * duna_dx[0][0] +
                              0.5 * dphi_dx[1][i] * duna_dx[0][1] +
                              0.5 * dphi_dx[2][i] * duna_dx[0][2] +
                              0.5 * dphi_dx[1][i] * duna_dx[1][0] +
                              0.5 * dphi_dx[2][i] * duna_dx[2][0];
                double H1Uy = 0.5 * dphi_dx[0][i] * duna_dx[0][1] +
                              dphi_dx[1][i] * duna_dx[1][1] +
                              0.5 * dphi_dx[0][i] * duna_dx[1][0] +
                              0.5 * dphi_dx[2][i] * duna_dx[1][2] +
                              0.5 * dphi_dx[2][i] * duna_dx[2][1];
                double H1Uz = 0.5 * dphi_dx[0][i] * duna_dx[0][2] +
                              0.5 * dphi_dx[1][i] * duna_dx[1][2] +
                              0.5 * dphi_dx[0][i] * duna_dx[2][0] +
                              0.5 * dphi_dx[1][i] * duna_dx[2][1] +
                              dphi_dx[2][i] * duna_dx[2][2];

                rhsVector2[3 * i] -= (phi_[i] * una_[0] * k1 + H1Ux * k2) * WJ;
                rhsVector2[3 * i + 1] -= (phi_[i] * una_[1] * k1 + H1Uy * k2) * WJ;
                rhsVector2[3 * i + 2] -= (phi_[i] * una_[2] * k1 + H1Uz * k2) * WJ;
            };

            for (int i = 0; i < 27; i++)
            {

                double H1Lx = dphiC_dx[0][i] * dlambda_dx[0][0] +
                              0.5 * dphiC_dx[1][i] * dlambda_dx[0][1] +
                              0.5 * dphiC_dx[2][i] * dlambda_dx[0][2] +
                              0.5 * dphiC_dx[1][i] * dlambda_dx[1][0] +
                              0.5 * dphiC_dx[2][i] * dlambda_dx[2][0];
                double H1Ly = 0.5 * dphiC_dx[0][i] * dlambda_dx[0][1] +
                              dphiC_dx[1][i] * dlambda_dx[1][1] +
                              0.5 * dphiC_dx[0][i] * dlambda_dx[1][0] +
                              0.5 * dphiC_dx[2][i] * dlambda_dx[1][2] +
                              0.5 * dphiC_dx[2][i] * dlambda_dx[2][1];
                double H1Lz = 0.5 * dphiC_dx[0][i] * dlambda_dx[0][2] +
                              0.5 * dphiC_dx[1][i] * dlambda_dx[1][2] +
                              0.5 * dphiC_dx[0][i] * dlambda_dx[2][0] +
                              0.5 * dphiC_dx[1][i] * dlambda_dx[2][1] +
                              dphiC_dx[2][i] * dlambda_dx[2][2];

                rhsVector1[3 * i] -= (phiC_[i] * lambda_[0] * k1 + H1Lx * k2) * WJ;
                rhsVector1[3 * i + 1] -= (phiC_[i] * lambda_[1] * k1 + H1Ly * k2) * WJ;
                rhsVector1[3 * i + 2] -= (phiC_[i] * lambda_[2] * k1 + H1Lz * k2) * WJ;
            };

            // deallocating memory
            for (int i = 0; i < 3; ++i)
                delete[] du_dx[i];
            delete[] du_dx;

            for (int i = 0; i < 3; ++i)
                delete[] duPrev_dx[i];
            delete[] duPrev_dx;

            for (int i = 0; i < 3; ++i)
                delete[] dlambda_dx[i];
            delete[] dlambda_dx;

            for (int i = 0; i < 3; ++i)
                delete[] dphi_dx[i];
            delete[] dphi_dx;

            for (int i = 0; i < 3; ++i)
                delete[] ainv_[i];
            delete[] ainv_;

            for (int i = 0; i < 3; ++i)
                delete[] Jac[i];
            delete[] Jac;

            for (int i = 0; i < 3; ++i)
                delete[] dphiC_dx[i];
            delete[] dphiC_dx;

            for (int i = 0; i < 3; ++i)
                delete[] ainvC_[i];
            delete[] ainvC_;
        };

        index++;
    };

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                           std::vector<IParameters *> &iparamC, int &ielem,
                                                                           double **jacobianNRMatrix, double *rhsVector)
{

    int LNN = 4 * DIM - 2;
    int LNNC = 18 * DIM - 27;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[LNNC], phiC_[LNNC], xsiC[DIM];
    for (int i = 0; i < LNNC; i++)
        wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[LNNC - 1]]->getINC();

    double **dphiC_dx;
    dphiC_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphiC_dx[i] = new double[LNNC];

    double **ainvC_;
    ainvC_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainvC_[i] = new double[DIM];

    // Data for FEM fine mesh computation (Lagrange field)
    double xsi[DIM], phi_[LNN];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    double tSUPG_, tPSPG_, tLSIC_;
    double djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

            // Computes stabilization parameters
            getNewParameterSUPG_FEM(tSUPG_, tPSPG_, tLSIC_, Jac, phi_, dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC, ainvC_, dphiC_dx);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG(djac_, weight_, tSUPG_, tPSPG_, LNN, phi_, LNNC,
                                                         phiC_, dphiC_dx, jacobianNRMatrix, rhsVector);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < DIM; ++i)
        delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainvC_[i];
    delete[] ainvC_;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMeshArlqStabLaplace_FEM_FEM(std::vector<Nodes *> &nodesCoarse_, int *connecC, int &ielem,
                                                                              double **arlequinStabD, double *arlequinStabVectorD,
                                                                              double **arlequinStab0, double *arlequinStabVector0)
{

    int LNN = 4 * DIM - 2;
    int LNNC = 4 * DIM - 2;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    double phiC_[LNNC], xsiC[DIM];

    double **dphiC_dx;
    dphiC_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphiC_dx[i] = new double[LNN];

    double **ainvC_;
    ainvC_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainvC_[i] = new double[DIM];

    double ***ddphiC_dx;
    ddphiC_dx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphiC_dx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphiC_dx[i][j] = new double[LNN];
    };

    // Data for FEM fine mesh computation (Lagrange field)
    double xsi[DIM], phi_[LNN];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double ***ddphi_dx;
    ddphi_dx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphi_dx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphi_dx[i][j] = new double[LNN];
    };

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    double tARLQ0_, tARLQ1_, djac_, weight_;
    double &alpha_f = parameters->getAlphaF();

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

            // Computes spatial second derivatives
            getSecondSpatialDerivatives_FEM(ainv_, ddphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsiC, phiC_);

            // Computes the coarse jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC, ainvC_);

            // Computes coarse spatial derivatives
            getSpatialDerivatives_FEM(xsiC, ainvC_, dphiC_dx);

            // Computes coarse spatial second derivatives
            getSecondSpatialDerivatives_FEM(ainvC_, ddphiC_dx);

            // get Arlequin stabilization parameter
            double wna_ = intPointWeightFunctionSpecial_FEM[index];
            getParameterArlequinLaplaceMN(LNN, LNNC, wna_, djac_, weight_, tARLQ0_, tARLQ1_,
                                          phi_, phiC_, dphi_dx, dphiC_dx, ddphi_dx, ddphiC_dx);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMeshArlqStabLaplace(LNN, LNNC, wna_, djac_,
                                                            weight_, tARLQ0_, tARLQ1_,
                                                            phi_, dphi_dx, dphiC_dx, ddphiC_dx,
                                                            arlequinStabD, arlequinStabVectorD,
                                                            arlequinStab0, arlequinStabVector0);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; ++j)
            delete[] ddphi_dx[i][j];
        delete[] ddphi_dx[i];
    };
    delete[] ddphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; ++j)
            delete[] ddphiC_dx[i][j];
        delete[] ddphiC_dx[i];
    };
    delete[] ddphiC_dx;

    return;
};

template <int DIM>
void Element<DIM>::getLagrangeMultipliersDifferentMeshArlqStab_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                       std::vector<IParameters *> &iparamC, int &ielem,
                                                                       double **arlequinStabD, double *arlequinStabVectorD,
                                                                       double **arlequinStab0, double *arlequinStabVector0)
{
    int LNN = 4 * DIM - 2;
    int LNNC = 18 * DIM - 27;

    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapFunction shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[LNNC], phiC_[LNNC], xsiC[DIM];
    for (int i = 0; i < LNNC; i++)
        wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[LNNC - 1]]->getINC();

    double **dphiC_dx;
    dphiC_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphiC_dx[i] = new double[LNNC];

    double ***ddphiC_dx;
    ddphiC_dx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphiC_dx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphiC_dx[i][j] = new double[LNNC];
    };

    double **ainvC_;
    ainvC_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainvC_[i] = new double[DIM];

    // Data for FEM fine mesh computation (Lagrange field)
    double xsi[DIM], phi_[LNN];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double ***ddphi_dx;
    ddphi_dx = new double **[DIM];
    for (int i = 0; i < DIM; ++i)
    {
        ddphi_dx[i] = new double *[DIM];
        for (int j = 0; j < DIM; j++)
            ddphi_dx[i][j] = new double[LNN];
    };

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    double tARLQ_, djac_, weight_;
    double &alpha_f = parameters->getAlphaF();

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++)
                xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

            // Computes second spatial derivatives
            getSecondSpatialDerivatives_FEM(ainv_, ddphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < DIM; k++)
                xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC, ainvC_, dphiC_dx);

            // Computes second spatial derivatives
            getSecondSpatialDerivatives_COARSE_ISO(xsiC, ainvC_, ddphiC_dx);

            // get Arlequin stabilization parameter
            double wna_ = alpha_f * intPointWeightFunctionSpecial_FEM[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev_FEM[index];
            getParameterArlequinMN_COARSE(LNN, LNNC, wna_, djac_, weight_, tARLQ_, phi_, phiC_,
                                          dphi_dx, dphiC_dx, ddphi_dx, ddphiC_dx);

            // Computes Matrix and vectors
            getMatrixAndVectorsDifferentMeshArlqStab(LNN, LNNC, wna_, djac_, weight_, tARLQ_,
                                                     dphi_dx, phiC_, dphiC_dx, ddphiC_dx,
                                                     arlequinStabD, arlequinStabVectorD,
                                                     arlequinStab0, arlequinStabVector0);
        };
        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < DIM; ++i)
        delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j < DIM; ++j)
        {
            ddphiC_dx[i][j];
            ddphi_dx[i][j];
        }
        delete[] ddphiC_dx[i];
        delete[] ddphi_dx[i];
    }
    delete[] ddphiC_dx;
    delete[] ddphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainvC_[i];
    delete[] ainvC_;

    return;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//----------------MOUNTS EACH TYPE OF INCOMPRESSIBLE FLOW PROBEM----------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------

template <int DIM>
void Element<DIM>::getSteadyElasticity_FEM(double **jacobianNRMatrix)
{

    int LNN = 4 * DIM - 2;
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    double phi_[LNN], xsi[DIM];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++) xsi[i] = nQuad.PointListFem(index, i);

        // Returns the quadrature integration weight
        double weight_ = nQuad.WeightListFem(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi, phi_);

        // Computes the jacobian matrix
        double djac_;
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);
        
        // Set element mesh moving parameters
        double jM0 = 1.e+10;
        double chi = 1.;
        double meshMovPar = (jM0/djac_);

        // Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

        // Computes the element matrix
        getElemMatrixElas(LNN, djac_, weight_, meshMovPar, dphi_dx,jacobianNRMatrix);

        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;
};

template <int DIM>
void Element<DIM>::setDirichletConstrain(int &LNN, double **jacobianNRMatrix)
{

    for (int i = 0; i < LNN; i++)
    {
        if ((*nodes_)[connect_[i]]->getConstrains(0) == 1)
        {

            for (int j = 0; j < LNN; j++)
            {
                jacobianNRMatrix[i][j] = 0.;
            };
            jacobianNRMatrix[i][i] = 1.;
        };
    };
};

template <int DIM>
void Element<DIM>::getTNS_FEM(int &iTime, double **jacobianNRMatrix, double *rhsVector)
{

    int LNN = 4 * DIM - 2;
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // variables
    double &alpha_f = parameters->getAlphaF();
    double tSUPG_, tPSPG_, tLSIC_;
    double djac_, weight_, xsi[DIM], phi_[LNN];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    int index = 0;
    ip_ = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListFem(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListFem(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi, phi_);

        // Computes the jacobian matrixgetElem
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

        // Compute Stabilization Parameters
        //getNewParameterSUPG_FEM(tSUPG_, tPSPG_, tLSIC_, Jac, phi_, dphi_dx);
        getParameterSUPG_FEM(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

        // Computes the element matrix
        double wna_ = alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];

        getElemMatrixTNS(LNN, wna_, djac_, weight_, tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx, jacobianNRMatrix);

        // Computes the RHS vector
        getResidualVectorTNS(iTime, LNN, wna_, djac_, weight_, tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx, rhsVector);

        index++;
        ip_++;
    };

    for (int i = 0; i < DIM; ++i){
        delete[] dphi_dx[i];
        delete[] ainv_[i];
        delete[] Jac[i];
    }
        
    delete[] dphi_dx;
    delete[] ainv_;
    delete[] Jac;
};

template <int DIM>
void Element<DIM>::computeNormErrorLaplace_FEM(double &error)
{

    int LNN = 4 * DIM - 2;
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // variables
    double djac_, weight_, xsi[DIM], phi_[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    error = 0.0;

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListFem(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListFem(index);

        // Computes the shape functions
        shapeQuad.evaluateFem(xsi, phi_);

        // Physical coordinates
        double x_[DIM], xPrev_[DIM];
        getInterpCoord(LNN, phi_, x_, xPrev_);

        // analict solution
        double force = std::exp(x_[0]) * (-1. + x_[0]) * x_[0] * (-1. + x_[1]) * x_[1];

        // Computes the jacobian matrix
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

        // Interpolated solution
        double vel[DIM] = {};
        for (int i = 0; i < LNN; i++)
            for (int j = 0; j < DIM; j++)
                vel[j] += (*nodes_)[connect_[i]]->getVelocityArlequin(j) * phi_[i];

        error += djac_ * weight_ * (vel[0] - force);

        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;
};

template <int DIM>
void Element<DIM>::computeNormError_FEM(double &errorU, double &errorP)
{

    int LNN = 4 * DIM - 2;
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // variables
    double djac_, weight_, xsi[DIM], phi_[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    errorU = 0.0;
    errorP = 0.0;

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListFem(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListFem(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi, phi_);

        // physical coordinates
        double coord[DIM] = {};
        for (int i = 0; i < LNN; i++)
            for (int j = 0; j < DIM; j++)
                coord[j] += (*nodes_)[connect_[i]]->getCoordinateValue(j) * phi_[j];

        // analitic solution
        double uA_ = cos(2. * pi * coord[0]) * sin(2. * pi * coord[1]);
        double vA_ = -sin(2. * pi * coord[0]) * cos(2. * pi * coord[1]);
        double pA_ = (1. - 8. * pi * pi * 0.01) * cos(2. * pi * coord[0]) * sin(2. * pi * coord[1]);

        // Computes the jacobian matrix
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

        // Interpolated vel and press
        double vel[DIM] = {};
        for (int i = 0; i < LNN; i++)
            for (int j = 0; j < DIM; j++)
                vel[j] += (*nodes_)[connect_[i]]->getVelocityArlequin(j) * phi_[i];

        double press = 0.0;
        for (int i = 0; i < LNN; i++)
            press += (*nodes_)[connect_[i]]->getPressureArlequin() * phi_[i];

        double normU = sqrt((vel[0] - uA_) * (vel[0] - uA_) + (vel[1] - vA_) * (vel[1] - vA_));
        errorU += djac_ * weight_ * normU * normU;
        errorP += djac_ * weight_ * (press - pA_) * (press - pA_);

        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;
};

template <int DIM>
void Element<DIM>::getLaplace_FEM(double **jacobianNRMatrix, double *rhsVector)
{

    int LNN = 4 * DIM - 2;
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // variables
    double djac_, weight_, xsi[DIM], phi_[LNN];

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[DIM];

    double **Jac;
    Jac = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        Jac[i] = new double[DIM];

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListFem(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListFem(index);

        // Computes the basis functions
        shapeQuad.evaluateFem(xsi, phi_);

        // Computes the jacobian matrixgetElem
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

        // Computes the element matrix
        double wna_ = intPointWeightFunction_FEM[index];
        getElemMatrixLaplace(LNN, wna_, djac_, weight_, dphi_dx, jacobianNRMatrix);

        // Computes the RHS vector
        getResidualVectorLaplace(LNN, wna_, djac_, weight_, phi_, dphi_dx, rhsVector);

        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] Jac[i];
    delete[] Jac;
};

template <int DIM>
void Element<DIM>::getLaplace_ISO(double **jacobianNRMatrix, double *rhsVector)
{
    int LNN = 18 * DIM - 27;
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // variables
    double djac_, weight_;

    double wpc[LNN], phi_[LNN], xsi[DIM];
    for (int i = 0; i < LNN; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[LNN];

    double **quadJacMat;
    quadJacMat = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        quadJacMat[i] = new double[DIM];

    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++)
            xsi[i] = nQuad.PointListIso(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        // Computes the jacobian matrix
        getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

        // Computes the element matrix
        double wna_ = intPointWeightFunction_ISO[index];
        getElemMatrixLaplace(LNN, wna_, djac_, weight_, dphi_dx, jacobianNRMatrix);

        index++;
    };

    // Computes the RHS vector
    // getResidualVectorLaplace(LNN,rhsVector);

    // setDirichletConstrain(LNN, jacobianNRMatrix);

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] quadJacMat[i];
    delete[] quadJacMat;
};

template <int DIM>
void Element<DIM>::getTNS_ISO(int &iTime, double **jacobianNRMatrix, double *rhsVector)
{

    int LNN = 18 * DIM - 27;
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapFunction shapeQuad;

    // variables
    double tSUPG_, tPSPG_, tLSIC_;
    double djac_, weight_;
    double &alpha_f = parameters->getAlphaF();

    double wpc[LNN], phi_[LNN], xsi[DIM];
    for (int i = 0; i < LNN; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[LNN - 1]]->getINC();

    double **dphi_dx;
    dphi_dx = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        dphi_dx[i] = new double[LNN];

    double **ainv_;
    ainv_ = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv_[i] = new double[LNN];

    double **quadJacMat;
    quadJacMat = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        quadJacMat[i] = new double[DIM];

    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        // Defines the integration points adimentional coordinates
        for (int i = 0; i < DIM; i++) xsi[i] = nQuad.PointListIso(index, i);

        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        // Computes the jacobian matrix
        getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

        // Compute Stabilization Parameters
        getNewParameterSUPG_ISO(tSUPG_, tPSPG_, tLSIC_, quadJacMat, phi_, dphi_dx);

        // Computes the element matrix
        double wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];
        getElemMatrixTNS(LNN, wna_, djac_, weight_, tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx, jacobianNRMatrix);

        // Computes the RHS vector
        getResidualVectorTNS(iTime, LNN, wna_, djac_, weight_, tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx, rhsVector);

        index++;
    };

    for (int i = 0; i < DIM; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < DIM; ++i)
        delete[] quadJacMat[i];
    delete[] quadJacMat;
};

//------------------------------------------------------------------------------
//-----------------------COMPUTE DRAG AND LIFT FORCES --------------------------
//------------------------------------------------------------------------------
template <>
void Element<2>::computeDragAndLiftForces_FEM(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce,
                                              double &dForce, double &lForce, double &aux_Mom, double &aux_Per)
{
    int nodesb_[3];

    if (side == 0)
    {
        for (int i = 0; i < 2; i++)
        {
            nodesb_[0]= connect_[1];        
            nodesb_[1]= connect_[4];
            nodesb_[2]= connect_[2];
        };
    }
    else
    {
        if (side == 1)
        {
            for (int i = 0; i < 2; i++)
            {
                nodesb_[0]= connect_[2];
                nodesb_[1]= connect_[5];
                nodesb_[2]= connect_[0];
            };
        }
        else
        {
            for (int i = 0; i < 2; i++)
            {
                nodesb_[0]= connect_[0];
                nodesb_[1]= connect_[3];
                nodesb_[2]= connect_[1];
            };
        };
    };

    BoundaryIntegQuadrature<2> bQuad;
    QuadShapeFunction<2> shapeQuad;
    BoundShapeFunction<2> shapeBound;
    double &visc_ = parameters->getViscosity();

    double n_vector[2] = {};
    double shearStress[2][2] = {};
    double ident[2][2];
    double load_friction[2] = {};
    double load_pressure[2] = {};

    ident[0][0] = 1.;
    ident[1][1] = 1.;

    double phi_[6], phiB_[3], xsi[2], xsiB[1], weightB;

    double **dphiB_;
    dphiB_ = new double *[1];
    for (int i = 0; i < 1; ++i)
        dphiB_[i] = new double[3];

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i)
        dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i)
        ainv_[i] = new double[2];

    double **Jac;
    Jac = new double *[2];
    for (int i = 0; i < 2; ++i)
        Jac[i] = new double[2];

    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i)
        du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i)
        duPrev_dx[i] = new double[2];

    double moment = 0.;
    double per = 0.;
    int index = 0;

    for (double *it = bQuad.begin(); it != bQuad.end(); it++)
    {

        xsiB[0] = bQuad.PointList(index, 0);
        weightB = bQuad.WeightList(index);

        if (side == 2)
        {
            xsi[0] = (-xsiB[0] + 1.) / 2.;
            xsi[1] = 0.;
        };
        if (side == 1)
        {
            xsi[1] = (xsiB[0] + 1.) / 2.;
            xsi[0] = 0.;
        };
        if (side == 0)
        {
            xsi[0] = (xsiB[0] + 1.) / 2.;
            xsi[1] = 1. - xsi[0];
        };

        // Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi, phi_);

        // Computes the jacobian matrix
        double djac_;
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

        int LNN = 6;
        getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);

        double press_;
        getInterpPress(LNN, phi_, press_);

        shapeBound.evaluateBoundaryFem(xsiB, phiB_);
        shapeBound.evaluateGradientBoundaryFem(xsiB, dphiB_);

        // computes the normal vector in the linear element
        double Tx = 0.;
        double Ty = 0.;
        for (int i = 0; i < 3; i++)
        {
            Tx += (*nodes_)[nodesb_[i]]->getCoordinateValue(0) * dphiB_[0][i];
            Ty += (*nodes_)[nodesb_[i]]->getCoordinateValue(1) * dphiB_[0][i];
        };

        double jacb_ = sqrt(Tx * Tx + Ty * Ty);

        n_vector[0] = Ty / jacb_;
        n_vector[1] = -Tx / jacb_;

        shearStress[0][0] = 2. * visc_ * du_dx[0][0];
        shearStress[0][1] = visc_ * (du_dx[0][1] + du_dx[1][0]);
        shearStress[1][0] = visc_ * (du_dx[0][1] + du_dx[1][0]);
        shearStress[1][1] = 2. * visc_ * du_dx[1][1];

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                load_pressure[i] += -press_ * ident[i][j] * n_vector[j] * jacb_ * weightB;
                load_friction[i] += shearStress[i][j] * n_vector[j] * jacb_ * weightB;
            }
        }

        // moment += ((-p_ + shearStress[0][0] + shearStress[1][0]) * y_
        //          -(-p_ + shearStress[0][1] + shearStress[1][1]) * (x_ - 0.248792267683901))
        //          * jacb_ * weightB;
        per += jacb_ * weightB;

        index++;

    };

    aux_Per = per;
    aux_Mom = moment;

    pDForce = -load_pressure[0];
    pLForce = -load_pressure[1];

    fDForce = -load_friction[0];
    fLForce = -load_friction[1];

    dForce = pDForce + fDForce;
    lForce = pLForce + fLForce;

    for (int i = 0; i < 2; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 1; ++i)
        delete[] dphiB_[i];
    delete[] dphiB_;

    for (int i = 0; i < 2; ++i)
        delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i)
        delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < 2; ++i)
        delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i)
        delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    return;
};

template <>
void Element<3>::computeDragAndLiftForces_FEM(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce,
                                              double &dForce, double &lForce, double &aux_Mom, double &aux_Per)
{

    double &visc_ = parameters->getViscosity();
    double localNodesBoundary_[6][3];

    if (side == 0)
    {
        for (int i = 0; i < 3; i++)
        {
            localNodesBoundary_[0][i] = (*nodes_)[connect_[3]]->getCoordinateValue(i);
            localNodesBoundary_[1][i] = (*nodes_)[connect_[1]]->getCoordinateValue(i);
            localNodesBoundary_[2][i] = (*nodes_)[connect_[2]]->getCoordinateValue(i);
            localNodesBoundary_[3][i] = (*nodes_)[connect_[9]]->getCoordinateValue(i);
            localNodesBoundary_[4][i] = (*nodes_)[connect_[5]]->getCoordinateValue(i);
            localNodesBoundary_[5][i] = (*nodes_)[connect_[8]]->getCoordinateValue(i);
        };
    };

    if (side == 1)
    {
        for (int i = 0; i < 3; i++)
        {
            localNodesBoundary_[0][i] = (*nodes_)[connect_[3]]->getCoordinateValue(i);
            localNodesBoundary_[1][i] = (*nodes_)[connect_[2]]->getCoordinateValue(i);
            localNodesBoundary_[2][i] = (*nodes_)[connect_[0]]->getCoordinateValue(i);
            localNodesBoundary_[3][i] = (*nodes_)[connect_[8]]->getCoordinateValue(i);
            localNodesBoundary_[4][i] = (*nodes_)[connect_[6]]->getCoordinateValue(i);
            localNodesBoundary_[5][i] = (*nodes_)[connect_[7]]->getCoordinateValue(i);
        };
    };

    if (side == 2)
    {
        for (int i = 0; i < 3; i++)
        {
            localNodesBoundary_[0][i] = (*nodes_)[connect_[3]]->getCoordinateValue(i);
            localNodesBoundary_[1][i] = (*nodes_)[connect_[0]]->getCoordinateValue(i);
            localNodesBoundary_[2][i] = (*nodes_)[connect_[1]]->getCoordinateValue(i);
            localNodesBoundary_[3][i] = (*nodes_)[connect_[7]]->getCoordinateValue(i);
            localNodesBoundary_[4][i] = (*nodes_)[connect_[4]]->getCoordinateValue(i);
            localNodesBoundary_[5][i] = (*nodes_)[connect_[9]]->getCoordinateValue(i);
        };
    };

    if (side == 3)
    {
        for (int i = 0; i < 3; i++)
        {
            localNodesBoundary_[0][i] = (*nodes_)[connect_[2]]->getCoordinateValue(i);
            localNodesBoundary_[1][i] = (*nodes_)[connect_[1]]->getCoordinateValue(i);
            localNodesBoundary_[2][i] = (*nodes_)[connect_[0]]->getCoordinateValue(i);
            localNodesBoundary_[3][i] = (*nodes_)[connect_[5]]->getCoordinateValue(i);
            localNodesBoundary_[4][i] = (*nodes_)[connect_[4]]->getCoordinateValue(i);
            localNodesBoundary_[5][i] = (*nodes_)[connect_[6]]->getCoordinateValue(i);
        };
    };

    BoundaryIntegQuadrature<3> bQuad;
    QuadShapeFunction<3> shapeQuad;
    BoundShapeFunction<3> shapeBound;

    double n_vector[3] = {};
    double shearStress[3][3] = {};
    double ident[3][3] = {};
    double load_friction[3] = {};
    double load_pressure[3] = {};

    ident[0][0] = 1.;
    ident[1][1] = 1.;
    ident[2][2] = 1.;

    double phi_[10], phiB_[6], xsi[3], xsiB[2], weightB, djac_, press_;

    double **dphiB_;
    dphiB_ = new double *[2];
    for (int i = 0; i < 2; ++i)
        dphiB_[i] = new double[6];

    double **dphi_dx;
    dphi_dx = new double *[3];
    for (int i = 0; i < 3; ++i)
        dphi_dx[i] = new double[10];

    double **ainv_;
    ainv_ = new double *[3];
    for (int i = 0; i < 3; ++i)
        ainv_[i] = new double[3];

    double **Jac;
    Jac = new double *[3];
    for (int i = 0; i < 3; ++i)
        Jac[i] = new double[3];

    double **du_dx;
    du_dx = new double *[3];
    for (int i = 0; i < 3; ++i)
        du_dx[i] = new double[3];

    double **duPrev_dx;
    duPrev_dx = new double *[3];
    for (int i = 0; i < 3; ++i)
        duPrev_dx[i] = new double[3];

    double moment = 0.;
    double per = 0.;
    int index = 0;

    for (double *it = bQuad.begin(); it != bQuad.end(); it++)
    {

        weightB = bQuad.WeightList(index);
        xsiB[0] = bQuad.PointList(index, 0);
        xsiB[1] = bQuad.PointList(index, 1);

        if (side == 0)
        {
            xsi[0] = 0.;
            xsi[1] = xsiB[0];
            xsi[2] = xsiB[1];
        };
        if (side == 1)
        {
            xsi[0] = xsiB[1];
            xsi[1] = 0.;
            xsi[2] = xsiB[0];
        };
        if (side == 2)
        {
            xsi[0] = xsiB[0];
            xsi[1] = xsiB[1];
            xsi[2] = 0.;
        };
        if (side == 3)
        {
            xsi[0] = xsiB[1];
            xsi[1] = xsiB[0];
            xsi[2] = 1. - xsiB[0] - xsiB[1];
        };

        // Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi, phi_);

        // Computes the jacobian matrix
        getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

        int LNN = 10;
        getInterpVelDer(LNN, dphi_dx, du_dx, duPrev_dx);

        getInterpPress(LNN, phi_, press_);

        shapeBound.evaluateBoundaryFem(xsiB, phiB_);
        shapeBound.evaluateGradientBoundaryFem(xsiB, dphiB_);

        double dx_dxsi1 = 0.;
        double dx_dxsi2 = 0.;
        double dy_dxsi1 = 0.;
        double dy_dxsi2 = 0.;
        double dz_dxsi1 = 0.;
        double dz_dxsi2 = 0.;

        for (int i = 0; i < 6; i++)
        {

            dx_dxsi1 += localNodesBoundary_[i][0] * dphiB_[0][i];
            dx_dxsi2 += localNodesBoundary_[i][0] * dphiB_[1][i];
            dy_dxsi1 += localNodesBoundary_[i][1] * dphiB_[0][i];
            dy_dxsi2 += localNodesBoundary_[i][1] * dphiB_[1][i];
            dz_dxsi1 += localNodesBoundary_[i][2] * dphiB_[0][i];
            dz_dxsi2 += localNodesBoundary_[i][2] * dphiB_[1][i];
        };

        double Tx = dy_dxsi1 * dz_dxsi2 - dy_dxsi2 * dz_dxsi1;
        double Ty = dz_dxsi1 * dx_dxsi2 - dz_dxsi2 * dx_dxsi1;
        double Tz = dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1;

        double jacs_ = sqrt(Tx * Tx + Ty * Ty + Tz * Tz);

        n_vector[0] = Tx / jacs_;
        n_vector[1] = Ty / jacs_;
        n_vector[2] = Tz / jacs_;

        shearStress[0][0] = 2. * visc_ * du_dx[0][0];
        shearStress[0][1] = visc_ * (du_dx[0][1] + du_dx[1][0]);
        shearStress[0][2] = visc_ * (du_dx[0][2] + du_dx[2][0]);
        shearStress[1][0] = visc_ * (du_dx[0][1] + du_dx[1][0]);
        shearStress[1][1] = 2. * visc_ * du_dx[1][1];
        shearStress[1][2] = visc_ * (du_dx[2][1] + du_dx[1][2]);
        shearStress[2][0] = visc_ * (du_dx[0][2] + du_dx[2][0]);
        shearStress[2][1] = visc_ * (du_dx[2][1] + du_dx[1][2]);
        shearStress[2][2] = 2. * visc_ * du_dx[2][2];

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                load_pressure[i] += -press_ * ident[i][j] * n_vector[j] * jacs_ * weightB;
                load_friction[i] += shearStress[i][j] * n_vector[j] * jacs_ * weightB;
            };
        };

        // moment +=  (-(load_pressure[0]+load_friction[0]) * (y_- 4.)
        // 		   + (load_pressure[0]+load_friction[0]) * (x_ - 4.)) * jacs_ * weightB[0];

        per += jacs_ * weightB;

        index++;
    };

    // perimeter is this case is cylinder surface area
    aux_Per = per;
    aux_Mom = moment;

    pDForce = -load_pressure[0];
    pLForce = -load_pressure[1];

    fDForce = -load_friction[0];
    fLForce = -load_friction[1];

    dForce = pDForce + fDForce;
    lForce = pLForce + fLForce;

    for (int i = 0; i < 3; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;
    for (int i = 0; i < 2; ++i)
        delete[] dphiB_[i];
    delete[] dphiB_;
    for (int i = 0; i < 3; ++i)
        delete[] ainv_[i];
    delete[] ainv_;
    for (int i = 0; i < 3; ++i)
        delete[] Jac[i];
    delete[] Jac;

    return;
};

template class Element<2>;
template class Element<3>;
