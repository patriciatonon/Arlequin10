#include "Element.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------WEIGHTED INTEGRATION POINTS-------------------------
//------------------------------------------------------------------------------
template <>
void Element<2>::setIntegPointWeightFunction_FEM()
{

    int dim = 2;
    QuadShapeFunction<2> shapeQuad;

    // NORMAL QUADRATURE
    NormalQuad nQuad = NormalQuad();
    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {
        intPointWeightFunctionPrev_FEM[index] = intPointWeightFunction_FEM[index];
        intPointWeightFunction_FEM[index] = 0.;
        index++;
    };

    double xsi[dim], phi_[6];
    index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        for (int i = 0; i < dim; i++) xsi[i] = nQuad.PointListFem(index, i);
        shapeQuad.evaluateFem(xsi, phi_);

        for (int i = 0; i < 6; i++)
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

        for (int i = 0; i < dim; i++) xsi[i] = sQuad.PointListFem(index, i);
        shapeQuad.evaluateFem(xsi, phi_);

        for (int i = 0; i < 6; i++)
        {
            intPointWeightFunctionSpecial_FEM[index] += (*nodes_)[connect_[i]]->getWeightFunction() * phi_[i];
        };

        index++;
    };
};

template <>
void Element<2>::setIntegPointWeightFunction_ISO()
{

    int dim = 2;
    QuadShapeFunction<2> shapeQuad;

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
    double wpc[9], phi_[9], xsi[dim];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        for (int i = 0; i < dim; i++) xsi[i] = nQuad.PointListIso(index, i);
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        for (int i = 0; i < 9; i++)
        {
            intPointWeightFunction_ISO[index] += (*nodes_)[connect_[i]]->getWeightFunction() * phi_[i];
        };

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

        for (int i = 0; i < dim; i++) xsi[i] = sQuad.PointListIso(index, i);
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        for (int i = 0; i < 9; i++)
        {
            intPointWeightFunctionSpecial_ISO[index] += (*nodes_)[connect_[i]]->getWeightFunction() * phi_[i];
        };

        index++;
    };
};

//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template <>
void Element<2>::getJacobianMatrix_FEM(double &djac_, double *xsi, double **Jac, double **ainv_)
{

    // Computes the spatial Jacobian matrix and its inverse
    QuadShapeFunction<2> shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[6];

    shapeQuad.evaluateGradientFem(xsi, dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    for (int i = 0; i < 6; i++)
    {

        double xna_ = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(0) +
                      (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(1) +
                      (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(1);

        dx_dxsi1 += xna_ * dphi[0][i];
        dx_dxsi2 += xna_ * dphi[1][i];
        dy_dxsi1 += yna_ * dphi[0][i];
        dy_dxsi2 += yna_ * dphi[1][i];
    };

    // Defining the Jacobian matrix
    Jac[0][0] = dx_dxsi1;
    Jac[0][1] = dx_dxsi2;
    Jac[1][0] = dy_dxsi1;
    Jac[1][1] = dy_dxsi2;

    // Computing the jacobian determinant
    djac_ = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];

    // Computing Jacobian inverse (transposed)
    // dxsi1/dx dxsi2/dx
    // dxsi1/dy dxsi2/dy
    ainv_[0][0] = dy_dxsi2 / djac_;
    ainv_[0][1] = -dy_dxsi1 / djac_;
    ainv_[1][0] = -dx_dxsi2 / djac_;
    ainv_[1][1] = dx_dxsi1 / djac_;

    for (int i = 0; i < 2; ++i) delete[] dphi[i];
    delete[] dphi;

    return;
};

template <>
void Element<2>::getJacobianMatrix_COARSE_FEM(double *xsi, double **Jac, double **ainv_)
{

    // Computes the spatial Jacobian matrix and its inverse
    QuadShapeFunction<2> shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[6];

    shapeQuad.evaluateGradientFem(xsi, dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    for (int i = 0; i < 6; i++)
    {

        double xna_ = alpha_f * (*nodesC_)[connectC_[i]]->getCoordinateValue(0) +
                      (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodesC_)[connectC_[i]]->getCoordinateValue(1) +
                      (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousCoordinateValue(1);

        dx_dxsi1 += xna_ * dphi[0][i];
        dx_dxsi2 += xna_ * dphi[1][i];
        dy_dxsi1 += yna_ * dphi[0][i];
        dy_dxsi2 += yna_ * dphi[1][i];
    };

    // Defining the Jacobian matrix
    Jac[0][0] = dx_dxsi1;
    Jac[0][1] = dx_dxsi2;
    Jac[1][0] = dy_dxsi1;
    Jac[1][1] = dy_dxsi2;

    // Computing the jacobian determinant
    double djacC_ = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];

    // Computing Jacobian inverse (transposed)
    ainv_[0][0] = dy_dxsi2 / djacC_;
    ainv_[0][1] = -dy_dxsi1 / djacC_;
    ainv_[1][0] = -dx_dxsi2 / djacC_;
    ainv_[1][1] = dx_dxsi1 / djacC_;

    for (int i = 0; i < 2; ++i) delete[] dphi[i];
    delete[] dphi;

    return;
};

template <>
void Element<2>::getJacobianMatrix_ISO(double &djac_, double *xsi, double **quadJacMat, double **ainv_)
{

    QuadShapeFunction<2> shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];

    int *inc_ = (*nodes_)[connect_[8]]->getINC();
    double wpc[9];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparameters), Npatch_);

    // Computes the Jacobian matrix - dx/dxsi
    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    for (int i = 0; i < 9; i++)
    {

        double xna_ = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(0) +
                      (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(1) +
                      (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(1);

        dx_dxsi1 += xna_ * dphi[0][i] / wpc[i];
        dx_dxsi2 += xna_ * dphi[1][i] / wpc[i];
        dy_dxsi1 += yna_ * dphi[0][i] / wpc[i];
        dy_dxsi2 += yna_ * dphi[1][i] / wpc[i];
    };

    // Computing the jacobian determinant
    double djacpp = dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1;

    // Computing jacobian inverse matrix (transposed) - dxsi/dx
    ainv_[0][0] = dy_dxsi2 / djacpp;
    ainv_[0][1] = -dy_dxsi1 / djacpp;
    ainv_[1][0] = -dx_dxsi2 / djacpp;
    ainv_[1][1] = dx_dxsi1 / djacpp;

    // Computing the quadrature jacobian matrix - dx/dxqsi
    int degm = (*iparameters)[Npatch_]->getDegree(0);
    int degn = (*iparameters)[Npatch_]->getDegree(1);
    int npcm = (*iparameters)[Npatch_]->getNcp(0);
    int npcn = (*iparameters)[Npatch_]->getNcp(1);
    double *uknot_ = (*iparameters)[Npatch_]->getuKnot();
    double *vknot_ = (*iparameters)[Npatch_]->getvKnot();

    int uind = inc_[0];
    int vind = inc_[1];

    double u1 = uknot_[uind];
    double u2 = uknot_[uind + 1];
    double v1 = vknot_[vind];
    double v2 = vknot_[vind + 1];

    double dxsi1_dqxsi1 = (u2 - u1) * 0.5;
    double dxsi1_dqxsi2 = 0.;
    double dxsi2_dqxsi1 = 0.;
    double dxsi2_dqxsi2 = (v2 - v1) * 0.5;

    quadJacMat[0][0] = dxsi1_dqxsi1 * dx_dxsi1 + dxsi2_dqxsi1 * dx_dxsi2;
    quadJacMat[0][1] = dxsi1_dqxsi2 * dx_dxsi1 + dxsi2_dqxsi2 * dx_dxsi2;
    quadJacMat[1][0] = dxsi1_dqxsi1 * dy_dxsi1 + dxsi2_dqxsi1 * dy_dxsi2;
    quadJacMat[1][1] = dxsi1_dqxsi2 * dy_dxsi1 + dxsi2_dqxsi2 * dy_dxsi2;

    djac_ = quadJacMat[0][0] * quadJacMat[1][1] - quadJacMat[0][1] * quadJacMat[1][0];

    for (int i = 0; i < 2; ++i) delete[] dphi[i];
    delete[] dphi;

    return;
};

template <>
void Element<2>::getJacobianMatrix_COARSE_ISO(double *xsi, double **quadJacMat, double **ainv_)
{

    QuadShapeFunction<2> shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    double **dphi;
    dphi = new double *[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];

    int *inc_ = (*nodesC_)[connectC_[8]]->getINC();
    double wpc[9];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodesC_)[connectC_[i]]->getWeightPC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparametersC), NpatchC_);

    // Computes the Jacobian matrix - dx/dxsi
    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    for (int i = 0; i < 9; i++)
    {

        double xna_ = alpha_f * (*nodesC_)[connectC_[i]]->getCoordinateValue(0) +
                      (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodesC_)[connectC_[i]]->getCoordinateValue(1) +
                      (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousCoordinateValue(1);

        dx_dxsi1 += xna_ * dphi[0][i] / wpc[i];
        dx_dxsi2 += xna_ * dphi[1][i] / wpc[i];
        dy_dxsi1 += yna_ * dphi[0][i] / wpc[i];
        dy_dxsi2 += yna_ * dphi[1][i] / wpc[i];
    };

    // Computing the jacobian determinant
    double djacpp = dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1;

    // Computing jacobian inverse matrix (transposed) - dxsi/dx
    ainv_[0][0] = dy_dxsi2 / djacpp;
    ainv_[0][1] = -dy_dxsi1 / djacpp;
    ainv_[1][0] = -dx_dxsi2 / djacpp;
    ainv_[1][1] = dx_dxsi1 / djacpp;

    // Computing the quadrature jacobian matrix - dx/dxqsi
    int degm = (*iparametersC)[NpatchC_]->getDegree(0);
    int degn = (*iparametersC)[NpatchC_]->getDegree(1);
    int npcm = (*iparametersC)[NpatchC_]->getNcp(0);
    int npcn = (*iparametersC)[NpatchC_]->getNcp(1);
    double *uknot_ = (*iparametersC)[NpatchC_]->getuKnot();
    double *vknot_ = (*iparametersC)[NpatchC_]->getvKnot();

    int uind = inc_[0];
    int vind = inc_[1];

    double u1 = uknot_[uind];
    double u2 = uknot_[uind + 1];
    double v1 = vknot_[vind];
    double v2 = vknot_[vind + 1];

    double dxsi1_dqxsi1 = (u2 - u1) * 0.5;
    double dxsi1_dqxsi2 = 0.;
    double dxsi2_dqxsi1 = 0.;
    double dxsi2_dqxsi2 = (v2 - v1) * 0.5;

    quadJacMat[0][0] = dxsi1_dqxsi1 * dx_dxsi1 + dxsi2_dqxsi1 * dx_dxsi2;
    quadJacMat[0][1] = dxsi1_dqxsi2 * dx_dxsi1 + dxsi2_dqxsi2 * dx_dxsi2;
    quadJacMat[1][0] = dxsi1_dqxsi1 * dy_dxsi1 + dxsi2_dqxsi1 * dy_dxsi2;
    quadJacMat[1][1] = dxsi1_dqxsi2 * dy_dxsi1 + dxsi2_dqxsi2 * dy_dxsi2;

    for (int i = 0; i < 2; ++i) delete[] dphi[i];
    delete[] dphi;

    return;
};

template <>
void Element<2>::getQuadJacobianMatrix_ISO(double *xsi, double **quadJacMatInv)
{

    QuadShapeFunction<2> shapeQuad;
    double &alpha_f = parameters->getAlphaF();

    // data from iga element
    int *inc_ = (*nodes_)[connect_[8]]->getINC();
    double wpc[9] = {};
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();

    double **dphi;
    dphi = new double *[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparameters), Npatch_);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    // Computing Jacobian matrix (physical to parametric)- dx/dxsi
    for (int i = 0; i < 9; i++)
    {

        double xna_ = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(0) +
                      (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]]->getCoordinateValue(1) +
                      (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousCoordinateValue(1);

        dx_dxsi1 += xna_ * dphi[0][i] / wpc[i];
        dx_dxsi2 += xna_ * dphi[1][i] / wpc[i];
        dy_dxsi1 += yna_ * dphi[0][i] / wpc[i];
        dy_dxsi2 += yna_ * dphi[1][i] / wpc[i];
    };

    // Computing the Jacobian matrix (physical to quadrature) - dx/dxqsi
    int degm = (*iparameters)[Npatch_]->getDegree(0);
    int degn = (*iparameters)[Npatch_]->getDegree(1);
    int npcm = (*iparameters)[Npatch_]->getNcp(0);
    int npcn = (*iparameters)[Npatch_]->getNcp(1);

    double *uknot_ = (*iparameters)[Npatch_]->getuKnot();
    double *vknot_ = (*iparameters)[Npatch_]->getvKnot();

    double u1 = uknot_[inc_[0]];
    double u2 = uknot_[inc_[0] + 1];
    double v1 = vknot_[inc_[1]];
    double v2 = vknot_[inc_[1] + 1];

    double dxsi1_dqxsi1 = (u2 - u1) * 0.5;
    double dxsi1_dqxsi2 = 0.;
    double dxsi2_dqxsi1 = 0.;
    double dxsi2_dqxsi2 = (v2 - v1) * 0.5;

    double quadJacMat[2][2];

    quadJacMat[0][0] = dxsi1_dqxsi1 * dx_dxsi1 + dxsi2_dqxsi1 * dx_dxsi2;
    quadJacMat[0][1] = dxsi1_dqxsi2 * dx_dxsi1 + dxsi2_dqxsi2 * dx_dxsi2;
    quadJacMat[1][0] = dxsi1_dqxsi1 * dy_dxsi1 + dxsi2_dqxsi1 * dy_dxsi2;
    quadJacMat[1][1] = dxsi1_dqxsi2 * dy_dxsi1 + dxsi2_dqxsi2 * dy_dxsi2;

    double djac_ = quadJacMat[0][0] * quadJacMat[1][1] - quadJacMat[0][1] * quadJacMat[1][0];

    // Computing inverse Jacobian (quadrature to physical) - dxqsi/dx
    quadJacMatInv[0][0] = quadJacMat[1][1] / djac_;
    quadJacMatInv[0][1] = -quadJacMat[0][1] / djac_;
    quadJacMatInv[1][0] = -quadJacMat[1][0] / djac_;
    quadJacMatInv[1][1] = quadJacMat[0][0] / djac_;

    for (int i = 0; i < 2; ++i) delete[] dphi[i];
    delete[] dphi;
}

//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template <>
void Element<2>::getSpatialDerivatives_FEM(double *xsi, double **ainv_, double **dphi_dx)
{

    QuadShapeFunction<2> shapeQuad;

    double **dphi;
    dphi = new double *[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[6];

    shapeQuad.evaluateGradientFem(xsi, dphi);

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 6; ++j)
            dphi_dx[i][j] = 0.;

    // Quadratic shape functions spatial first derivatives
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 6; k++)
                dphi_dx[i][k] += ainv_[i][j] * dphi[j][k];

    for (int i = 0; i < 2; ++i) delete[] dphi[i];
    delete[] dphi;

    return;
};

template <>
void Element<2>::getSpatialDerivatives_ISO(double *xsi, double **ainv_, double **dphi_dx)
{

    QuadShapeFunction<2> shapeQuad;

    double **dphi;
    dphi = new double *[2];
    for (int i = 0; i < 2; ++i)
        dphi[i] = new double[9];

    double wpc[9] = {};
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparameters), Npatch_);

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 9; ++j)
            dphi_dx[i][j] = 0.;

    // Quadratic shape functions spatial first derivatives
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 9; k++)
                dphi_dx[i][k] += ainv_[i][j] * dphi[j][k];

    for (int i = 0; i < 2; ++i) delete[] dphi[i];
    delete[] dphi;
    return;
};

template <>
void Element<2>::getSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double **dphi_dx)
{

    QuadShapeFunction<2> shapeQuad;

    double **dphi;
    dphi = new double *[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];

    double wpc[9] = {};
    for (int i = 0; i < 9; i++) wpc[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *inc_ = (*nodesC_)[connectC_[8]]->getINC();

    shapeQuad.evaluateGradientIso(xsi, dphi, wpc, inc_, (*iparametersC), NpatchC_);

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 9; ++j)
            dphi_dx[i][j] = 0.;

    // Quadratic shape functions spatial first derivatives
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 9; k++)
                dphi_dx[i][k] += ainv_[i][j] * dphi[j][k];

    for (int i = 0; i < 2; ++i) delete[] dphi[i];
    delete[] dphi;
    return;
};

template <>
void Element<2>::getSecondSpatialDerivatives_FEM(double **ainv_, double ***ddphi_dx)
{

    QuadShapeFunction<2> shapeQuad;
    double inter[2][2][6] = {};
    double ainvT_[2][2] = {};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            ainvT_[i][j] = ainv_[j][i];
        };
    };

    double ***ddphi;
    ddphi = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddphi[i] = new double *[2];
        for (int j = 0; j < 2; j++)
            ddphi[i][j] = new double[6];
    }

    shapeQuad.evaluateHessianFem(ddphi);

    // Quadratic shape functions spatial second derivatives
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 6; nf++)
                    inter[i][j][nf] += ainv_[i][k] * ddphi[k][j][nf];

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 6; ++k)
                ddphi_dx[i][j][k] = 0.;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 6; nf++)
                    ddphi_dx[i][j][nf] += inter[i][k][nf] * ainvT_[k][j];

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] ddphi[i][j];
        };
        delete[] ddphi[i];
    };
    delete[] ddphi;

    return;
};

template <>
void Element<2>::getSecondSpatialDerivatives_ISO(double *xsi, double **ainv_, double ***ddphi_dx)
{

    QuadShapeFunction<2> shapeQuad;
    double ainvT_[2][2] = {};
    double inter[2][2][9] = {};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            ainvT_[i][j] = ainv_[j][i];
        };
    };

    double ***ddphi;
    ddphi = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddphi[i] = new double *[2];
        for (int j = 0; j < 2; j++)
            ddphi[i][j] = new double[9];
    }

    double wpc[9] = {};
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    shapeQuad.evaluateHessianIso(xsi, ddphi, wpc, inc_, (*iparameters), Npatch_);

    // Quadratic shape functions spatial second derivatives
    
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 9; nf++)
                    inter[i][j][nf] += ainv_[i][k] * ddphi[k][j][nf];

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 9; ++k)
                ddphi_dx[i][j][k] = 0.;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 9; nf++)
                    ddphi_dx[i][j][nf] += inter[i][k][nf] * ainvT_[k][j];

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] ddphi[i][j];
        };
        delete[] ddphi[i];
    };
    delete[] ddphi;

    return;
};

template <>
void Element<2>::getSecondSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double ***ddphi_dx)
{

    QuadShapeFunction<2> shapeQuad;

    double ainvT_[2][2] = {};
    double inter[2][2][9] = {};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            ainvT_[i][j] = ainv_[j][i];
        };
    };

    double ***ddphi;
    ddphi = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddphi[i] = new double *[2];
        for (int j = 0; j < 2; j++)
            ddphi[i][j] = new double[9];
    }

    double wpc[9] = {};
    for (int i = 0; i < 9; i++)
        wpc[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *inc_ = (*nodesC_)[connectC_[8]]->getINC();

    shapeQuad.evaluateHessianIso(xsi, ddphi, wpc, inc_, (*iparametersC), NpatchC_);

    // Quadratic shape functions spatial second derivatives
    
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 9; nf++)
                    inter[i][j][nf] += ainv_[i][k] * ddphi[k][j][nf];

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 9; ++k)
                ddphi_dx[i][j][k] = 0.;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 9; nf++)
                    ddphi_dx[i][j][nf] += inter[i][k][nf] * ainvT_[k][j];

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
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
void Element<DIM>::getInterpCoord(int &na, double *phi_, double *x_, double *xPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        x_[i] = 0.;
        xPrev_[i] = 0.;
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            x_[j] += (*nodes_)[connect_[i]]->getCoordinateValue(j) * phi_[i];
            xPrev_[j] += (*nodes_)[connect_[i]]->getPreviousCoordinateValue(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpCoord_ISO(int &na, double *phi_, double *x_, double *xPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        x_[i] = 0.;
        xPrev_[i] = 0.;
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            x_[j] += (*nodes_)[connect_[i]]->getCoordinateValue(j) * phi_[i] / (*nodes_)[connect_[i]] -> getWeightPC();
            xPrev_[j] += (*nodes_)[connect_[i]]->getPreviousCoordinateValue(j) * phi_[i] / (*nodes_)[connect_[i]] -> getWeightPC();
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpVel(int &na, double *phi_, double *u_, double *uPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        u_[i] = 0.;
        uPrev_[i] = 0.;
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            u_[j] += (*nodes_)[connect_[i]]->getVelocity(j) * phi_[i];
            uPrev_[j] += (*nodes_)[connect_[i]]->getPreviousVelocity(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpVelDer(int &na, double **dphi_dx, double **du_dx, double **duPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            du_dx[i][j] = 0.;
            duPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
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
void Element<DIM>::getInterpSecondVelDer(int &na, double ***ddphi_dx, double ***ddu_dxdx, double ***dduPrev_dxdx)
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

    for (int i = 0; i < na; i++)
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
void Element<DIM>::getInterpVelCoarse(int &na, double *phi_, double *u_, double *uPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        u_[i] = 0.;
        uPrev_[i] = 0.;
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            u_[j] += (*nodesC_)[connectC_[i]]->getVelocity(j) * phi_[i];
            uPrev_[j] += (*nodesC_)[connectC_[i]]->getPreviousVelocity(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpVelDerCoarse(int &na, double **dphi_dx, double **du_dx, double **duPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            du_dx[i][j] = 0.;
            duPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
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
void Element<DIM>::getInterpSecondVelDerCoarse(int &na, double ***ddphi_dx, double ***ddu_dxdx, double ***dduPrev_dxdx)
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

    for (int i = 0; i < na; i++)
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
void Element<DIM>::getInterpMeshVel(int &na, double *phi_, double *uMesh_, double *uMeshPrev_)
{

   for (int i = 0; i < DIM; i++)
    {
        uMesh_[i] = 0.;
        uMeshPrev_[i] = 0.;
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            uMesh_[j] += (*nodes_)[connect_[i]]->getMeshVelocity(j) * phi_[i];
            uMeshPrev_[j] += (*nodes_)[connect_[i]]->getPreviousMeshVelocity(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpMeshVelDer(int &na, double **dphi_dx, double **duMesh_dx, double **duMeshPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            duMesh_dx[i][j] = 0.;
            duMeshPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
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
void Element<DIM>::getInterpMeshVelCoarse(int &na, double *phi_, double *uMesh_, double *uMeshPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        uMesh_[i] = 0.;
        uMeshPrev_[i] = 0.;
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            uMesh_[j] += (*nodesC_)[connectC_[i]]->getMeshVelocity(j) * phi_[i];
            uMeshPrev_[j] += (*nodesC_)[connectC_[i]]->getPreviousMeshVelocity(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpMeshVelDerCoarse(int &na, double **dphi_dx, double **duMesh_dx, double **duMeshPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            duMesh_dx[i][j] = 0.;
            duMeshPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
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
void Element<DIM>::getInterpAccel(int &na, double *phi_, double *accel_, double *accelPrev_)
{

    for (int i = 0; i < DIM; i++)
    {
        accel_[i] = 0.;
        accelPrev_[i] = 0.;
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            accel_[j] += (*nodes_)[connect_[i]]->getAcceleration(j) * phi_[i];
            accelPrev_[j] += (*nodes_)[connect_[i]]->getPreviousAcceleration(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpAccelDer(int &na, double **dphi_dx, double **daccel_dx, double **daccelPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            daccel_dx[i][j] = 0.;
            daccelPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
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
void Element<DIM>::getInterpAccelDerCoarse(int &na, double **dphi_dx, double **daccel_dx, double **daccelPrev_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            daccel_dx[i][j] = 0.;
            daccelPrev_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
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
void Element<DIM>::getInterpPress(int &na, double *phi_, double &press_)
{
    press_ = 0.;

    for (int i = 0; i < na; i++)
    {
        press_ += (*nodes_)[connect_[i]]->getPressure() * phi_[i];
    };
};

template <int DIM>
void Element<DIM>::getInterpPressDer(int &na, double **dphi_dx, double *dpress_dx)
{

    for (int i = 0; i < DIM; i++)
    {
        dpress_dx[i] = 0.;
   
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            dpress_dx[j] += (*nodes_)[connect_[i]]->getPressure() * dphi_dx[j][i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpSecondPressDer(int &na, double ***ddphi_dx, double **ddpress_dxdx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            ddpress_dxdx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                ddpress_dxdx[j][k] += (*nodes_)[connect_[i]]->getPressure() * ddphi_dx[j][k][i];
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpSecondPressDerCoarse(int &na, double ***ddphi_dx, double **ddpress_dxdx)
{

    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            ddpress_dxdx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                ddpress_dxdx[j][k] += (*nodesC_)[connectC_[i]]->getPressure() * ddphi_dx[j][k][i];
            };
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpLambda(int &na, double *phi_, double *lambda_)
{

    for (int i = 0; i < DIM; i++)
    {
        lambda_[i] = 0.;

    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            lambda_[j] += (*nodes_)[connect_[i]]->getLagrangeMultiplier(j) * phi_[i];
        };
    };
};

template <int DIM>
void Element<DIM>::getInterpLambdaDer(int &na, double **dphi_dx, double **dlambda_dx)
{

   for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            dlambda_dx[i][j] = 0.;
        };
    };

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                dlambda_dx[j][k] += (*nodes_)[connect_[i]]->getLagrangeMultiplier(j) * dphi_dx[k][i];
            };
        };
    };
};

//------------------------------------------------------------------------------
//-----------------------COMPUTE DRAG AND LIFT FORCES --------------------------
//------------------------------------------------------------------------------
template <>
void Element<2>::computeDragAndLiftForces_ISO(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce,
                                              double &dForce, double &lForce, double &aux_Mom, double &aux_Per)
{

    double localNodesBoundaryIso_[3][2];
    double pcWeightl[3];
    double pcWeight[9];
    int *incL_;

    for (int i = 0; i < 9; i++) pcWeight[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    if (side == 0)
    {

        pcWeightl[0] = (*nodes_)[connect_[0]]->getWeightPC();
        pcWeightl[1] = (*nodes_)[connect_[1]]->getWeightPC();
        pcWeightl[2] = (*nodes_)[connect_[2]]->getWeightPC();

        for (int i = 0; i < 2; i++)
        {
            localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[0]]->getCoordinateValue(i)) / pcWeightl[0];
            localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[1]]->getCoordinateValue(i)) / pcWeightl[1];
            localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[2]]->getCoordinateValue(i)) / pcWeightl[2];
        };

        incL_ = (*nodes_)[connect_[2]]->getINC();
    }
    else
    {

        if (side == 1)
        {

            pcWeightl[0] = (*nodes_)[connect_[2]]->getWeightPC();
            pcWeightl[1] = (*nodes_)[connect_[5]]->getWeightPC();
            pcWeightl[2] = (*nodes_)[connect_[8]]->getWeightPC();

            for (int i = 0; i < 2; i++)
            {
                localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[2]]->getCoordinateValue(i)) / pcWeightl[0];
                localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[5]]->getCoordinateValue(i)) / pcWeightl[1];
                localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[8]]->getCoordinateValue(i)) / pcWeightl[2];
            };

            incL_ = (*nodes_)[connect_[8]]->getINC();

            if (side == 2)
            {

                pcWeightl[0] = (*nodes_)[connect_[8]]->getWeightPC();
                pcWeightl[1] = (*nodes_)[connect_[7]]->getWeightPC();
                pcWeightl[2] = (*nodes_)[connect_[6]]->getWeightPC();

                for (int i = 0; i < 2; i++)
                {
                    localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[8]]->getCoordinateValue(i)) / pcWeightl[0];
                    localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[7]]->getCoordinateValue(i)) / pcWeightl[1];
                    localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[6]]->getCoordinateValue(i)) / pcWeightl[2];
                };

                incL_ = (*nodes_)[connect_[8]]->getINC();
            };
        }
        else
        {

            pcWeightl[0] = (*nodes_)[connect_[6]]->getWeightPC();
            pcWeightl[1] = (*nodes_)[connect_[3]]->getWeightPC();
            pcWeightl[2] = (*nodes_)[connect_[0]]->getWeightPC();

            for (int i = 0; i < 2; i++)
            {
                localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[6]]->getCoordinateValue(i)) / pcWeightl[0];
                localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[3]]->getCoordinateValue(i)) / pcWeightl[1];
                localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[0]]->getCoordinateValue(i)) / pcWeightl[2];
            };

            incL_ = (*nodes_)[connect_[6]]->getINC();
        };
    };

    BoundaryIntegQuadratureIso<2> bQuad;
    QuadShapeFunction<2> shapeQuad;
    BoundShapeFunction<2> shapeBound;
    double &visc_ = parameters->getViscosity();

    double n_vector[2] = {};
    double shearStress[2][2] = {};
    double load_friction[2] = {};
    double load_pressure[2] = {};

    double ident[2][2];
    ident[0][0] = 1.;
    ident[0][1] = 0.;
    ident[1][0] = 0.;
    ident[1][1] = 1.;

    double phi_[9], phiCB_[3], xsi[2], xsiB[1], weightB;

    double **dphiCB_;
    dphiCB_ = new double *[1];
    for (int i = 0; i < 1; ++i) dphiCB_[i] = new double[3];

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double *[2];
    for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];

    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double moment = 0.;
    double per = 0.;

    int index = 0;
    for (double *it = bQuad.begin(); it != bQuad.end(); it++)
    {

        xsiB[0] = bQuad.PointList(index, 0);
        weightB = bQuad.WeightList(index);

        if (side == 0)
        {
            xsi[0] = xsiB[0];
            xsi[1] = -1;
        };
        if (side == 1)
        {
            xsi[0] = -1;
            xsi[1] = xsiB[0];
        };
        if (side == 2)
        {
            xsi[0] = xsiB[0];
            xsi[1] = 1;
        };
        if (side == 3)
        {
            xsi[0] = -1;
            xsi[1] = xsiB[0];
        };

        // Computes the shape functions
        shapeQuad.evaluateIso(xsi, phi_, pcWeight, inc_, (*iparameters), Npatch_);

        // Computes the jacobian matrix
        double djac_;
        getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

        // Interpolates velocity and its derivatives values
        int na = 9;
        getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);

        double press_;
        getInterpPress(na, phi_, press_);

        // Get Functions and Derivatives on the boundary
        shapeBound.evaluateBoundaryIso(xsiB, phiCB_, pcWeightl, side, incL_, (*iparameters), Npatch_);
        shapeBound.evaluateGradientBoundaryIso(xsiB, dphiCB_, pcWeightl, side, incL_, (*iparameters), Npatch_);

        // (Jacobian parametric space to physical one)
        double Tx = 0.;
        double Ty = 0.;
        for (int i = 0; i < 3; i++)
        {
            Tx += localNodesBoundaryIso_[i][0] * dphiCB_[0][i];
            Ty += localNodesBoundaryIso_[i][1] * dphiCB_[0][i];
        };

        // (Jacobian parental space to the physical one)
        if ((side == 0) || (side == 2))
        {

            int degm = (*iparameters)[Npatch_]->getDegree(0);
            int npcm = (*iparameters)[Npatch_]->getNcp(0);
            double *uknot_ = (*iparameters)[Npatch_]->getuKnot();

            int uind = incL_[0];

            double u1 = uknot_[uind];
            double u2 = uknot_[uind + 1];

            double dxsi1_dqxsi1 = (u2 - u1) * 0.5;

            Tx *= dxsi1_dqxsi1;
            Ty *= dxsi1_dqxsi1;
        }
        else
        {

            int degn = (*iparameters)[Npatch_]->getDegree(1);
            int npcn = (*iparameters)[Npatch_]->getNcp(1);
            double *vknot_ = (*iparameters)[Npatch_]->getvKnot();

            int vind = incL_[1];

            double v1 = vknot_[vind];
            double v2 = vknot_[vind + 1];

            double dxsi2_dqxsi2 = (v2 - v1) * 0.5;

            Tx *= dxsi2_dqxsi2;
            Ty *= dxsi2_dqxsi2;
        }

        double jacb_ = sqrt(Tx * Tx + Ty * Ty);

        n_vector[0] = Ty / jacb_;
        n_vector[1] = -Tx / jacb_;

        shearStress[0][0] = 2. * visc_ * du_dx[0][0];
        shearStress[0][1] = visc_ * (du_dx[0][1] + du_dx[1][0]);
        shearStress[1][0] = visc_ * (du_dx[0][1] + du_dx[1][1]);
        shearStress[1][1] = 2. * visc_ * du_dx[1][1];

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                load_pressure[i] += -press_ * ident[i][j] * n_vector[j] * jacb_ * weightB;
                load_friction[i] += shearStress[i][j] * n_vector[j] * jacb_ * weightB;
            }
        }

        // moment +=  (-(load_pressure[0]+load_friction[0]) * (y_- 4.)
        //            + (load_pressure[0]+load_friction[0]) * (x_ - 4.)) * jacb_ * weightB[0];
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

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;
    for (int i = 0; i < 1; ++i) delete[] dphiCB_[i];
    delete[] dphiCB_;
    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;
    for (int i = 0; i < 2; ++i) delete[] quadJacMat[i];
    delete[] quadJacMat;
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;
    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    return;
};

template <>
void Element<2>::computeDragAndLiftForces_FEM(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce,
                                              double &dForce, double &lForce, double &aux_Mom, double &aux_Per)
{

    double localNodesBoundary_[3][2];

    if (side == 0)
    {
        for (int i = 0; i < 2; i++)
        {
            localNodesBoundary_[0][i] = (*nodes_)[connect_[1]]->getCoordinateValue(i);
            localNodesBoundary_[1][i] = (*nodes_)[connect_[4]]->getCoordinateValue(i);
            localNodesBoundary_[2][i] = (*nodes_)[connect_[2]]->getCoordinateValue(i);
        };
    }
    else
    {
        if (side == 1)
        {
            for (int i = 0; i < 2; i++)
            {
                localNodesBoundary_[0][i] = (*nodes_)[connect_[2]]->getCoordinateValue(i);
                localNodesBoundary_[1][i] = (*nodes_)[connect_[5]]->getCoordinateValue(i);
                localNodesBoundary_[2][i] = (*nodes_)[connect_[0]]->getCoordinateValue(i);
            };
        }
        else
        {
            for (int i = 0; i < 2; i++)
            {
                localNodesBoundary_[0][i] = (*nodes_)[connect_[0]]->getCoordinateValue(i);
                localNodesBoundary_[1][i] = (*nodes_)[connect_[3]]->getCoordinateValue(i);
                localNodesBoundary_[2][i] = (*nodes_)[connect_[1]]->getCoordinateValue(i);
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
    ident[0][1] = 0.;
    ident[1][0] = 0.;
    ident[1][1] = 1.;

    double phi_[6], phiB_[3], xsi[2], xsiB[1], weightB;

    double **dphiB_;
    dphiB_ = new double *[1];
    for (int i = 0; i < 1; ++i) dphiB_[i] = new double[3];

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double *[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

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

        int na = 6;
        getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);

        double press_;
        getInterpPress(na, phi_, press_);

        shapeBound.evaluateBoundaryFem(xsiB, phiB_);
        shapeBound.evaluateGradientBoundaryFem(xsiB, dphiB_);

        // computes the normal vector in the linear element
        double Tx = 0.;
        double Ty = 0.;
        for (int i = 0; i < 3; i++)
        {
            Tx += localNodesBoundary_[i][0] * dphiB_[0][i];
            Ty += localNodesBoundary_[i][1] * dphiB_[0][i];
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

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 1; ++i) delete[] dphiB_[i];
    delete[] dphiB_;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;
    
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;
    
    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    return;
};

//------------------------------------------------------------------------------
//--------------------COMPUTE BEZIER TRANSFORMATION MATRIX----------------------
//------------------------------------------------------------------------------
template <>
void Element<2>::getMatrixC(double **MatrixC)
{

    int degm = (*iparameters)[Npatch_]->getDegree(0); // Degree in u direction
    int degn = (*iparameters)[Npatch_]->getDegree(1); // Degree in v direction
    int npcm = (*iparameters)[Npatch_]->getNcp(0);    // Number of control points in u direction
    int npcn = (*iparameters)[Npatch_]->getNcp(1);    // Number of control points in v direction
    double *uknot_ = (*iparameters)[Npatch_]->getuKnot();
    double *vknot_ = (*iparameters)[Npatch_]->getvKnot(); // Knots v direction
    int dimu = degm + npcm + 1;
    int dimv = degn + npcn + 1;
    double MatrixCu[3][3] = {};
    double MatrixCuNext[3][3] = {};
    double MatrixCv[3][3] = {};
    double MatrixCvNext[3][3] = {};
    double alphasu[3] = {};
    double alphasv[3] = {};

    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    int uind = inc_[0];
    int vind = inc_[1];

    // BEZIER EXTRACTION C U DIRECTION
    int au_ = degm;
    int bu_ = au_ + 1;
    int nb = 0;               // number of Bézier Elements
    int i_, mult, r, s, save; // auxiliar
    double numer, alpha;

    for (int i = 0; i <= degm; i++)
    {
        MatrixCu[i][i] = 1.;
        MatrixCuNext[i][i] = 1.;
    };
    while (bu_ <= uind + 1)
    {
        for (int i = 0; i <= degm; i++)
        {
            for (int j = 0; j <= degm; j++)
            {
                MatrixCu[i][j] = MatrixCuNext[i][j];
            };
        };
        for (int i = 0; i <= degm; i++)
        {
            MatrixCuNext[i][i] = 1.;
        };
        i_ = bu_;
        while ((bu_ < dimu - 1) && (uknot_[bu_ + 1] == uknot_[bu_]))
            bu_++;
        mult = bu_ - i_ + 1;
        if (mult < degm)
        {
            numer = uknot_[bu_] - uknot_[au_];
            for (int j = degm; j > mult; j--)
            {
                alphasu[j - mult - 1] = numer / (uknot_[au_ + j] - uknot_[au_]);
            };
            r = degm - mult; // Insert knot r times
            for (int j = 1; j <= r; j++)
            {
                save = r - j;
                s = mult + j;
                for (int k = degm; k >= s; k--)
                {
                    alpha = alphasu[k - s];
                    for (int i = 0; i <= degm; i++)
                    {
                        MatrixCu[i][k] = alpha * MatrixCu[i][k] + (1.0 - alpha) * MatrixCu[i][k - 1];
                    };
                };
                int cc = degm - j;
                for (int i = save; i <= j + save; i++)
                {
                    MatrixCuNext[i][save] = MatrixCu[cc][degm];
                    cc++;
                };
            };
            nb = nb + 1;
            if (bu_ <= uind + 1)
            {
                au_ = bu_;
                bu_++;
            };
        }
        else
        {
            nb = nb + 1;
            if (bu_ <= uind + 1)
            {
                au_ = bu_;
            };
        };
    };

    // BEZIER EXTRACTION C V DIRECTION
    au_ = degn;
    bu_ = au_ + 1;
    nb = 0; // number of Bézier Elements
    for (int i = 0; i <= degn; i++)
    {
        MatrixCv[i][i] = 1.;
        MatrixCvNext[i][i] = 1.;
    };
    while (bu_ <= vind + 1)
    {
        for (int i = 0; i <= degn; i++)
        {
            for (int j = 0; j <= degn; j++)
            {
                MatrixCv[i][j] = MatrixCvNext[i][j];
            };
        };
        for (int i = 0; i <= degn; i++)
        {
            MatrixCvNext[i][i] = 1.;
        };
        i_ = bu_;
        while ((bu_ < dimv - 1) && (vknot_[bu_ + 1] == vknot_[bu_]))
            bu_++;
        mult = bu_ - i_ + 1;
        if (mult < degn)
        {
            numer = vknot_[bu_] - vknot_[au_];
            for (int j = degn; j > mult; j--)
            {
                alphasv[j - mult - 1] = numer / (vknot_[au_ + j] - vknot_[au_]);
            };
            r = degn - mult; // Insert knot r times
            for (int j = 1; j <= r; j++)
            {
                save = r - j;
                s = mult + j;
                for (int k = degn; k >= s; k--)
                {
                    alpha = alphasv[k - s];
                    for (int i = 0; i <= degn; i++)
                    {
                        MatrixCv[i][k] = alpha * MatrixCv[i][k] + (1.0 - alpha) * MatrixCv[i][k - 1];
                    };
                };
                int cc = degn - j;
                for (int i = save; i <= j + save; i++)
                {
                    MatrixCvNext[i][save] = MatrixCv[cc][degm];
                    cc++;
                };
            };
            nb = nb + 1;
            if (bu_ <= vind + 1)
            {
                au_ = bu_;
                bu_++;
            };
        }
        else
        {
            nb = nb + 1;
            if (bu_ <= uind + 1)
            {
                au_ = bu_;
            };
        };
    };

    for (int k = 0; k <= degn; k++)
    {
        for (int l = 0; l <= degn; l++)
        {
            for (int i = 0; i <= degm; i++)
            {
                for (int j = 0; j <= degm; j++)
                {
                    MatrixC[i + k * (degn + 1)][j + l * (degm + 1)] = MatrixCu[i][j] * MatrixCv[k][l];
                };
            };
        };
    };
    return;
}

template <>
void Element<2>::getInvMatrixC(double **MatrixCuInv, double **MatrixCvInv)
{

    // BEZIER EXTRACTION MATRIXC - PARAMETRIC DIRECTION u
    int degm = (*iparameters)[Npatch_]->getDegree(0); // Degree
    int npcm = (*iparameters)[Npatch_]->getNcp(0);    // Number of control points
    int dimu = degm + npcm + 1;                       // Number of Knots
    double *uknot_ = (*iparameters)[Npatch_]->getuKnot();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();
    int uind = inc_[0]; // index of base fuction (parametric direction u) that starts in the index_ element
    double alphasu[3];
    double MatrixCu[3][3] = {};
    double MatrixCuNext[3][3] = {};

    // Algorithm
    int au_ = degm;
    int bu_ = au_ + 1;
    int nb = 0;               // number of Bézier Elements
    int i_, mult, r, s, save; // auxiliar
    double numer, alpha;

    for (int i = 0; i <= degm; i++)
    {
        MatrixCu[i][i] = 1.;
        MatrixCuNext[i][i] = 1.;
    };
    while (bu_ <= uind + 1)
    {
        for (int i = 0; i <= degm; i++)
        {
            for (int j = 0; j <= degm; j++)
            {
                MatrixCu[i][j] = MatrixCuNext[i][j];
            };
        };
        for (int i = 0; i <= degm; i++)
        {
            MatrixCuNext[i][i] = 1.;
        };
        i_ = bu_;
        while ((bu_ < dimu - 1) && (uknot_[bu_ + 1] == uknot_[bu_]))
            bu_++;
        mult = bu_ - i_ + 1;
        if (mult < degm)
        {
            numer = uknot_[bu_] - uknot_[au_];
            for (int j = degm; j > mult; j--)
            {
                alphasu[j - mult - 1] = numer / (uknot_[au_ + j] - uknot_[au_]);
            };
            r = degm - mult; // Insert knot r times
            for (int j = 1; j <= r; j++)
            {
                save = r - j;
                s = mult + j;
                for (int k = degm; k >= s; k--)
                {
                    alpha = alphasu[k - s];
                    for (int i = 0; i <= degm; i++)
                    {
                        MatrixCu[i][k] = alpha * MatrixCu[i][k] + (1.0 - alpha) * MatrixCu[i][k - 1];
                    };
                };
                int cc = degm - j;
                for (int i = save; i <= j + save; i++)
                {
                    MatrixCuNext[i][save] = MatrixCu[cc][degm];
                    cc++;
                };
            };
            nb = nb + 1;
            if (bu_ <= uind + 1)
            {

                au_ = bu_;
                bu_++;
            };
        }
        else
        {
            nb = nb + 1;
            if (bu_ <= uind + 1)
            {
                au_ = bu_;
            };
        };
    };

    // Computing Matrix Cu determinant
    double det_ = MatrixCu[0][0] * MatrixCu[1][1] * MatrixCu[2][2] + MatrixCu[0][1] * MatrixCu[1][2] * MatrixCu[2][0] + MatrixCu[0][2] * MatrixCu[1][0] * MatrixCu[2][1] -
                  MatrixCu[0][1] * MatrixCu[1][0] * MatrixCu[2][2] - MatrixCu[0][0] * MatrixCu[1][2] * MatrixCu[2][1] - MatrixCu[0][2] * MatrixCu[1][1] * MatrixCu[2][0];

    // Inverting Matrix Cu
    MatrixCuInv[0][0] = (MatrixCu[1][1] * MatrixCu[2][2] - MatrixCu[1][2] * MatrixCu[2][1]) / det_;
    MatrixCuInv[1][0] = -(MatrixCu[1][0] * MatrixCu[2][2] - MatrixCu[1][2] * MatrixCu[2][0]) / det_;
    MatrixCuInv[2][0] = (MatrixCu[1][0] * MatrixCu[2][1] - MatrixCu[1][1] * MatrixCu[2][0]) / det_;
    MatrixCuInv[0][1] = -(MatrixCu[0][1] * MatrixCu[2][2] - MatrixCu[0][2] * MatrixCu[2][1]) / det_;
    MatrixCuInv[1][1] = (MatrixCu[0][0] * MatrixCu[2][2] - MatrixCu[0][2] * MatrixCu[2][0]) / det_;
    MatrixCuInv[2][1] = -(MatrixCu[0][0] * MatrixCu[2][1] - MatrixCu[0][1] * MatrixCu[2][0]) / det_;
    MatrixCuInv[0][2] = (MatrixCu[0][1] * MatrixCu[1][2] - MatrixCu[0][2] * MatrixCu[1][1]) / det_;
    MatrixCuInv[1][2] = -(MatrixCu[0][0] * MatrixCu[1][2] - MatrixCu[0][2] * MatrixCu[1][0]) / det_;
    MatrixCuInv[2][2] = (MatrixCu[0][0] * MatrixCu[1][1] - MatrixCu[0][1] * MatrixCu[1][0]) / det_;

    // BEZIER EXTRACTION MATRIXC - PARAMETRIC DIRECTION v
    int degn = (*iparameters)[Npatch_]->getDegree(1); // Degree
    int npcn = (*iparameters)[Npatch_]->getNcp(1);    // Number of control points
    int dimv = degn + npcn + 1;
    double *vknot_ = (*iparameters)[Npatch_]->getvKnot();
    int vind = inc_[1]; // index of base fuction (parametric direction v) that starts in the index_ element
    double alphasv[3] = {};
    double MatrixCv[3][3] = {};
    double MatrixCvNext[3][3] = {};

    // BEZIER EXTRACTION C V DIRECTION
    au_ = degn;
    bu_ = au_ + 1;
    nb = 0; // number of Bézier Elements
    for (int i = 0; i <= degn; i++)
    {
        MatrixCv[i][i] = 1.;
        MatrixCvNext[i][i] = 1.;
    };
    while (bu_ <= vind + 1)
    {
        for (int i = 0; i <= degn; i++)
        {
            for (int j = 0; j <= degn; j++)
            {
                MatrixCv[i][j] = MatrixCvNext[i][j];
            };
        };
        for (int i = 0; i <= degn; i++)
        {
            MatrixCvNext[i][i] = 1.;
        };
        i_ = bu_;
        while ((bu_ < dimv - 1) && (vknot_[bu_ + 1] == vknot_[bu_]))
            bu_++;
        mult = bu_ - i_ + 1;
        if (mult < degn)
        {
            numer = vknot_[bu_] - vknot_[au_];
            for (int j = degn; j > mult; j--)
            {
                alphasv[j - mult - 1] = numer / (vknot_[au_ + j] - vknot_[au_]);
            };
            r = degn - mult; // Insert knot r times
            for (int j = 1; j <= r; j++)
            {
                save = r - j;
                s = mult + j;
                for (int k = degn; k >= s; k--)
                {
                    alpha = alphasv[k - s];
                    for (int i = 0; i <= degn; i++)
                    {
                        MatrixCv[i][k] = alpha * MatrixCv[i][k] + (1.0 - alpha) * MatrixCv[i][k - 1];
                    };
                };
                int cc = degn - j;
                for (int i = save; i <= j + save; i++)
                {
                    MatrixCvNext[i][save] = MatrixCv[cc][degm];
                    cc++;
                };
            };
            nb = nb + 1;
            if (bu_ <= vind + 1)
            {
                au_ = bu_;
                bu_++;
            };
        }
        else
        {
            nb = nb + 1;
            if (bu_ <= uind + 1)
            {
                au_ = bu_;
            };
        };
    };

    // Computing Matrix Cv determinant
    det_ = MatrixCv[0][0] * MatrixCv[1][1] * MatrixCv[2][2] + MatrixCv[0][1] * MatrixCv[1][2] * MatrixCv[2][0] + MatrixCv[0][2] * MatrixCv[1][0] * MatrixCv[2][1] -
           MatrixCv[0][1] * MatrixCv[1][0] * MatrixCv[2][2] - MatrixCv[0][0] * MatrixCv[1][2] * MatrixCv[2][1] - MatrixCv[0][2] * MatrixCv[1][1] * MatrixCv[2][0];

    // Inverting Matrix Cv
    MatrixCvInv[0][0] = (MatrixCv[1][1] * MatrixCv[2][2] - MatrixCv[1][2] * MatrixCv[2][1]) / det_;
    MatrixCvInv[1][0] = -(MatrixCv[1][0] * MatrixCv[2][2] - MatrixCv[1][2] * MatrixCv[2][0]) / det_;
    MatrixCvInv[2][0] = (MatrixCv[1][0] * MatrixCv[2][1] - MatrixCv[1][1] * MatrixCv[2][0]) / det_;
    MatrixCvInv[0][1] = -(MatrixCv[0][1] * MatrixCv[2][2] - MatrixCv[0][2] * MatrixCv[2][1]) / det_;
    MatrixCvInv[1][1] = (MatrixCv[0][0] * MatrixCv[2][2] - MatrixCv[0][2] * MatrixCv[2][0]) / det_;
    MatrixCvInv[2][1] = -(MatrixCv[0][0] * MatrixCv[2][1] - MatrixCv[0][1] * MatrixCv[2][0]) / det_;
    MatrixCvInv[0][2] = (MatrixCv[0][1] * MatrixCv[1][2] - MatrixCv[0][2] * MatrixCv[1][1]) / det_;
    MatrixCvInv[1][2] = -(MatrixCv[0][0] * MatrixCv[1][2] - MatrixCv[0][2] * MatrixCv[1][0]) / det_;
    MatrixCvInv[2][2] = (MatrixCv[0][0] * MatrixCv[1][1] - MatrixCv[0][1] * MatrixCv[1][0]) / det_;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template <>
void Element<2>::getNewParameterSUPG_FEM(double &tSUPG_, double &tPSPG_, double &tLSIC_, double **Jac, double *phi_, double **dphi_dx)
{

    double MatrixD[2][2], MatrixInvD[2][2], MatrixQh[2][2], MatrixInvQh[2][2], MatrixG[2][2], rr[2][2];
    double r[2] = {};
    double ua_[2] = {};
    double hrqd, tSUGN1_, tSUGN2_, tSUGN3_;

    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();

    // Matrix D = polynomial order of the base functions
    MatrixD[0][0] = 2.;
    MatrixD[1][1] = 2.;
    MatrixD[1][0] = 0.;
    MatrixD[0][1] = 0.;

    // Inverse Matrix D
    double detMatrixD = MatrixD[0][0] * MatrixD[1][1] - MatrixD[0][1] * MatrixD[1][0];
    MatrixInvD[0][0] = (1. / detMatrixD) * MatrixD[1][1];
    MatrixInvD[1][1] = (1. / detMatrixD) * MatrixD[0][0];
    MatrixInvD[0][1] = -(1. / detMatrixD) * MatrixD[0][1];
    MatrixInvD[1][0] = -(1. / detMatrixD) * MatrixD[1][0];

    // Matrix Q "hat"
    MatrixQh[0][0] = Jac[0][0] * MatrixInvD[0][0] + Jac[0][1] * MatrixInvD[1][0];
    MatrixQh[0][1] = Jac[0][0] * MatrixInvD[0][1] + Jac[0][1] * MatrixInvD[1][1];
    MatrixQh[1][0] = Jac[1][0] * MatrixInvD[0][0] + Jac[1][1] * MatrixInvD[1][0];
    MatrixQh[1][1] = Jac[1][0] * MatrixInvD[0][1] + Jac[1][1] * MatrixInvD[1][1];

    // Matrix Inverse Q "hat"
    double detQh = MatrixQh[0][0] * MatrixQh[1][1] - MatrixQh[0][1] * MatrixQh[1][0];

    MatrixInvQh[0][0] = (1. / detQh) * MatrixQh[1][1];
    MatrixInvQh[1][1] = (1. / detQh) * MatrixQh[0][0];
    MatrixInvQh[0][1] = -(1. / detQh) * MatrixQh[0][1];
    MatrixInvQh[1][0] = -(1. / detQh) * MatrixQh[1][0];

    // Matrix G
    MatrixG[0][0] = MatrixInvQh[0][0] * MatrixInvQh[0][0] + MatrixInvQh[1][0] * MatrixInvQh[1][0];
    MatrixG[0][1] = MatrixInvQh[0][0] * MatrixInvQh[0][1] + MatrixInvQh[1][0] * MatrixInvQh[1][1];
    MatrixG[1][0] = MatrixInvQh[0][1] * MatrixInvQh[0][0] + MatrixInvQh[1][1] * MatrixInvQh[1][0];
    MatrixG[1][1] = MatrixInvQh[0][1] * MatrixInvQh[0][1] + MatrixInvQh[1][1] * MatrixInvQh[1][1];

    // Calculation hrqd
    for (int i = 0; i < 6; i++)
    {

        double ua = alpha_f * (*nodes_)[connect_[i]]->getVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(0);
        double va = alpha_f * (*nodes_)[connect_[i]]->getVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(1);

        double uma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(0);
        double vma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(1);

        ua -= uma;
        va -= vma;

        ua_[0] += ua * phi_[i];
        ua_[1] += va * phi_[i];

        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
    };

    double rNorm = sqrt(r[0] * r[0] + r[1] * r[1]);
    ;
    double referpar = 0.000000000000000001;

    // physical unitary direction from gradient velocity
    r[0] /= (rNorm + referpar);
    r[1] /= (rNorm + referpar);

    rr[0][0] = r[0] * r[0];
    rr[0][1] = r[0] * r[1];
    rr[1][0] = r[1] * r[0];
    rr[1][1] = r[1] * r[1];

    hrqd = 2. / ((sqrt(rr[0][0] * MatrixG[0][0] + rr[0][1] * MatrixG[0][1] + rr[1][0] * MatrixG[1][0] + rr[1][1] * MatrixG[1][1])) + referpar);

    // hmin e hmax
    // calculating the eigen-values from G
    double a_ = 1.;
    double b_ = -(MatrixG[0][0] + MatrixG[1][1]);
    double c_ = (MatrixG[0][0] * MatrixG[1][1] - MatrixG[0][1] * MatrixG[1][0]);

    double lambdaMax = (-b_ + sqrt(b_ * b_ - 4. * a_ * c_)) / (2. * a_);
    double lambdaMin = (-b_ - sqrt(b_ * b_ - 4. * a_ * c_)) / (2. * a_);

    double hmin = 2. / sqrt(lambdaMax);
    double hmax = 2. / sqrt(lambdaMin);

    if (hrqd < hmin)
        hrqd = hmin;
    if (hrqd > hmax)
        hrqd = hmax;

    // tSUGN1_ (-2)
    tSUGN1_ = (ua_[0] * ua_[0]) * MatrixG[0][0] +
              (ua_[0] * ua_[1]) * MatrixG[0][1] +
              (ua_[1] * ua_[0]) * MatrixG[1][0] +
              (ua_[1] * ua_[1]) * MatrixG[1][1];

    tSUGN2_ = dTime_ / 2.;

    // tSUGN3_ (-1)
    tSUGN3_ = (visc_ / dens_) * (4. / (hrqd * hrqd));

    // Computing tSUPG parameter
    tSUPG_ = 1. / (sqrt(tSUGN1_ + (1. / (tSUGN2_ * tSUGN2_)) + (tSUGN3_ * tSUGN3_)));

    tPSPG_ = tSUPG_;

    tLSIC_ = hrqd * hrqd / tSUPG_;

    return;
};

template <>
void Element<2>::getNewParameterSUPG_ISO(double &tSUPG_, double &tPSPG_, double &tLSIC_, double **quadJacMat, double *phi_, double **dphi_dx)
{

    double MatrixD[2][2] = {};
    double MatrixInvD[2][2], MatrixQh[2][2], MatrixInvQh[2][2], MatrixG[2][2], rr[2][2];
    double r[2] = {};
    double ua_[2] = {};
    double hrqd, tSUGN1_, tSUGN2_, tSUGN3_;
    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();

    // Computation of matrix D
    double **MatrixCuInv;
    MatrixCuInv = new double *[3];
    for (int i = 0; i < 3; ++i) MatrixCuInv[i] = new double[3];

    double **MatrixCvInv;
    MatrixCvInv = new double *[3];
    for (int i = 0; i < 3; ++i) MatrixCvInv[i] = new double[3];

    int degm = (*iparameters)[Npatch_]->getDegree(0);
    int degn = (*iparameters)[Npatch_]->getDegree(1);

    // Bezier parametric element lenght
    double dimxsi = 2.;
    double dimeta = 2.;

    getInvMatrixC(MatrixCuInv, MatrixCvInv);

    std::vector<double> Dxsi(2), Deta(2);

    Dxsi.clear();
    Deta.clear();

    for (int i = 0; i <= degm; i++)
    {
        Dxsi[0] += i * (MatrixCuInv[i][1] - MatrixCuInv[i][0]);
        Dxsi[1] += i * (MatrixCuInv[i][2] - MatrixCuInv[i][1]);
    }

    Dxsi[0] *= (dimxsi / degm);
    Dxsi[1] *= (dimxsi / degm);

    for (int i = 0; i <= degn; i++)
    {
        Deta[0] += i * (MatrixCvInv[i][1] - MatrixCvInv[i][0]);
        Deta[1] += i * (MatrixCvInv[i][2] - MatrixCvInv[i][1]);
    }

    Deta[0] *= (dimeta / degn);
    Deta[1] *= (dimeta / degn);

    auto it = std::minmax_element(Dxsi.begin(), Dxsi.end());
    auto it2 = std::minmax_element(Deta.begin(), Deta.end());

    // RQD - MAX
    MatrixD[0][0] = *it.second;
    MatrixD[1][1] = *it2.second;

    // Inverse Matrix D
    double detMatrixD = MatrixD[0][0] * MatrixD[1][1] - MatrixD[0][1] * MatrixD[1][0];
    MatrixInvD[0][0] = (1. / detMatrixD) * MatrixD[1][1];
    MatrixInvD[1][1] = (1. / detMatrixD) * MatrixD[0][0];
    MatrixInvD[0][1] = -(1. / detMatrixD) * MatrixD[0][1];
    MatrixInvD[1][0] = -(1. / detMatrixD) * MatrixD[1][0];

    // Matrix Q "hat"
    MatrixQh[0][0] = quadJacMat[0][0] * MatrixInvD[0][0] + quadJacMat[0][1] * MatrixInvD[1][0];
    MatrixQh[0][1] = quadJacMat[0][0] * MatrixInvD[0][1] + quadJacMat[0][1] * MatrixInvD[1][1];
    MatrixQh[1][0] = quadJacMat[1][0] * MatrixInvD[0][0] + quadJacMat[1][1] * MatrixInvD[1][0];
    MatrixQh[1][1] = quadJacMat[1][0] * MatrixInvD[0][1] + quadJacMat[1][1] * MatrixInvD[1][1];

    // Matrix Inverse Q "hat"
    double detQh = MatrixQh[0][0] * MatrixQh[1][1] - MatrixQh[0][1] * MatrixQh[1][0];

    MatrixInvQh[0][0] = (1. / detQh) * MatrixQh[1][1];
    MatrixInvQh[1][1] = (1. / detQh) * MatrixQh[0][0];
    MatrixInvQh[0][1] = -(1. / detQh) * MatrixQh[0][1];
    MatrixInvQh[1][0] = -(1. / detQh) * MatrixQh[1][0];

    // Matrix G
    MatrixG[0][0] = MatrixInvQh[0][0] * MatrixInvQh[0][0] + MatrixInvQh[1][0] * MatrixInvQh[1][0];
    MatrixG[0][1] = MatrixInvQh[0][0] * MatrixInvQh[0][1] + MatrixInvQh[1][0] * MatrixInvQh[1][1];
    MatrixG[1][0] = MatrixInvQh[0][1] * MatrixInvQh[0][0] + MatrixInvQh[1][1] * MatrixInvQh[1][0];
    MatrixG[1][1] = MatrixInvQh[0][1] * MatrixInvQh[0][1] + MatrixInvQh[1][1] * MatrixInvQh[1][1];

    // Calculation hrqd
    for (int i = 0; i < 9; i++)
    {

        double ua = alpha_f * (*nodes_)[connect_[i]]->getVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(0);
        double va = alpha_f * (*nodes_)[connect_[i]]->getVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(1);

        double uma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(0);
        double vma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(1);

        ua -= uma;
        va -= vma;

        ua_[0] += ua * phi_[i];
        ua_[1] += va * phi_[i];

        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
    };

    double rNorm = sqrt(r[0] * r[0] + r[1] * r[1]);
    double referpar = 0.000000000000000001;

    r[0] /= (rNorm + referpar);
    r[1] /= (rNorm + referpar);

    rr[0][0] = r[0] * r[0];
    rr[0][1] = r[0] * r[1];
    rr[1][0] = r[1] * r[0];
    rr[1][1] = r[1] * r[1];

    hrqd = 2. / ((sqrt(rr[0][0] * MatrixG[0][0] + rr[0][1] * MatrixG[0][1] + rr[1][0] * MatrixG[1][0] + rr[1][1] * MatrixG[1][1])) + referpar);

    // hmin e hmax
    // calculating the eigen-values from G
    double a_ = 1.;
    double b_ = -(MatrixG[0][0] + MatrixG[1][1]);
    double c_ = (MatrixG[0][0] * MatrixG[1][1] - MatrixG[0][1] * MatrixG[1][0]);

    double square = (b_ * b_ - 4. * a_ * c_);

    if (square < 0.000001)
        square = 0.;

    double lambdaMax = (-b_ + sqrt(square)) / (2. * a_);
    double lambdaMin = (-b_ - sqrt(square)) / (2. * a_);

    double hmin = 2. / sqrt(lambdaMax);
    double hmax = 2. / sqrt(lambdaMin);

    if (hrqd < hmin)
        hrqd = hmin;
    if (hrqd > hmax)
        hrqd = hmax;

    // tSUGN1_ (-2)
    tSUGN1_ = (ua_[0] * ua_[0]) * MatrixG[0][0] +
              (ua_[0] * ua_[1]) * MatrixG[0][1] +
              (ua_[1] * ua_[0]) * MatrixG[1][0] +
              (ua_[1] * ua_[1]) * MatrixG[1][1];

    tSUGN2_ = dTime_ / 2.;

    // tSUGN3_ (-1)
    tSUGN3_ = (visc_ / dens_) * (4. / (hrqd * hrqd));

    // Computing tSUPG parameter
    tSUPG_ = 1. / (sqrt(tSUGN1_ + (1. / (tSUGN2_ * tSUGN2_)) + (tSUGN3_ * tSUGN3_)));

    tPSPG_ = tSUPG_;

    tLSIC_ = hrqd * hrqd / tSUPG_;

    for (int i = 0; i < 3; ++i) delete[] MatrixCuInv[i];
    delete[] MatrixCuInv;

    for (int i = 0; i < 3; ++i)delete[] MatrixCvInv[i];
    delete[] MatrixCvInv;

    return;
};

template <>
void Element<2>::getParameterArlequin_FEM(double &tARLQ_, double *phi_, double **dphi_dx)
{

    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();

    double hRGN_ = 0.;
    double hUGN_ = 0.;
    double tSUGN1_, tSUGN2_, tSUGN3_;

    double r[2] = {};
    double s[2] = {};

    double u_ = 0.;
    double v_ = 0.;

    for (int i = 0; i < 6; i++)
    {

        double ua = alpha_f * (*nodes_)[connect_[i]]->getVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(0);
        double va = alpha_f * (*nodes_)[connect_[i]]->getVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(1);

        double uma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(0);
        double vma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(1);

        ua -= uma;
        va -= vma;

        u_ += ua * phi_[i];
        v_ += va * phi_[i];

        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
    };

    double uNorm = sqrt(u_ * u_ + v_ * v_);
    if (uNorm > 1.e-10)
    {
        s[0] = u_ / uNorm;
        s[1] = v_ / uNorm;
    }
    else
    {
        s[0] = 1. / sqrt(2.);
        s[1] = 1. / sqrt(2.);
    };

    double rNorm = sqrt(r[0] * r[0] + r[1] * r[1]);
    if (rNorm >= 1.e-10)
    {
        r[0] /= rNorm;
        r[1] /= rNorm;
    }
    else
    {
        r[0] = 1. / sqrt(2.);
        r[1] = 1. / sqrt(2.);
    };

    for (int i = 0; i < 6; i++)
    {
        hRGN_ += fabs(r[0] * dphi_dx[0][i] + r[1] * dphi_dx[1][i]);
        hUGN_ += fabs(s[0] * dphi_dx[0][i] + s[1] * dphi_dx[1][i]);
    };

    if (hRGN_ >= 1.e-10)
    {
        hRGN_ = 2. / hRGN_;
    }
    else
    {
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10)
    {
        hUGN_ = 2. / hUGN_;
    }
    else
    {
        hUGN_ = 2. / 1.e-10;
    };

    if (uNorm >= 1.e-10)
    {
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }
    else
    {
        tSUGN1_ = hUGN_ / 2.e-10;
    };

    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

    if (fabs(tSUGN1_) <= 1.e-10)
        tSUGN1_ = 1.e-10;
    if (fabs(tSUGN3_) <= 1.e-10)
        tSUGN3_ = 1.e-10;

    // Computing tARLQ parameter
    double tsupg = 1. / sqrt(1. / (tSUGN1_ * tSUGN1_) +
                             1. / (tSUGN2_ * tSUGN2_) +
                             1. / (tSUGN3_ * tSUGN3_));

    tARLQ_ = 1. * k1 * tsupg * 1.e-2;

    return;
}

template <>
void Element<2>::getParameterArlequin2_FEM(int &index, double &djac_, double &weight_, double &tARLQ_, double *phi_,
                                           double *phiC_, double **dphi_dx, double ***ddphi_dx)
{

    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &visc_ = parameters->getViscosity();
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    int na = 6;

    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVel(na, phi_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVel(na, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    //Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMesh_dx[i] = new double[2];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMeshPrev_dx[i] = new double[2];

    double dumeshna_dx[2][2];
    getInterpMeshVelDer(na, dphi_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];
        };
    };

    //Acceleration derivatives
    double **daccel_dx;
    daccel_dx = new double *[2];
    for (int i = 0; i < 2; ++i) daccel_dx[i] = new double[2];

    double **daccelPrev_dx;
    daccelPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) daccelPrev_dx[i] = new double[2];

    double daccelm_dx[2][2];
    getInterpAccelDer(na, dphi_dx, daccel_dx, daccelPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            daccelm_dx[i][j] = alpha_m * daccel_dx[i][j] + (1. - alpha_m) * daccelPrev_dx[i][j];
        };
    };

    //Lagrange Multipliers derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dlambda_dx[i] = new double[2];
    getInterpLambdaDer(na, dphi_dx, dlambda_dx);


    //Second pressure derivatives
    double **ddpress_dxdx;
    ddpress_dxdx = new double *[2];
    for (int i = 0; i < 2; ++i) ddpress_dxdx[i] = new double[2];
    getInterpSecondPressDer(na, ddphi_dx, ddpress_dxdx);


    //Second velocity derivatives
    double ***ddu_dxdx;
    ddu_dxdx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddu_dxdx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
        {
            ddu_dxdx[i][j] = new double[2];
        };
    };

    double ***dduPrev_dxdx;
    dduPrev_dxdx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        dduPrev_dxdx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
        {
            dduPrev_dxdx[i][j] = new double[2];
        };
    };

    double dduna_dxdx[2][2][2];
    getInterpSecondVelDer(na, ddphi_dx, ddu_dxdx, dduPrev_dxdx);
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 2; k++){
                dduna_dxdx[i][j][k] = alpha_f * ddu_dxdx[i][j][k] + (1. - alpha_f) * dduPrev_dxdx[i][j][k];
            };
        };
    };


    double lambda1[12][12] = {};
    double lambda0[12][18] = {};
    double convec[12] = {};
    double iner[12] = {};
    double visc[12] = {};
    double press[12] = {};
    double gradlambda[12] = {};

    double wna_ = 1.;//alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];
    
    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < 6; i++)
    {

        // convection
        convec[2*i] += (((dphi_dx[0][i] * (dduna_dxdx[0][0][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[0][0][1] * (una_[1] - umeshna_[1]))) +
                         (dphi_dx[1][i] * (dduna_dxdx[0][1][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[0][1][1] * (una_[1] - umeshna_[1])))) +
                        ((dphi_dx[0][i] * (duna_dx[0][0] * (duna_dx[0][0] - dumeshna_dx[0][0]) + duna_dx[0][1] * (duna_dx[1][0] - dumeshna_dx[1][0]))) +
                         (dphi_dx[1][i] * (duna_dx[0][0] * (duna_dx[0][1] - dumeshna_dx[0][1]) + duna_dx[0][1] * (duna_dx[1][1] - dumeshna_dx[1][1]))))) * WJ;

        convec[2*i+1] += (((dphi_dx[0][i] * (dduna_dxdx[1][0][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[1][0][1] * (una_[1] - umeshna_[1]))) +
                           (dphi_dx[1][i] * (dduna_dxdx[1][1][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[1][1][1] * (una_[1] - umeshna_[1])))) +
                          ((dphi_dx[0][i] * (duna_dx[1][0] * (duna_dx[0][0] - dumeshna_dx[0][0]) + duna_dx[1][1] * (duna_dx[1][0] - dumeshna_dx[1][0]))) +
                           (dphi_dx[1][i] * (duna_dx[1][0] * (duna_dx[0][1] - dumeshna_dx[0][1] + duna_dx[1][1] * (duna_dx[1][1] - dumeshna_dx[1][1])))))) * WJ;

        // inercia
        iner[2*i] += (dphi_dx[0][i] * daccelm_dx[0][0] + dphi_dx[1][i] * daccelm_dx[0][1]) * WJ;
        iner[2*i+1] += (dphi_dx[0][i] * daccelm_dx[1][0] + dphi_dx[1][i] * daccelm_dx[1][1]) * WJ;

        // viscosity
        visc[2*i] += ((ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * dduna_dxdx[0][0][0] + dduna_dxdx[0][1][1]) +
                      (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * dduna_dxdx[1][0][1]) * WJ * visc_;
        visc[2*i+1] += ((ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * dduna_dxdx[0][1][0] +
                        (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * dduna_dxdx[1][1][1] + dduna_dxdx[1][0][0])) * WJ * visc_;

        // pressure
        press[2*i] += (dphi_dx[0][i] * ddpress_dxdx[0][0] + dphi_dx[1][i] * ddpress_dxdx[1][0]) * WJ;
        press[2*i+1] += (dphi_dx[0][i] * ddpress_dxdx[0][1] + dphi_dx[1][i] * ddpress_dxdx[1][1]) * WJ;

        gradlambda[2*i] += (dphi_dx[0][i] * dlambda_dx[0][0] + dphi_dx[1][i] * dlambda_dx[0][1]) * k1 * WJ / wna_;
        gradlambda[2*i+1] += (dphi_dx[0][i] * dlambda_dx[1][0] + dphi_dx[1][i] * dlambda_dx[1][1]) * k1 * WJ / wna_;

        for (int j = 0; j < 6; j++)
        {
            // lagrange multipliers
            lambda1[2*i][2*j] += phi_[i] * phi_[j] * k1 * WJ/wna_;
            lambda1[2*i+1][2*j+1] += phi_[i] * phi_[j] * k1 * WJ/wna_;
        };
    };

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            // lagrange multipliers
            lambda0[2*i][2*j] += phi_[i] * phiC_[j] * k1 * WJ/wna_;
            lambda0[2*i+1][2*j+1] += phi_[i] * phiC_[j] * k1 * WJ/wna_;
        };
    };

    double normLambda1[12] = {};
    double normLambda0[12] = {};

    double uL[12], uG[18];
    for (int i = 0; i < 6; i++)
    {
        uL[2*i] = (*nodes_)[connect_[i]]->getVelocity(0);
        uL[2*i+1] = (*nodes_)[connect_[i]]->getVelocity(1);
    };
    for (int i = 0; i < 9; i++)
    {
        uG[2*i] = (*nodesC_)[connectC_[i]]->getVelocity(0);
        uG[2*i+1] = (*nodesC_)[connectC_[i]]->getVelocity(1);
    };

    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 12; j++)
        {
            normLambda1[i] += lambda1[i][j] * uL[j];
        };
    };

    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 18; j++)
        {
            normLambda0[i] += lambda0[i][j] * uG[j];
        };
    };

    double vecnormLambda0 = 0.;
    double vecnormLambda1 = 0.;
    double normConvec = 0.;
    double normIner = 0.;
    double normVisc = 0.;
    double normPress = 0.;
    double normGradLambda = 0.;

    for (int i = 0; i < 12; i++)
    {
        vecnormLambda0 += normLambda0[i] * normLambda0[i];
        vecnormLambda1 += normLambda1[i] * normLambda1[i];
        normConvec += convec[i] * convec[i];
        normIner += iner[i] * iner[i];
        normVisc += visc[i] * visc[i];
        normPress += press[i] * press[i];
        normGradLambda += gradlambda[i] * gradlambda[i];
    };

    vecnormLambda0 = sqrt(vecnormLambda0);
    vecnormLambda1 = sqrt(vecnormLambda1);
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

    double tA1 = std::fabs(vecnormLambda1) / std::fabs(normConvec);
    double tA2 = std::fabs(vecnormLambda0) / std::fabs(normConvec);
    double tB1 = std::fabs(vecnormLambda1) / std::fabs(normIner);
    double tB2 = std::fabs(vecnormLambda0) / std::fabs(normIner);
    double tC1 = std::fabs(vecnormLambda1) / std::fabs(normVisc);
    double tC2 = std::fabs(vecnormLambda0) / std::fabs(normVisc);
    double tD1 = std::fabs(vecnormLambda1) / std::fabs(normPress);
    double tD2 = std::fabs(vecnormLambda0) / std::fabs(normPress);
    double tE1 = std::fabs(vecnormLambda1) / std::fabs(normGradLambda);
    double tE2 = std::fabs(vecnormLambda0) / std::fabs(normGradLambda);
    if(tE1>1.)tE1=0.0000000001;
    if(tE2>1.)tE2=0.0000000001;

    double tA = 1. / sqrt(1. / (tA1 * tA1) + 1. / (tA2 * tA2));

    double tB = 1. / sqrt(1. / (tB1 * tB1) + 1. / (tB2 * tB2));

    double tC = 1. / sqrt(1. / (tC1 * tC1) + 1. / (tC2 * tC2));

    double tD = 1. / sqrt(1. / (tD1 * tD1) + 1. / (tD2 * tD2));

    double tE = 1. / sqrt(1. / (tE1 * tE1) + 1. / (tE2 * tE2));

    tARLQ_ = 1. / sqrt(1. / (tA * tA) +
                        1. / (tB * tB) +
                        1. / (tC * tC) +
                        1. / (tE * tE) +
                        1. / (tD * tD));

    //  tARLQ_ = 1./sqrt(1. / (tE1 * tE1)+1./(tE2 * tE2));
    

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < 2; ++i) delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] daccel_dx[i];
    delete[] daccel_dx;

    for (int i = 0; i < 2; ++i) delete[] daccelPrev_dx[i];
    delete[] daccelPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] dlambda_dx[i];
    delete[] dlambda_dx;

    for (int i = 0; i < 2; ++i) delete[] ddpress_dxdx[i];
    delete[] ddpress_dxdx;


    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] ddu_dxdx[i][j];
        };
        delete[] ddu_dxdx[i];
    };
    delete[] ddu_dxdx;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] dduPrev_dxdx[i][j];
        };
        delete[] dduPrev_dxdx[i];
    };
    delete[] dduPrev_dxdx;

};


template <>
void Element<2>::getParameterArlequinMN_FEM(int &index, double &djac_, double &weight_, double &tARLQ_, double *phi_,
                                           double *phiC_, double **dphi_dx, double ***ddphi_dx)
{

    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters-> getGamma();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &dTime_ = parameters->getTimeStep();    
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    int na = 6;

    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVel(na, phi_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVel(na, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    //Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMesh_dx[i] = new double[2];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMeshPrev_dx[i] = new double[2];

    double dumeshna_dx[2][2];
    getInterpMeshVelDer(na, dphi_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];
        };
    };

    double lambda1[12][12] = {};
    double lambda0[12][18] = {};
    double convec[12][12] = {};
    double iner[12][12] = {};
    double visc[12][12] = {};
    double press[12][12] = {};
    double gradlambda[12][12] = {};

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {

           //inercia matrix
            iner[2*i][2*j] += (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * alpha_m;
            iner[2*i+1][2*j+1] += (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * alpha_m;

            //convecction matrix
            double convec1 = (dphi_dx[0][i] * ((una_[0] - umeshna_[0]) * ddphi_dx[0][0][j] + (duna_dx[0][0] - dumeshna_dx[0][0]) * dphi_dx[0][j] + 
                                               (una_[1] - umeshna_[1]) * ddphi_dx[1][0][j] + (duna_dx[1][0] - dumeshna_dx[1][0]) * dphi_dx[1][j]) +  
                              dphi_dx[1][i] * ((una_[0] - umeshna_[0]) * ddphi_dx[0][1][j] + (duna_dx[0][1] - dumeshna_dx[0][1]) * dphi_dx[0][j] + 
                                               (una_[1] - umeshna_[1]) * ddphi_dx[1][1][j] + (duna_dx[1][1] - dumeshna_dx[1][1]) * dphi_dx[1][j]));

            convec[2*i][2*j] += convec1 * alpha_f * gamma * dTime_;
            convec[2*i+1][2*j+1] += convec1 * alpha_f * gamma * dTime_;

            //pressure term
            press[2*i][j] += (dphi_dx[0][i] * ddphi_dx[0][0][j] + dphi_dx[1][i] * ddphi_dx[1][0][j])/dens_;
            press[2*i+1][j] += (dphi_dx[0][i] * ddphi_dx[0][1][j] + dphi_dx[1][i] * ddphi_dx[1][1][j])/dens_;

            gradlambda[2*i][2*j] += (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j])/dens_*k1;
            gradlambda[2*i+1][2*j+1] += (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j])/dens_*k1;
        
            // viscosity
            visc[2*i][2*j]     += (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * ddphi_dx[0][0][j] + ddphi_dx[1][1][j]) * visc_ * alpha_f * gamma * dTime_;
            visc[2*i][2*j+1]   += (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * ddphi_dx[0][1][j] * visc_ * alpha_f * gamma * dTime_;
            visc[2*i+1][2*j]   += (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * ddphi_dx[1][0][j] * visc_ * alpha_f * gamma * dTime_;
            visc[2*i+1][2*j+1] += (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * ddphi_dx[1][1][j] + ddphi_dx[0][0][j]) * visc_ * alpha_f * gamma * dTime_;

            // lagrange multipliers
            lambda1[2*i][2*j] += phi_[i] * phi_[j] * k1;
            lambda1[2*i+1][2*j+1] += phi_[i] * phi_[j] * k1;
        };

        for (int j = 0; j < 9; j++)
        {
            // lagrange multipliers
            lambda0[2*i][2*j] += phi_[i] * phiC_[j] * k1;
            lambda0[2*i+1][2*j+1] += phi_[i] * phiC_[j] * k1;
        };

    };

    double normLambda0 = 0.;
    double normLambda1 = 0.;
    double normConvec = 0.;
    double normIner = 0.;
    double normVisc = 0.;
    double normPress = 0.;
    double normGradLambda = 0.;

    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 12; j++){
            normLambda1 += lambda1[i][j] * lambda1[i][j];
            normConvec += convec[i][j] * convec[i][j];
            normIner += iner[i][j] * iner[i][j];
            normVisc += visc[i][j] * visc[i][j];
            normPress += press[i][j] * press[i][j];
            normGradLambda += gradlambda[i][j] * gradlambda[i][j];
        }
        for (int j = 0; j < 18; j++){
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
    if(tE1>1.)tE1=0.0000000001;
    if(tE2>1.)tE2=0.0000000001;

    double tA = 1. / sqrt(1. / (tA1 * tA1) + 1. / (tA2 * tA2));

    double tB = 1. / sqrt(1. / (tB1 * tB1) + 1. / (tB2 * tB2));

    double tC = 1. / sqrt(1. / (tC1 * tC1) + 1. / (tC2 * tC2));

    double tD = 1. / sqrt(1. / (tD1 * tD1) + 1. / (tD2 * tD2));

    double tE = 1. / sqrt(1. / (tE1 * tE1) + 1. / (tE2 * tE2));

    tARLQ_ = 1. / sqrt(1. / (tA * tA) +
                        1. / (tB * tB) +
                        1. / (tC * tC) +
                        1. / (tE * tE) +
                        1. / (tD * tD));  

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < 2; ++i) delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;

};


template <>
void Element<2>::getParameterArlequin2_COARSE_ISO(double &djac_, double &weight_, double &tARLQ_, double *phi_,
                                                  double *phiC_, double **dphi_dx, double **dphiC_dx, 
                                                  double ***ddphi_dx, double ***ddphiC_dx)
{

    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &visc_ = parameters->getViscosity();
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    int na = 9;

    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVelCoarse(na, phiC_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVelCoarse(na, phiC_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDerCoarse(na, dphiC_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    //Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMesh_dx[i] = new double[2];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMeshPrev_dx[i] = new double[2];

    double dumeshna_dx[2][2];
    getInterpMeshVelDerCoarse(na, dphiC_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];
        };
    };

    //Acceleration derivatives
    double **daccel_dx;
    daccel_dx = new double *[2];
    for (int i = 0; i < 2; ++i) daccel_dx[i] = new double[2];

    double **daccelPrev_dx;
    daccelPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) daccelPrev_dx[i] = new double[2];

    double daccelm_dx[2][2];
    getInterpAccelDerCoarse(na, dphiC_dx, daccel_dx, daccelPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            daccelm_dx[i][j] = alpha_m * daccel_dx[i][j] + (1. - alpha_m) * daccelPrev_dx[i][j];
        };
    };

    //Second pressure derivatives
    double **ddpress_dxdx;
    ddpress_dxdx = new double *[2];
    for (int i = 0; i < 2; ++i) ddpress_dxdx[i] = new double[2];
    getInterpSecondPressDerCoarse(na, ddphiC_dx, ddpress_dxdx);


    //Second velocity derivatives
    double ***ddu_dxdx;
    ddu_dxdx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddu_dxdx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
        {
            ddu_dxdx[i][j] = new double[2];
        };
    };

    double ***dduPrev_dxdx;
    dduPrev_dxdx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        dduPrev_dxdx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
        {
            dduPrev_dxdx[i][j] = new double[2];
        };
    };

    double dduna_dxdx[2][2][2];
    getInterpSecondVelDerCoarse(na, ddphiC_dx, ddu_dxdx, dduPrev_dxdx);
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 2; k++){
                dduna_dxdx[i][j][k] = alpha_f * ddu_dxdx[i][j][k] + (1. - alpha_f) * dduPrev_dxdx[i][j][k];
            };
        };
    };

    na = 6;
    //Lagrange Multipliers derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[2];
    for (int i = 0; i < 2; ++i)
        dlambda_dx[i] = new double[2];
    getInterpLambdaDer(na, dphi_dx, dlambda_dx);


    double lambda1[12][12] = {};
    double lambda0[12][18] = {};
    double convec[12] = {};
    double iner[12] = {};
    double visc[12] = {};
    double press[12] = {};
    double gradlambda[12] = {};

    double WJ = weight_ * djac_;

    for (int i = 0; i < 6; i++)
    {

        // convection
        convec[2*i] += (((dphi_dx[0][i] * (dduna_dxdx[0][0][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[0][0][1] * (una_[1] - umeshna_[1]))) +
                         (dphi_dx[1][i] * (dduna_dxdx[0][1][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[0][1][1] * (una_[1] - umeshna_[1])))) +
                        ((dphi_dx[0][i] * (duna_dx[0][0] * (duna_dx[0][0] - dumeshna_dx[0][0]) + duna_dx[0][1] * (duna_dx[1][0] - dumeshna_dx[1][0]))) +
                         (dphi_dx[1][i] * (duna_dx[0][0] * (duna_dx[0][1] - dumeshna_dx[0][1]) + duna_dx[0][1] * (duna_dx[1][1] - dumeshna_dx[1][1]))))) * WJ;

        convec[2*i+1] += (((dphi_dx[0][i] * (dduna_dxdx[1][0][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[1][0][1] * (una_[1] - umeshna_[1]))) +
                           (dphi_dx[1][i] * (dduna_dxdx[1][1][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[1][1][1] * (una_[1] - umeshna_[1])))) +
                          ((dphi_dx[0][i] * (duna_dx[1][0] * (duna_dx[0][0] - dumeshna_dx[0][0]) + duna_dx[1][1] * (duna_dx[1][0] - dumeshna_dx[1][0]))) +
                           (dphi_dx[1][i] * (duna_dx[1][0] * (duna_dx[0][1] - dumeshna_dx[0][1] + duna_dx[1][1] * (duna_dx[1][1] - dumeshna_dx[1][1])))))) * WJ; 

        // inercia
        iner[2*i] += (dphi_dx[0][i] * daccelm_dx[0][0] + dphi_dx[1][i] * daccelm_dx[0][1]) * WJ;
        iner[2*i+1] += (dphi_dx[0][i] * daccelm_dx[1][0] + dphi_dx[1][i] * daccelm_dx[1][1]) * WJ;

        // viscosity
        visc[2*i] += ((ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * dduna_dxdx[0][0][0] + dduna_dxdx[0][1][1]) +
                      (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * dduna_dxdx[1][0][1]) * WJ * visc_;
        visc[2*i+1] += ((ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * dduna_dxdx[0][1][0] +
                        (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * dduna_dxdx[1][1][1] + dduna_dxdx[1][0][0])) * WJ * visc_;

        // pressure
        press[2*i] += (dphi_dx[0][i] * ddpress_dxdx[0][0] + dphi_dx[1][i] * ddpress_dxdx[1][0]) * WJ;
        press[2*i+1] += (dphi_dx[0][i] * ddpress_dxdx[0][1] + dphi_dx[1][i] * ddpress_dxdx[1][1]) * WJ;

        gradlambda[2*i] += (dphi_dx[0][i] * dlambda_dx[0][0] + dphi_dx[1][i] * dlambda_dx[0][1]) * WJ * k1;
        gradlambda[2*i+1] += (dphi_dx[0][i] * dlambda_dx[1][0] + dphi_dx[1][i] * dlambda_dx[1][1]) * WJ * k1;

        for (int j = 0; j < 6; j++)
        {
            // lagrange multipliers
            lambda1[2*i][2*j] += phi_[i] * phi_[j] * WJ * k1;
            lambda1[2*i+1][2*j+1] += phi_[i] * phi_[j] * WJ * k1;
        };
    };

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            // lagrange multipliers
            lambda0[2*i][2*j] += phi_[i] * phiC_[j] * WJ * k1;
            lambda0[2*i+1][2*j+1] += phi_[i] * phiC_[j] * WJ * k1;
        };
    };

    double normLambda1[12] = {};
    double normLambda0[12] = {};

    double uL[12], uG[18];
    for (int i = 0; i < 6; i++)
    {
        uL[2*i] = (*nodes_)[connect_[i]]->getVelocity(0);
        uL[2*i+1] = (*nodes_)[connect_[i]]->getVelocity(1);
    };
    for (int i = 0; i < 9; i++)
    {
        uG[2*i] = (*nodesC_)[connectC_[i]]->getVelocity(0);
        uG[2*i+1] = (*nodesC_)[connectC_[i]]->getVelocity(1);
    };

    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 12; j++)
        {
            normLambda1[i] += lambda1[i][j] * uL[j];
        };
    };

    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 18; j++)
        {
            normLambda0[i] += lambda0[i][j] * uG[j];
        };
    };

    double vecnormLambda0 = 0.;
    double vecnormLambda1 = 0.;
    double normConvec = 0.;
    double normIner = 0.;
    double normVisc = 0.;
    double normPress = 0.;
    double normGradLambda = 0.;

    for (int i = 0; i < 12; i++)
    {
        vecnormLambda0 += normLambda0[i] * normLambda0[i];
        vecnormLambda1 += normLambda1[i] * normLambda1[i];
        normConvec += convec[i] * convec[i];
        normIner += iner[i] * iner[i];
        normVisc += visc[i] * visc[i];
        normPress += press[i] * press[i];
        normGradLambda += gradlambda[i] * gradlambda[i];
    };

    vecnormLambda0 = sqrt(vecnormLambda0);
    vecnormLambda1 = sqrt(vecnormLambda1);
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

    double tA1 = std::fabs(vecnormLambda1) / std::fabs(normConvec);
    double tA2 = std::fabs(vecnormLambda0) / std::fabs(normConvec);
    double tB1 = std::fabs(vecnormLambda1) / std::fabs(normIner);
    double tB2 = std::fabs(vecnormLambda0) / std::fabs(normIner);
    double tC1 = std::fabs(vecnormLambda1) / std::fabs(normVisc);
    double tC2 = std::fabs(vecnormLambda0) / std::fabs(normVisc);
    double tD1 = std::fabs(vecnormLambda1) / std::fabs(normPress);
    double tD2 = std::fabs(vecnormLambda0) / std::fabs(normPress);
    double tE1 = std::fabs(vecnormLambda1) / std::fabs(normGradLambda);
    double tE2 = std::fabs(vecnormLambda0) / std::fabs(normGradLambda);
    if(tE1>1.)tE1=0.0000000001;
    if(tE2>1.)tE2=0.0000000001;

    double tA = 1. / sqrt(1. / (tA1 * tA1) +
                          1. / (tA2 * tA2));

    double tB = 1. / sqrt(1. / (tB1 * tB1) +
                          1. / (tB2 * tB2));

    double tC = 1. / sqrt(1. / (tC1 * tC1) +
                          1. / (tC2 * tC2));

    double tD = 1. / sqrt(1. / (tD1 * tD1) +
                          1. / (tD2 * tD2));

    double tE = 1. / sqrt(1. / (tE1 * tE1) +
                          1. / (tE2 * tE2));

    tARLQ_ = 1. / sqrt(1. / (tA * tA) +
                        1. / (tB * tB) +
                        1. / (tC * tC) +
                        1. / (tE * tE) +
                        1. / (tD * tD));


    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] duMesh_dx[i];
    delete[] duMesh_dx;
    
    for (int i = 0; i < 2; ++i) delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] daccel_dx[i];
    delete[] daccel_dx;
    
    for (int i = 0; i < 2; ++i) delete[] daccelPrev_dx[i];
    delete[] daccelPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] dlambda_dx[i];
    delete[] dlambda_dx;

    for (int i = 0; i < 2; ++i) delete[] ddpress_dxdx[i];
    delete[] ddpress_dxdx;


    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] ddu_dxdx[i][j];
        };
        delete[] ddu_dxdx[i];
    };
    delete[] ddu_dxdx;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] dduPrev_dxdx[i][j];
        };
        delete[] dduPrev_dxdx[i];
    };
    delete[] dduPrev_dxdx;

};

template <>
void Element<2>::getParameterArlequin_ISO(double &tARLQ_, double *phi_, double **dphi_dx)
{

    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();

    double hRGN_ = 0.;
    double hUGN_ = 0.;
    double tSUGN1_, tSUGN2_, tSUGN3_;

    double r[2] = {};
    double s[2] = {};

    double u_ = 0.;
    double v_ = 0.;

    for (int i = 0; i < 9; i++)
    {

        double ua = alpha_f * (*nodes_)[connect_[i]]->getVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(0);
        double va = alpha_f * (*nodes_)[connect_[i]]->getVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousVelocity(1);

        double uma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(0);
        double vma = alpha_f * (*nodes_)[connect_[i]]->getMeshVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]]->getPreviousMeshVelocity(1);

        ua -= uma;
        va -= vma;

        u_ += ua * phi_[i];
        v_ += va * phi_[i];

        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
    };

    double uNorm = sqrt(u_ * u_ + v_ * v_);
    if (uNorm > 1.e-10)
    {
        s[0] = u_ / uNorm;
        s[1] = v_ / uNorm;
    }
    else
    {
        s[0] = 1. / sqrt(2.);
        s[1] = 1. / sqrt(2.);
    };

    double rNorm = sqrt(r[0] * r[0] + r[1] * r[1]);

    if (rNorm >= 1.e-10)
    {
        r[0] /= rNorm;
        r[1] /= rNorm;
    }
    else
    {
        r[0] = 1. / sqrt(2.);
        r[1] = 1. / sqrt(2.);
    };

    for (int i = 0; i < 9; i++)
    {
        hRGN_ += fabs(r[0] * dphi_dx[0][i] + r[1] * dphi_dx[1][i]);
        hUGN_ += fabs(s[0] * dphi_dx[0][i] + s[1] * dphi_dx[1][i]);
    };

    if (hRGN_ >= 1.e-10)
    {
        hRGN_ = 2. / hRGN_;
    }
    else
    {
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10)
    {
        hUGN_ = 2. / hUGN_;
    }
    else
    {
        hUGN_ = 2. / 1.e-10;
    };

    if (uNorm >= 1.e-10)
    {
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }
    else
    {
        tSUGN1_ = hUGN_ / 2.e-10;
    };

    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

    if (fabs(tSUGN1_) <= 1.e-10)
        tSUGN1_ = 1.e-10;
    if (fabs(tSUGN3_) <= 1.e-10)
        tSUGN3_ = 1.e-10;

    // Computing tARLQ parameter
    double tsupg = 1. / sqrt(1. / (tSUGN1_ * tSUGN1_) +
                             1. / (tSUGN2_ * tSUGN2_) +
                             1. / (tSUGN3_ * tSUGN3_));

    tARLQ_ = -1. * k1 * tsupg * 1.e-2;

    return;
}


template <>
void Element<2>::getParameterArlequin_COARSE_ISO(double &tARLQ_, double *phi_, double **dphi_dx)
{

    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();

    double hRGN_ = 0.;
    double hUGN_ = 0.;
    double tSUGN1_, tSUGN2_, tSUGN3_;

    double r[2] = {};
    double s[2] = {};

    double u_ = 0.;
    double v_ = 0.;

    for (int i = 0; i < 9; i++)
    {

        double ua = alpha_f * (*nodesC_)[connectC_[i]]->getVelocity(0) + (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousVelocity(0);
        double va = alpha_f * (*nodesC_)[connectC_[i]]->getVelocity(1) + (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousVelocity(1);

        double uma = alpha_f * (*nodesC_)[connectC_[i]]->getMeshVelocity(0) + (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousMeshVelocity(0);
        double vma = alpha_f * (*nodesC_)[connectC_[i]]->getMeshVelocity(1) + (1. - alpha_f) * (*nodesC_)[connectC_[i]]->getPreviousMeshVelocity(1);

        ua -= uma;
        va -= vma;

        u_ += ua * phi_[i];
        v_ += va * phi_[i];

        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
    };

    double uNorm = sqrt(u_ * u_ + v_ * v_);
    if (uNorm > 1.e-10)
    {
        s[0] = u_ / uNorm;
        s[1] = v_ / uNorm;
    }
    else
    {
        s[0] = 1. / sqrt(2.);
        s[1] = 1. / sqrt(2.);
    };

    double rNorm = sqrt(r[0] * r[0] + r[1] * r[1]);

    if (rNorm >= 1.e-10)
    {
        r[0] /= rNorm;
        r[1] /= rNorm;
    }
    else
    {
        r[0] = 1. / sqrt(2.);
        r[1] = 1. / sqrt(2.);
    };

    for (int i = 0; i < 9; i++)
    {
        hRGN_ += fabs(r[0] * dphi_dx[0][i] + r[1] * dphi_dx[1][i]);
        hUGN_ += fabs(s[0] * dphi_dx[0][i] + s[1] * dphi_dx[1][i]);
    };

    if (hRGN_ >= 1.e-10)
    {
        hRGN_ = 2. / hRGN_;
    }
    else
    {
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10)
    {
        hUGN_ = 2. / hUGN_;
    }
    else
    {
        hUGN_ = 2. / 1.e-10;
    };

    if (uNorm >= 1.e-10)
    {
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }
    else
    {
        tSUGN1_ = hUGN_ / 2.e-10;
    };

    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

    if (fabs(tSUGN1_) <= 1.e-10)
        tSUGN1_ = 1.e-10;
    if (fabs(tSUGN3_) <= 1.e-10)
        tSUGN3_ = 1.e-10;

    // Computing tARLQ parameter
    double tsupg = 1. / sqrt(1. / (tSUGN1_ * tSUGN1_) +
                             1. / (tSUGN2_ * tSUGN2_) +
                             1. / (tSUGN3_ * tSUGN3_));

    tARLQ_ = -1. * k1 * tsupg * 1.e-2;

    return;
}

//------------------------------------------------------------------------------
//-----------------------------ELEMENT LOCAL MATRIX-----------------------------
//------------------------------------------------------------------------------

template <>
void Element<2>::getElemMatrix(int &na, double &djac_, double &weight_, double &tSUPG_, double &tPSPG_, double &tLSIC_, int &index,
                               double *phi_, double **dphi_dx, double **jacobianNRMatrix)
{

    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();

    // Interpolated variables
    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVel(na, phi_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVel(na, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    double wna_;
    if (na == 6)wna_ = alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];
    if (na == 9)wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];
    
    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < na; j++)
        {

            double wSUPGi = (una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i];
            double wSUPGj = (una_[0] - umeshna_[0]) * dphi_dx[0][j] + (una_[1] - umeshna_[1]) * dphi_dx[1][j];

            // Mass matrix (use for both directions)
            double mM = phi_[i] * phi_[j] * dens_ * alpha_m +
                        wSUPGi * phi_[j] * tSUPG_ * dens_ * alpha_m;

            // Difusion matrix (viscosity)
            double Kxx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * visc_ * alpha_f * gamma * dTime_;
            double Kxy = dphi_dx[1][i] * dphi_dx[0][j] * visc_ * alpha_f * gamma * dTime_;
            double Kyx = dphi_dx[0][i] * dphi_dx[1][j] * visc_ * alpha_f * gamma * dTime_;
            double Kyy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + dphi_dx[0][i] * dphi_dx[0][j]) * visc_ * alpha_f * gamma * dTime_;

            // Convection matrixes
            double Cxx = wSUPGj * phi_[i] * dens_ * alpha_f * gamma * dTime_ +
                         wSUPGi * wSUPGj * tSUPG_ * dens_ * alpha_f * gamma * dTime_;
            double Cyy = Cxx;

            double Cuu = phi_[i] * duna_dx[0][0] * phi_[j] * dens_ * alpha_f * gamma * dTime_;
            double Cuv = phi_[i] * duna_dx[0][1] * phi_[j] * dens_ * alpha_f * gamma * dTime_;
            double Cvu = phi_[i] * duna_dx[1][0] * phi_[j] * dens_ * alpha_f * gamma * dTime_;
            double Cvv = phi_[i] * duna_dx[1][1] * phi_[j] * dens_ * alpha_f * gamma * dTime_;

            // Stabilization LSIC matrix
            double KLSxx = dphi_dx[0][i] * dphi_dx[0][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSxy = dphi_dx[0][i] * dphi_dx[1][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSyx = dphi_dx[1][i] * dphi_dx[0][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSyy = dphi_dx[1][i] * dphi_dx[1][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;

            jacobianNRMatrix[2*i][2*j] += (mM + Kxx + KLSxx + Cxx + Cuu) * WJ;
            jacobianNRMatrix[2*i+1][2*j+1] += (mM + Kyy + KLSyy + Cyy + Cvv) * WJ;
            jacobianNRMatrix[2*i][2*j+1] += (Kxy + Cuv + KLSxy) * WJ;
            jacobianNRMatrix[2*i+1][2*j] += (Kyx + Cvu + KLSyx) * WJ;

            // multipy pressure direction x and y
            double QSUPGx = -(dphi_dx[0][i] * phi_[j]) +
                            wSUPGi * dphi_dx[0][j] * tSUPG_;
            double QSUPGy = -(dphi_dx[1][i] * phi_[j]) +
                            wSUPGi * dphi_dx[1][j] * tSUPG_;

            // multiply velocity direction x and y
            double Qx = dphi_dx[0][j] * phi_[i] * alpha_f * gamma * dTime_;
            double Qy = dphi_dx[1][j] * phi_[i] * alpha_f * gamma * dTime_;

            jacobianNRMatrix[2*na+i][2*j] += Qx * WJ;
            jacobianNRMatrix[2*na+i][2*j+1] += Qy * WJ;
            jacobianNRMatrix[2*i][2*na+j] += QSUPGx * WJ;
            jacobianNRMatrix[2*i+1][2*na+j] += QSUPGy * WJ;

            // PSPG stabilization matrixes
            double Hx = dphi_dx[0][i] * phi_[j] * tPSPG_ * alpha_m;
            double Hy = dphi_dx[1][i] * phi_[j] * tPSPG_ * alpha_m;

            double Gx = dphi_dx[0][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;
            double Gy = dphi_dx[1][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;

            double Q = (dphi_dx[0][i] * dphi_dx[0][j] +
                        dphi_dx[1][i] * dphi_dx[1][j]) *
                       tPSPG_ / (dens_);

            jacobianNRMatrix[2*na+i][2*j] += (Hx + Gx) * WJ;
            jacobianNRMatrix[2*na+i][2*j+1] += (Hy + Gy) * WJ;
            jacobianNRMatrix[2*na+i][2*na+j] += Q * WJ;
        };
    };

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;
    
    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    return;
};


template <>
void Element<2>::getMatrixAndVectorsSameMesh(int &na, double &djac_, double &weight_, double *phi_, double **dphi_dx,
                                             double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    // Fluid Data
    double &alpha_f = parameters->getAlphaF();
    double &dens_ = parameters->getDensity();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    // Interpolated variables

    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVel(na, phi_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    //Lagrange Multiplieres
    double lambda_[2];
    getInterpLambda(na, phi_, lambda_);

    //Lagrange Multiplieres derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dlambda_dx[i] = new double[2];
    getInterpLambdaDer(na, dphi_dx, dlambda_dx);;


    double WJ = weight_ * djac_;

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < na; j++)
        {
            // L2 operator
            double L2 = -phi_[i] * phi_[j] * k1;
            lagrMultMatrix[2*i][2*j] += L2 * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += L2 * WJ;

            // H1 operator
            double H1xx = -(dphi_dx[0][i] * dphi_dx[0][j] +
                           0.5 * dphi_dx[1][i] * dphi_dx[1][j]) * k2;
            double H1xy = -0.5 * dphi_dx[1][i] * dphi_dx[0][j] * k2;
            double H1yx = -0.5 * dphi_dx[0][i] * dphi_dx[1][j] * k2;
            double H1yy = -(dphi_dx[1][i] * dphi_dx[1][j] +
                           0.5 * dphi_dx[0][i] * dphi_dx[0][j]) * k2;

            lagrMultMatrix[2*i][2*j] += H1xx * WJ;
            lagrMultMatrix[2*i+1][2*j] += H1yx * WJ;
            lagrMultMatrix[2*i][2*j+1] += H1xy * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += H1yy * WJ;
        };

        // L2 operator - Lagrange
        double l2x_ = phi_[i] * lambda_[0] * k1;
        double l2y_ = phi_[i] * lambda_[1] * k1;

        // H1 operator - Lagrange
        double h1x_ = (dphi_dx[0][i] * dlambda_dx[0][0] +
                       0.5 * dphi_dx[1][i] * dlambda_dx[0][1] +
                       0.5 * dphi_dx[1][i] * dlambda_dx[1][0]) * k2;

        double h1y_ = (0.5 * dphi_dx[0][i] * dlambda_dx[0][1] +
                       dphi_dx[1][i] * dlambda_dx[1][1] +
                       0.5 * dphi_dx[0][i] * dlambda_dx[1][0]) * k2;

        rhsVector1[2*i] += (l2x_ + h1x_) * WJ;
        rhsVector1[2*i+1] += (l2y_ + h1y_) * WJ;

        // L2 operator - Velocity
        double l2ux_ = phi_[i] * una_[0] * k1;
        double l2uy_ = phi_[i] * una_[1] * k1;

        // H1 operator - Velocity
        double h1ux_ = (dphi_dx[0][i] * duna_dx[0][0] +
                        0.5 * dphi_dx[1][i] * duna_dx[0][1] +
                        0.5 * dphi_dx[1][i] * duna_dx[1][0]) * k2;

        double h1uy_ = (0.5 * dphi_dx[0][i] * duna_dx[0][1] +
                        dphi_dx[1][i] * duna_dx[1][1] +
                        0.5 * dphi_dx[0][i] * duna_dx[1][0]) * k2;

        rhsVector2[2*i] += (l2ux_ + h1ux_) * WJ;
        rhsVector2[2*i+1] += (l2uy_ + h1uy_) * WJ;
    };

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;
    
    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] dlambda_dx[i];
    delete[] dlambda_dx;
};

template <>
void Element<2>::getMatrixAndVectorsSameMesh_tSUPG_tPSPG(int &na, double &djac_, double &weight_, double &tSUPG_, double &tPSPG_,
                                                         double *phi_, double **dphi_dx, double **jacobianNRMatrix, double *rhsVector)
{

    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVel(na, phi_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVel(na, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Lagrange Multiplieres
    double lambda_[2];
    getInterpLambda(na, phi_, lambda_);

    double WJ = weight_ * djac_;

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < na; j++)
        {

            double msupg = -((una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i]) * phi_[j] * tSUPG_ * k1;
            jacobianNRMatrix[2*i][2*j] += msupg * WJ;
            jacobianNRMatrix[2*i+1][2*j+1] += msupg * WJ;

            double mpspgx = -(dphi_dx[0][i] * phi_[j]) * tPSPG_ * k1 / dens_;
            double mpspgy = -(dphi_dx[1][i] * phi_[j]) * tPSPG_ * k1 / dens_;

            jacobianNRMatrix[2*na+i][2*j] += mpspgx * WJ;
            jacobianNRMatrix[2*na+i][2*j+1] += mpspgy * WJ;
        };

        double vsupgx = ((una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i]) * lambda_[0] * tSUPG_ * k1;
        double vsupgy = ((una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i]) * lambda_[1] * tSUPG_ * k1;

        double vpspg = (dphi_dx[0][i] * lambda_[0] + dphi_dx[1][i] * lambda_[1]) * tPSPG_ * k1/ dens_;

        rhsVector[2*i] += vsupgx * WJ;
        rhsVector[2*i+1] += vsupgy * WJ;
        rhsVector[2*na+i] += vpspg * WJ;
    };
};

template <>
void Element<2>::getMatrixAndVectorsSameMeshArlqStab(int &na, double &djac_, double &weight_, double &tARLQ_, int &index,
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
    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVel(na, phi_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVel(na, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    //Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMesh_dx[i] = new double[2];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMeshPrev_dx[i] = new double[2];

    double dumeshna_dx[2][2];
    getInterpMeshVelDer(na, dphi_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];
        };
    };

    //Acceleration derivatives
    double **daccel_dx;
    daccel_dx = new double *[2];
    for (int i = 0; i < 2; ++i) daccel_dx[i] = new double[2];

    double **daccelPrev_dx;
    daccelPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) daccelPrev_dx[i] = new double[2];

    double daccelm_dx[2][2];
    getInterpAccelDer(na, dphi_dx, daccel_dx, daccelPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            daccelm_dx[i][j] = alpha_m * daccel_dx[i][j] + (1. - alpha_m) * daccelPrev_dx[i][j];
        };
    };

    //Lagrange Multiplieres derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dlambda_dx[i] = new double[2];
    getInterpLambdaDer(na, dphi_dx, dlambda_dx);

    //Second pressure derivatives
    double **ddpress_dxdx;
    ddpress_dxdx = new double *[2];
    for (int i = 0; i < 2; ++i) ddpress_dxdx[i] = new double[2];
    getInterpSecondPressDer(na, ddphi_dx, ddpress_dxdx);

    //Second velocity derivatives
    double ***ddu_dxdx;
    ddu_dxdx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddu_dxdx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
        {
            ddu_dxdx[i][j] = new double[2];
        };
    };

    double ***dduPrev_dxdx;
    dduPrev_dxdx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        dduPrev_dxdx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
        {
            dduPrev_dxdx[i][j] = new double[2];
        };
    };

    double dduna_dxdx[2][2][2];
    getInterpSecondVelDer(na, ddphi_dx, ddu_dxdx, dduPrev_dxdx);    
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 2; k++){
                dduna_dxdx[i][j][k] = alpha_f * ddu_dxdx[i][j][k] + (1. - alpha_f) * dduPrev_dxdx[i][j][k];
            };
        };
    };

    double wna_;
    if (na == 6) wna_ = alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];
    if (na == 9) wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];
    
    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < na; j++)
        {

            // tarlq x mass matrix
            double mass = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_;
            arlequinStab1[2*i][2*j] += mass * WJ * alpha_m;
            arlequinStab1[2*i+1][2*j+1] += mass * WJ * alpha_m;

 
            // tarlq x convecction matrix
            double convec1 = (dphi_dx[0][i] * ((una_[0] - umeshna_[0]) * ddphi_dx[0][0][j] + (duna_dx[0][0] - dumeshna_dx[0][0]) * dphi_dx[0][j] + 
                                               (una_[1] - umeshna_[1]) * ddphi_dx[1][0][j] + (duna_dx[1][0] - dumeshna_dx[1][0]) * dphi_dx[1][j]) +  
                              dphi_dx[1][i] * ((una_[0] - umeshna_[0]) * ddphi_dx[0][1][j] + (duna_dx[0][1] - dumeshna_dx[0][1]) * dphi_dx[0][j] + 
                                               (una_[1] - umeshna_[1]) * ddphi_dx[1][1][j] + (duna_dx[1][1] - dumeshna_dx[1][1]) * dphi_dx[1][j])) * tARLQ_;

            arlequinStab1[2*i][2*j] += convec1 * WJ * alpha_f * gamma * dTime_;
            arlequinStab1[2*i+1][2*j+1] += convec1 * WJ * alpha_f * gamma * dTime_;

            //tarlq x pressure term
            double pressx = (dphi_dx[0][i] * ddphi_dx[0][0][j] + dphi_dx[1][i] * ddphi_dx[1][0][j]) * tARLQ_/dens_;
            double pressy = (dphi_dx[0][i] * ddphi_dx[0][1][j] + dphi_dx[1][i] * ddphi_dx[1][1][j]) * tARLQ_/dens_;
            arlequinStab1[2*i][2*na+j] += pressx * WJ;
            arlequinStab1[2*i+1][2*na+j] += pressy * WJ;

            // Diagonal matrix
            double LL = -(dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ * k1 / dens_;

            arlequinStabD[2*i][2*j] += LL * WJ / wna_;
            arlequinStabD[2*i+1][2*j+1] += LL * WJ / wna_;


           
        };

        // tarlq x mass matrix
        double massx = -(dphi_dx[0][i] * daccelm_dx[0][0] + dphi_dx[1][i] * daccelm_dx[0][1]) * tARLQ_;
        double massy = -(dphi_dx[0][i] * daccelm_dx[1][0] + dphi_dx[1][i] * daccelm_dx[1][1]) * tARLQ_;

        arlequinStabVector1[2*i] += massx * WJ;
        arlequinStabVector1[2*i+1] += massy * WJ;

        // tarlq x convecction matrix
        double convec1x = -((dphi_dx[0][i] * (dduna_dxdx[0][0][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[0][1][0] * (una_[1] - umeshna_[1]))) + 
                           (dphi_dx[1][i] * (dduna_dxdx[0][0][1] * (una_[0] - umeshna_[0]) + dduna_dxdx[0][1][1] * (una_[1] - umeshna_[1])))) * tARLQ_;

        double convec2x = -((dphi_dx[0][i] * (duna_dx[0][0] * (duna_dx[0][0] - dumeshna_dx[0][0]) + duna_dx[0][1] * (duna_dx[1][0] - dumeshna_dx[1][0]))) +   
                           (dphi_dx[1][i] * (duna_dx[0][0] * (duna_dx[0][1] - dumeshna_dx[0][1]) + duna_dx[0][1] * (duna_dx[1][1] - dumeshna_dx[1][1])))) * tARLQ_;

        double convec1y = -((dphi_dx[0][i] * (dduna_dxdx[1][0][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[1][1][0] * (una_[1] - umeshna_[1]))) + 
                           (dphi_dx[1][i] * (dduna_dxdx[1][0][1] * (una_[0] - umeshna_[0]) + dduna_dxdx[1][1][1] * (una_[1] - umeshna_[1])))) * tARLQ_;

        double convec2y = -((dphi_dx[0][i] * (duna_dx[1][0] * (duna_dx[0][0] - dumeshna_dx[0][0]) + duna_dx[1][1] * (duna_dx[1][0] - dumeshna_dx[1][0]))) +   
                           (dphi_dx[1][i] * (duna_dx[1][0] * (duna_dx[0][1] - dumeshna_dx[0][1]) + duna_dx[1][1] * (duna_dx[1][1] - dumeshna_dx[1][1])))) * tARLQ_;

        arlequinStabVector1[2*i] += (convec1x + convec2x) * WJ;
        arlequinStabVector1[2*i+1] += (convec1y + convec2y) * WJ;

        // tarlq x pressure term
         double pressx = -(dphi_dx[0][i] * ddpress_dxdx[0][0] + dphi_dx[1][i] * ddpress_dxdx[1][0]) * tARLQ_/dens_;
         double pressy = -(dphi_dx[0][i] * ddpress_dxdx[0][1] + dphi_dx[1][i] * ddpress_dxdx[1][1]) * tARLQ_/dens_;

        arlequinStabVector1[2*i] += pressx * WJ;
        arlequinStabVector1[2*i+1] += pressy * WJ;

        // Diagonal matrix
        double LLx = (dphi_dx[0][i] * dlambda_dx[0][0] + dphi_dx[1][i] * dlambda_dx[0][1]) * tARLQ_ * k1 / dens_;
        double LLy = (dphi_dx[0][i] * dlambda_dx[1][0] + dphi_dx[1][i] * dlambda_dx[1][1]) * tARLQ_ * k1 / dens_;

        arlequinStabVectorD[2*i] += LLx * WJ/ wna_;
        arlequinStabVectorD[2*i+1] += LLy * WJ/ wna_;
    };

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < 2; ++i) delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] daccel_dx[i];
    delete[] daccel_dx;

    for (int i = 0; i < 2; ++i) delete[] daccelPrev_dx[i];
    delete[] daccelPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] dlambda_dx[i];
    delete[] dlambda_dx;

    for (int i = 0; i < 2; ++i) delete[] ddpress_dxdx[i];
    delete[] ddpress_dxdx;


    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] ddu_dxdx[i][j];
        };
        delete[] ddu_dxdx[i];
    };
    delete[] ddu_dxdx;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] dduPrev_dxdx[i][j];
        };
        delete[] dduPrev_dxdx[i];
    };
    delete[] dduPrev_dxdx;

};

template <>
void Element<2>::getMatrixAndVectorsDifferentMesh(double &djac_, double &weight_, int &na, double *phi_,
                                                  double **dphi_dx, int &naC, double *phiC_, double **dphiC_dx,
                                                  double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)


{
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    // Interpolated variables
    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVelCoarse(naC, phiC_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDerCoarse(naC, dphiC_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    double lambda_[2];
    getInterpLambda(na, phi_, lambda_);

    //Lagrange Multipliers derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dlambda_dx[i] = new double[2];
    getInterpLambdaDer(na, dphi_dx, dlambda_dx);

    double WJ = weight_ * djac_;

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < naC; j++)
        {

            // L2 operator
            double L2 = phi_[i] * phiC_[j] * k1;

            lagrMultMatrix[2*i][2*j] += L2 * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += L2 * WJ;

            // H1 operator
            double H1xx = (dphi_dx[0][i] * dphiC_dx[0][j] +
                           0.5 * dphi_dx[1][i] * dphiC_dx[1][j]) * k2;
            double H1xy = 0.5 * dphi_dx[1][i] * dphiC_dx[0][j] * k2;
            double H1yx = 0.5 * dphi_dx[0][i] * dphiC_dx[1][j] * k2;
            double H1yy = (dphi_dx[1][i] * dphiC_dx[1][j] +
                           0.5 * dphi_dx[0][i] * dphiC_dx[0][j]) * k2;

            lagrMultMatrix[2*i][2*j] += H1xx * WJ;
            lagrMultMatrix[2*i+1][2*j] += H1yx * WJ;
            lagrMultMatrix[2*i][2*j+1] += H1xy * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += H1yy * WJ;
        };

        // L2 operator - Velocity
        double l2ux_ = -phi_[i] * una_[0] * k1;
        double l2uy_ = -phi_[i] * una_[1] * k1;

        // H1 operator - Velocity
        double h1ux_ = -(dphi_dx[0][i] * duna_dx[0][0] +
                        0.5 * dphi_dx[1][i] * duna_dx[0][1] +
                        0.5 * dphi_dx[1][i] * duna_dx[1][0]) * k2;

        double h1uy_ = -(0.5 * dphi_dx[0][i] * duna_dx[0][1] +
                        dphi_dx[1][i] * duna_dx[1][1] +
                        0.5 * dphi_dx[0][i] * duna_dx[1][0]) * k2;

        rhsVector2[2*i] += (l2ux_ + h1ux_) * WJ;
        rhsVector2[2*i+1] += (l2uy_ + h1uy_) * WJ;
    };

    for (int i = 0; i < naC; i++)
    {

        // L2 operator - Lagrange
        double l2x_ = -phiC_[i] * lambda_[0] * k1;
        double l2y_ = -phiC_[i] * lambda_[1] * k1;

        // H1 operator - Lagrange
        double h1x_ = -(dphiC_dx[0][i] * dlambda_dx[0][0] +
                       0.5 * dphiC_dx[1][i] * dlambda_dx[0][1] +
                       0.5 * dphiC_dx[1][i] * dlambda_dx[1][0]) * k2;

        double h1y_ = -(0.5 * dphiC_dx[0][i] * dlambda_dx[0][1] +
                       dphiC_dx[1][i] * dlambda_dx[1][1] +
                       0.5 * dphiC_dx[0][i] * dlambda_dx[1][0]) * k2;

        rhsVector1[2*i] += (l2x_ + h1x_) * WJ;
        rhsVector1[2*i+1] += (l2y_ + h1y_) * WJ;
    };

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] dlambda_dx[i];
    delete[] dlambda_dx;
   
};

template <>
void Element<2>::getMatrixAndVectorsDifferentMeshTemp(double &djac_, double &weight_, double &tARLQ_, int &na, double *phi_,
                                                      double **dphi_dx, int &naC, double *phiC_, double **dphiC_dx,
                                                      double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                      double **arlequinStab, double *arlequinStabVector)

{
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    // Interpolated variables
    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVelCoarse(naC, phiC_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDerCoarse(naC, dphiC_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    double lambda_[2];
    getInterpLambda(na, phi_, lambda_);

    //Lagrange Multipliers derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dlambda_dx[i] = new double[2];
    getInterpLambdaDer(na, dphi_dx, dlambda_dx);

    double WJ = weight_ * djac_;

    for (int i = 0; i < na; i++){
        for (int j = 0; j < naC; j++){

            // L2 operator
            double L2 = phi_[i] * phiC_[j] * k1;

            lagrMultMatrix[2*i][2*j] += L2 * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += L2 * WJ;

            // H1 operator
            double H1xx = (dphi_dx[0][i] * dphiC_dx[0][j] +
                           0.5 * dphi_dx[1][i] * dphiC_dx[1][j]) * k2;
            double H1xy = 0.5 * dphi_dx[1][i] * dphiC_dx[0][j] * k2;
            double H1yx = 0.5 * dphi_dx[0][i] * dphiC_dx[1][j] * k2;
            double H1yy = (dphi_dx[1][i] * dphiC_dx[1][j] +
                           0.5 * dphi_dx[0][i] * dphiC_dx[0][j]) * k2;

            lagrMultMatrix[2*i][2*j] += H1xx * WJ;
            lagrMultMatrix[2*i+1][2*j] += H1yx * WJ;
            lagrMultMatrix[2*i][2*j+1] += H1xy * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += H1yy * WJ;

             // Arlequin Stabilization term
            double LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ * k1/ dens_;

            arlequinStab[2*i][2*j] += LL * WJ;
            arlequinStab[2*i+1][2*j+1] += LL * WJ; 

        };


        // L2 operator - Velocity
        double l2ux_ = -phi_[i] * una_[0] * k1;
        double l2uy_ = -phi_[i] * una_[1] * k1;

        // H1 operator - Velocity
        double h1ux_ = -(dphi_dx[0][i] * duna_dx[0][0] +
                        0.5 * dphi_dx[1][i] * duna_dx[0][1] +
                        0.5 * dphi_dx[1][i] * duna_dx[1][0]) * k2;

        double h1uy_ = -(0.5 * dphi_dx[0][i] * duna_dx[0][1] +
                        dphi_dx[1][i] * duna_dx[1][1] +
                        0.5 * dphi_dx[0][i] * duna_dx[1][0]) * k2;

        rhsVector2[2*i] += (l2ux_ + h1ux_) * WJ;
        rhsVector2[2*i+1] += (l2uy_ + h1uy_) * WJ;


        // L2 operator - Lagrange
        double l2x_ = -phiC_[i] * lambda_[0] * k1;
        double l2y_ = -phiC_[i] * lambda_[1] * k1;

        // H1 operator - Lagrange
        double h1x_ = -(dphiC_dx[0][i] * dlambda_dx[0][0] +
                       0.5 * dphiC_dx[1][i] * dlambda_dx[0][1] +
                       0.5 * dphiC_dx[1][i] * dlambda_dx[1][0]) * k2;

        double h1y_ = -(0.5 * dphiC_dx[0][i] * dlambda_dx[0][1] +
                       dphiC_dx[1][i] * dlambda_dx[1][1] +
                       0.5 * dphiC_dx[0][i] * dlambda_dx[1][0]) * k2;

        rhsVector1[2*i] += (l2x_ + h1x_) * WJ;
        rhsVector1[2*i+1] += (l2y_ + h1y_) * WJ;

        //Arlequin Stabilization vector
        double LLx = (dphi_dx[0][i] * dlambda_dx[0][0] + dphi_dx[1][i] * dlambda_dx[0][1]) * tARLQ_ * k1/ dens_;
        double LLy = (dphi_dx[0][i] * dlambda_dx[1][0] + dphi_dx[1][i] * dlambda_dx[1][1]) * tARLQ_ * k1/ dens_;

        arlequinStabVector[2*i] -= LLx * WJ;
        arlequinStabVector[2*i+1] -= LLy * WJ;

    };

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] dlambda_dx[i];
    delete[] dlambda_dx;
   
};



template <>
void Element<2>::getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG(double &djac_, double &weight_, double &tSUPG_, double &tPSPG_,
                                                              int &na, double *phi_, int &naC, double *phiC_, double **dphiC_dx,
                                                              double **jacobianNRMatrix, double *rhsVector)
{

    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();

    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVelCoarse(naC, phiC_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVelCoarse(naC, phiC_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    double lambda_[2];
    getInterpLambda(na, phi_, lambda_);

    double WJ = weight_ * djac_;

    for (int i = 0; i < naC; i++)
    {
        for (int j = 0; j < na; j++)
        {

            double msupg = ((una_[0] - umeshna_[0]) * dphiC_dx[0][i] + (una_[1] - umeshna_[1]) * dphiC_dx[1][i]) * phi_[j] * tSUPG_ * k1;
            jacobianNRMatrix[2*i][2*j] += msupg * WJ;
            jacobianNRMatrix[2*i+1][2*j+1] += msupg * WJ;

            double mpspgx = (dphiC_dx[0][i] * phi_[j]) * tPSPG_ * k1/ dens_;
            double mpspgy = (dphiC_dx[1][i] * phi_[j]) * tPSPG_ * k1/ dens_;

            jacobianNRMatrix[2*naC+i][2*j] += mpspgx * WJ;
            jacobianNRMatrix[2*naC+i][2*j+1] += mpspgy * WJ;
        };

        double vsupgx = -((una_[0] - umeshna_[0]) * dphiC_dx[0][i] + (una_[1] - umeshna_[1]) * dphiC_dx[1][i]) * lambda_[0] * tSUPG_ * k1;
        double vsupgy = -((una_[0] - umeshna_[0]) * dphiC_dx[0][i] + (una_[1] - umeshna_[1]) * dphiC_dx[1][i]) * lambda_[1] * tSUPG_ * k1;

        double vpspg = -(dphiC_dx[0][i] * lambda_[0] + dphiC_dx[1][i] * lambda_[1]) * tPSPG_ * k1/ dens_;

        rhsVector[2*i] += vsupgx * WJ;
        rhsVector[2*i+1] += vsupgy * WJ;
        rhsVector[2*naC+i] += vpspg * WJ;
    };
};

template <>
void Element<2>::getMatrixAndVectorsDifferentMeshArlqStab(double &djac_, double &weight_, double &tARLQ_, int &index,
                                                          int &na, double **dphi_dx, int &naC, double *phiC_, double **dphiC_dx, double ***ddphiC_dx,
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
    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVelCoarse(naC, phiC_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVelCoarse(naC, phiC_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDerCoarse(naC, dphiC_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    //Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMesh_dx[i] = new double[2];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMeshPrev_dx[i] = new double[2];

    double dumeshna_dx[2][2];
    getInterpMeshVelDerCoarse(naC, dphiC_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];
        };
    };

    //Acceleration derivatives
    double **daccel_dx;
    daccel_dx = new double *[2];
    for (int i = 0; i < 2; ++i) daccel_dx[i] = new double[2];

    double **daccelPrev_dx;
    daccelPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) daccelPrev_dx[i] = new double[2];

    double daccelm_dx[2][2];
    getInterpAccelDerCoarse(naC, dphiC_dx, daccel_dx, daccelPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            daccelm_dx[i][j] = alpha_m * daccel_dx[i][j] + (1. - alpha_m) * daccelPrev_dx[i][j];
        };
    };

    //Second pressure derivatives
    double **ddpress_dxdx;
    ddpress_dxdx = new double *[2];
    for (int i = 0; i < 2; ++i) ddpress_dxdx[i] = new double[2];
    getInterpSecondPressDerCoarse(naC, ddphiC_dx, ddpress_dxdx);


    //Second velocity derivatives
    double ***ddu_dxdx;
    ddu_dxdx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddu_dxdx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
        {
            ddu_dxdx[i][j] = new double[2];
        };
    };

    double ***dduPrev_dxdx;
    dduPrev_dxdx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        dduPrev_dxdx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
        {
            dduPrev_dxdx[i][j] = new double[2];
        };
    };

    double dduna_dxdx[2][2][2];
    getInterpSecondVelDerCoarse(naC, ddphiC_dx, ddu_dxdx, dduPrev_dxdx);
    
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 2; k++){
                dduna_dxdx[i][j][k] = alpha_f * ddu_dxdx[i][j][k] + (1. - alpha_f) * dduPrev_dxdx[i][j][k];
            };
        };
    };

    //Lagrange Multiplieres derivatives
    double **dlambda_dx;
    dlambda_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dlambda_dx[i] = new double[2];    
    getInterpLambdaDer(na, dphi_dx, dlambda_dx);

    double wna_;
    if (na == 6) wna_ = alpha_f * intPointWeightFunctionSpecial_FEM[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev_FEM[index];
    if (na == 9) wna_ = alpha_f * intPointWeightFunctionSpecial_ISO[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev_ISO[index];
    
    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < naC; j++)
        {

            // tarlq x mass matrix
            double mass = -(dphi_dx[0][i] * dphiC_dx[0][j] + dphi_dx[1][i] * dphiC_dx[1][j]) * tARLQ_ ;

            arlequinStab0[2*i][2*j] += mass * WJ * alpha_m;
            arlequinStab0[2*i+1][2*j+1] += mass * WJ * alpha_m;

            // tarlq x convecction matrix
            double convec1 = -((dphi_dx[0][i] * (ddphiC_dx[0][0][j] * (una_[0] - umeshna_[0]) + ddphiC_dx[0][1][j] * (una_[1] - umeshna_[1]))) + 
                              (dphi_dx[1][i] * (ddphiC_dx[1][0][j] * (una_[0] - umeshna_[0]) + ddphiC_dx[1][1][j] * (una_[1] - umeshna_[1])))) * tARLQ_;
            double convec2 = -((dphi_dx[0][i] * (dphiC_dx[0][j] * (duna_dx[0][0] - dumeshna_dx[0][0]) + dphiC_dx[1][j] * (duna_dx[1][0] - dumeshna_dx[1][0]))) + 
                              (dphi_dx[1][i] * (dphiC_dx[0][j] * (duna_dx[0][1] - dumeshna_dx[0][1]) + dphiC_dx[1][j] * (duna_dx[1][1] - dumeshna_dx[1][1])))) * tARLQ_;

            arlequinStab0[2*i][2*j] += (convec1 + convec2) * WJ * alpha_f * gamma * dTime_;
            arlequinStab0[2*i+1][2*j+1] += (convec1 + convec2) * WJ * alpha_f * gamma * dTime_;

            // tarlq x pressure term
            double pressx = -(dphi_dx[0][i]* ddphiC_dx[0][0][j] + dphi_dx[1][i]* ddphiC_dx[0][1][j]) *  tARLQ_/dens_;
            double pressy = -(dphi_dx[0][i]* ddphiC_dx[1][0][j] + dphi_dx[1][i]* ddphiC_dx[1][1][j]) *  tARLQ_/dens_;
            arlequinStab0[2*i][2*naC+j] += pressx * WJ;
            arlequinStab0[2*i+1][2*naC+j] += pressy * WJ;
        };

        // tarlq x mass matrix
        double massx = (dphi_dx[0][i] * daccelm_dx[0][0] + dphi_dx[1][i] * daccelm_dx[0][1]) * tARLQ_;
        double massy = (dphi_dx[0][i] * daccelm_dx[1][0] + dphi_dx[1][i] * daccelm_dx[1][1]) * tARLQ_;

        arlequinStabVector0[2*i] += massx * WJ;
        arlequinStabVector0[2*i+1] += massy * WJ;

        // tarlq x convecction matrix
        double convec1x = ((dphi_dx[0][i] * (dduna_dxdx[0][0][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[0][0][1] * (una_[1] - umeshna_[1]))) + 
                           (dphi_dx[1][i] * (dduna_dxdx[0][1][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[0][1][1] * (una_[1] - umeshna_[1])))) * tARLQ_;
        double convec2x = ((dphi_dx[0][i] * (duna_dx[0][0] * (duna_dx[0][0] - dumeshna_dx[0][0]) + duna_dx[0][1] * (duna_dx[1][0] - dumeshna_dx[1][0]))) + 
                           (dphi_dx[1][i] * (duna_dx[0][0] * (duna_dx[0][1] - dumeshna_dx[0][1]) + duna_dx[0][1] * (duna_dx[1][1] - dumeshna_dx[1][1]))))* tARLQ_;

        double convec1y = ((dphi_dx[0][i] * (dduna_dxdx[1][0][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[1][0][1] * (una_[1] - umeshna_[1]))) + 
                           (dphi_dx[1][i] * (dduna_dxdx[1][1][0] * (una_[0] - umeshna_[0]) + dduna_dxdx[1][1][1] * (una_[1] - umeshna_[1]))))* tARLQ_;
        double convec2y = ((dphi_dx[0][i] * (duna_dx[1][0] * (duna_dx[0][0] - dumeshna_dx[0][0]) + duna_dx[1][1] * (duna_dx[1][0] - dumeshna_dx[1][0]))) + 
                           (dphi_dx[1][i] * (duna_dx[1][0] * (duna_dx[0][1] - dumeshna_dx[0][1]) + duna_dx[1][1] * (duna_dx[1][1] - dumeshna_dx[1][1]))))* tARLQ_;

        arlequinStabVector0[2*i] += (convec1x + convec2x) * WJ;
        arlequinStabVector0[2*i+1] += (convec1y + convec2y) * WJ;

        //tarlq x pressure term
        double pressx = (dphi_dx[0][i]* ddpress_dxdx[0][0] + dphi_dx[1][i]* ddpress_dxdx[1][0]) * tARLQ_/dens_;
        double pressy = (dphi_dx[0][i]* ddpress_dxdx[0][1] + dphi_dx[1][i]* ddpress_dxdx[1][1]) * tARLQ_/dens_;
        arlequinStabVector0[2*i] += pressx * WJ;
        arlequinStabVector0[2*i+1] += pressy * WJ;
    };

    // Arlequin Diagonal Stabilization terms
    for (int i = 0; i < na; i++)
    {
        for (int j = 0; j < na; j++)
        {

            double LL = -(dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ * k1/ dens_;

            arlequinStabD[2*i][2*j] += LL * WJ/wna_;
            arlequinStabD[2*i+1][2*j+1] += LL * WJ/wna_;
        };

        // Stabilization term
        double LLx = (dphi_dx[0][i] * dlambda_dx[0][0] + dphi_dx[1][i] * dlambda_dx[0][1]) * tARLQ_ * k1/ dens_;
        double LLy = (dphi_dx[0][i] * dlambda_dx[1][0] + dphi_dx[1][i] * dlambda_dx[1][1]) * tARLQ_ * k1/ dens_;

        arlequinStabVectorD[2*i] += LLx * WJ/wna_;
        arlequinStabVectorD[2*i+1] += LLy * WJ/wna_;
    };

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < 2; ++i) delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] daccel_dx[i];
    delete[] daccel_dx;
    
    for (int i = 0; i < 2; ++i) delete[] daccelPrev_dx[i];
    delete[] daccelPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] dlambda_dx[i];
    delete[] dlambda_dx;

    for (int i = 0; i < 2; ++i) delete[] ddpress_dxdx[i];
    delete[] ddpress_dxdx;


    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] ddu_dxdx[i][j];
        };
        delete[] ddu_dxdx[i];
    };
    delete[] ddu_dxdx;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
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
template <>
void Element<2>::getResidualVector(int &na,double &djac_, double &weight_, double &tSUPG_, double &tPSPG_, double &tLSIC_,
                                   int &index, double *phi_, double **dphi_dx, double *rhsVector)
{

    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();

    double ff[2];
    ff[0] = parameters->getFieldForce(0);
    ff[1] = parameters->getFieldForce(1);

    // Interpolated variables

    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVel(na, phi_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVel(na, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Acceleration
    double accel_[2], accelPrev_[2], accelm_[2];
    getInterpAccel(na, phi_, accel_, accelPrev_);
    for (int i = 0; i < 2; i++) accelm_[i] = alpha_m * accel_[i] + (1. - alpha_m) * accelPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    //Pressure
    double press_;
    getInterpPress(na, phi_, press_);

    //Pressure derivatives
    double dpress_dx[2];
    getInterpPressDer(na, dphi_dx, dpress_dx);

    double wna_;
    if (na == 6) wna_ = alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];
    if (na == 9) wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];

    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < na; i++)
    {

        double mx = phi_[i] * accelm_[0] * dens_ +
                    ((una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i]) * accelm_[0] * tSUPG_ * dens_;
        double my = phi_[i] * accelm_[1] * dens_ +
                    ((una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i]) * accelm_[1] * tSUPG_ * dens_;

        double Kx = 2. * dphi_dx[0][i] * duna_dx[0][0] * visc_ +
                    dphi_dx[1][i] * duna_dx[0][1] * visc_ +
                    dphi_dx[1][i] * duna_dx[1][0] * visc_;

        double Ky = dphi_dx[0][i] * duna_dx[0][1] * visc_ +
                    2. * dphi_dx[1][i] * duna_dx[1][1] * visc_ +
                    dphi_dx[0][i] * duna_dx[1][0] * visc_;

        double KLSx = dphi_dx[0][i] * (duna_dx[0][0] + duna_dx[1][1]) * tLSIC_ * dens_;
        double KLSy = dphi_dx[1][i] * (duna_dx[0][0] + duna_dx[1][1]) * tLSIC_ * dens_;

        double Cx = (duna_dx[0][0] * (una_[0] - umeshna_[0]) + duna_dx[0][1] * (una_[1] - umeshna_[1])) * phi_[i] * dens_ +
                    ((una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i]) *
                    ((una_[0] - umeshna_[0]) * duna_dx[0][0] + (una_[1] - umeshna_[1]) * duna_dx[0][1]) * tSUPG_ * dens_;
        double Cy = (duna_dx[1][0] * (una_[0] - umeshna_[1]) + duna_dx[1][1] * (una_[1] - umeshna_[1])) * phi_[i] * dens_ +
                    ((una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[0]) * dphi_dx[1][i]) *
                    ((una_[0] - umeshna_[0]) * duna_dx[1][0] + (una_[1] - umeshna_[1]) * duna_dx[1][1]) * tSUPG_ * dens_;

        double Px = -(dphi_dx[0][i] * press_) +
                    ((dphi_dx[0][i] * (una_[0] - umeshna_[0]) + dphi_dx[1][i] * (una_[1] - umeshna_[0])) * dpress_dx[0] * tSUPG_);
        double Py = -(dphi_dx[1][i] * press_) +
                    ((dphi_dx[0][i] * (una_[0] - umeshna_[0]) + dphi_dx[1][i] * (una_[1] - umeshna_[1])) * dpress_dx[1] * tSUPG_);

        double Q = ((duna_dx[0][0] + duna_dx[1][1]) * phi_[i]) +
                   (dphi_dx[0][i] * dpress_dx[0] + dphi_dx[1][i] * dpress_dx[1]) * tPSPG_ / dens_ +
                   dphi_dx[0][i] * ((una_[0] - umeshna_[0]) * duna_dx[0][0] + (una_[1] - umeshna_[1]) * duna_dx[0][1]) * tPSPG_ +
                   dphi_dx[1][i] * ((una_[0] - umeshna_[0]) * duna_dx[1][0] + (una_[1] - umeshna_[1]) * duna_dx[1][1]) * tPSPG_ +
                   dphi_dx[0][i] * accelm_[0] * tPSPG_ +
                   dphi_dx[1][i] * accelm_[1] * tPSPG_;

        double Ffvx = phi_[i] * dens_ * ff[0] +
                      tSUPG_ * dens_ * ((una_[0] - umeshna_[1]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i]) * ff[0];
        double Ffvy = phi_[i] * dens_ * ff[1] +
                      tSUPG_ * dens_ * ((una_[0] - umeshna_[0]) * dphi_dx[0][i] + (una_[1] - umeshna_[1]) * dphi_dx[1][i]) * ff[1];

        double Ffp = tPSPG_ * (dphi_dx[0][i] * ff[0] + dphi_dx[1][i] * ff[1]);

        rhsVector[2*i] += (-mx + Ffvx - Kx - Px - Cx - KLSx) * WJ;
        rhsVector[2*i+1] += (-my + Ffvy - Ky - Py - Cy - KLSy) * WJ;
        rhsVector[2*na+i] += (Ffp - Q) * WJ;
    };

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;
    
    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    return;
};



//------------------------------------------------------------------------------
//----------------------ARLEQUIN ELEMENT LOCAL MATRIX/VECTOR--------------------
//------------------------------------------------------------------------------
// template<>
// void Element<2>::getLagrangeMultipliersSameMesh_FEM(double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){

//     //quadrature and functions local classes
//     NormalQuad              nQuad = NormalQuad();
//     QuadShapeFunction<2>    shapeQuad;

//     //data for computation of IGA basis functions
//     double xsi[2], phi_[6];

//     double **dphi_dx;
//     dphi_dx = new double*[2];
//     for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

//     double **ainv_;
//     ainv_ = new double*[2];
//     for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

//     double **Jac;
//     Jac = new double*[2];
//     for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

//     int index = 0;
//     for(double* it = nQuad.beginFem(); it != nQuad.endFem(); it++){

//         //Defines the integration points adimentional coordinates
//         xsi[0] = nQuad.PointListFem(index,0);
//         xsi[1] = nQuad.PointListFem(index,1);
//         //Returns the quadrature integration weight
//         weight_ = nQuad.WeightListFem(index);

//         //Computes the velocity shape functions
//         shapeQuad.evaluateFem(xsi,phi_);

//         //Computes the jacobian matrix
//         getJacobianMatrix_FEM(xsi,Jac,ainv_);

//         //Computes spatial derivatives
//         getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

//         //Interpolates variables
//         getInterpolatedVariablesSameMesh_FEM(phi_,dphi_dx);

//         //Computes matrixes and vectors
//         getMatrixAndVectorsSameMesh_FEM(phi_,dphi_dx,lagrMultMatrix,
//                                         rhsVector1,rhsVector2);

//         index++;

//     };

//     for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
//     delete [] dphi_dx;
//     for (int i = 0; i < 2; ++i) delete [] ainv_[i];
//     delete [] ainv_;
//     for (int i = 0; i < 2; ++i) delete [] Jac[i];
//     delete [] Jac;

//     return;

// };

template <>
void Element<2>::getLagrangeMultipliersSameMesh_FEM(int &index, double **lagrMultMatrix,
                                                    double *rhsVector1, double *rhsVector2)

{
    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    double djac_, weight_;

    // data for computation of IGA basis functions
    double xsi[2], phi_[6];

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double *[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

    // Defines the integration points adimentional coordinates
    xsi[0] = nQuad.PointListFem(index, 0);
    xsi[1] = nQuad.PointListFem(index, 1);
    // Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    // Computes the velocity shape functions
    shapeQuad.evaluateFem(xsi, phi_);

    // Computes the jacobian matrix
    getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

    // Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

    // Computes matrixes and vectors
    int na = 6;
    getMatrixAndVectorsSameMesh(na, djac_, weight_, phi_, dphi_dx, lagrMultMatrix,
                                rhsVector1, rhsVector2);

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;
    
    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;
    
    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersSameMesh_tSUPG_tPSPG_FEM(int &index, double **jacobianNRMatrix, double *rhsVector)
{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // data for computation of IGA basis functions
    double xsi[2], phi_[6];

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
    double djac_, weight_;

    // Defines the integration points adimentional coordinates
    xsi[0] = nQuad.PointListFem(index, 0);
    xsi[1] = nQuad.PointListFem(index, 1);
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
    int na = 6;
    getMatrixAndVectorsSameMesh_tSUPG_tPSPG(na,djac_, weight_, tSUPG_, tPSPG_, phi_, dphi_dx,
                                            jacobianNRMatrix, rhsVector);

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
}

// template <>
// void Element<2>::getLagrangeMultipliersSameMesh_tSUPG_tPSPG_FEM(double **jacobianNRMatrix, double *rhsVector)
// {

//     // quadrature and functions local classes
//     NormalQuad nQuad = NormalQuad();
//     QuadShapeFunction<2> shapeQuad;

//     // data for computation of IGA basis functions
//     double xsi[2], phi_[6];

//     double **dphi_dx;
//     dphi_dx = new double *[2];
//     for (int i = 0; i < 2; ++i)
//         dphi_dx[i] = new double[6];

//     double **ainv_;
//     ainv_ = new double *[2];
//     for (int i = 0; i < 2; ++i)
//         ainv_[i] = new double[2];

//     double **Jac;
//     Jac = new double *[2];
//     for (int i = 0; i < 2; ++i)
//         Jac[i] = new double[2];

//     int index = 0;
//     for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
//     {

//         // Defines the integration points adimentional coordinates
//         xsi[0] = nQuad.PointListFem(index, 0);
//         xsi[1] = nQuad.PointListFem(index, 1);
//         // Returns the quadrature integration weight
//         weight_ = nQuad.WeightListFem(index);

//         // Computes the velocity shape functions
//         shapeQuad.evaluateFem(xsi, phi_);

//         // Computes the jacobian matrix
//         getJacobianMatrix_FEM(xsi, Jac, ainv_);

//         // Computes spatial derivatives
//         getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

//         // Computes stabilization parameters
//         getNewParameterSUPG_FEM(Jac, phi_, dphi_dx);

//         // Computes matrixes and vectors
//         getMatrixAndVectorsSameMesh_tSUPG_tPSPG_FEM(phi_, dphi_dx, jacobianNRMatrix, rhsVector);

//         index++;
//     };

//     for (int i = 0; i < 2; ++i)
//         delete[] dphi_dx[i];
//     delete[] dphi_dx;
//     for (int i = 0; i < 2; ++i)
//         delete[] ainv_[i];
//     delete[] ainv_;
//     for (int i = 0; i < 2; ++i)
//         delete[] Jac[i];
//     delete[] Jac;

//     return;
// }

template <>
void Element<2>::getLagrangeMultipliersSameMeshArlqStab_FEM(int &index, int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                            std::vector<IParameters *> &iparamC, double **arlequinStabD,
                                                            double *arlequinStabVectorD, double **arlequinStab1, double *arlequinStabVector1)
{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for IGA coarse mesh computation
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9], xsiC[2], phiC_[9];
    for (int i = 0; i < 9; i++) wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]]->getINC();

    // data for computation of FEM fine mesh computations
    double xsi[2], phi_[6];

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double ***ddphi_dx;
    ddphi_dx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddphi_dx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
            ddphi_dx[i][j] = new double[6];
    };

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double *[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

    double tARLQ_, djac_, weight_;

    // Defines the integration points adimentional coordinates
    xsi[0] = nQuad.PointListFem(index, 0);
    xsi[1] = nQuad.PointListFem(index, 1);
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
    for (int k = 0; k < 2; k++)
        xsiC[k] = intPointCorrespXsi_FEM[index][k];

    // Computes the coqrse shape functions
    shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

    
    // Computes tARLQ_
    //old way
     // getParameterArlequin_FEM(tARLQ_,phi_,dphi_dx);
   
    //Vector norm
    // getParameterArlequin2_FEM(index, djac_, weight_, tARLQ_, phi_, phiC_, dphi_dx, ddphi_dx);
   
    //matriz norm
    getParameterArlequinMN_FEM(index, djac_, weight_, tARLQ_, phi_, phiC_, dphi_dx, ddphi_dx);

    //by Elemnt matriz norm
    // tARLQ_ = tARLQEL_;

    // Computes matrixes and vectors
    int na = 6;
    getMatrixAndVectorsSameMeshArlqStab(na,djac_, weight_, tARLQ_, index, phi_, dphi_dx, ddphi_dx,
                                        arlequinStabD, arlequinStabVectorD,
                                        arlequinStab1, arlequinStabVector1);


    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
            delete[] ddphi_dx[i][j];
        delete[] ddphi_dx[i];
    };
    delete[] ddphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;
};

// template<>
// void Element<2>::getLagrangeMultipliersSameMeshArlqStab_FEM(double **arlequinStabD, double *arlequinStabVectorD,
//                                                             double **arlequinStab1, double *arlequinStabVector1){

//     //quadrature and functions local classes
// NormalQuad              nQuad = NormalQuad();
// QuadShapeFunction<2>    shapeQuad;

// //data for computation of IGA basis functions
// double xsi[2], phi_[6];

// double **dphi_dx;
// dphi_dx = new double*[2];
// for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

// double ***ddphi_dx;
// ddphi_dx = new double**[2];
// for (int i = 0; i < 2; ++i) {
//     ddphi_dx[i] = new double*[2];
//     for (int j = 0; j < 2; j++) ddphi_dx[i][j] = new double[6];
// };

// double **ainv_;
// ainv_ = new double*[2];
// for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

// double **Jac;
// Jac = new double*[2];
// for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

// int index = 0;
// for(double* it = nQuad.beginFem(); it != nQuad.endFem(); it++){

//     //Defines the integration points adimentional coordinates
//     xsi[0] = nQuad.PointListFem(index,0);
//     xsi[1] = nQuad.PointListFem(index,1);
//     //Returns the quadrature integration weight
//     weight_ = nQuad.WeightListFem(index);

//     //Computes the velocity shape functions
//     shapeQuad.evaluateFem(xsi,phi_);

//     //Computes the jacobian matrix
//     getJacobianMatrix_FEM(xsi,Jac,ainv_);

//     //Computes spatial derivatives
//     getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

//     //Computes spatial second derivatives
//     getSecondSpatialDerivatives_FEM(xsi,ainv_,ddphi_dx);

//     //Interpolates variables
//     getInterpolatedVariablesSameMeshArlqStab_FEM(phi_,dphi_dx,ddphi_dx);

//     // getParameterArlequin2_FEM(phi_,dphi_dx,ddphi_dx);
//     // //get Arlequin stabilization parameter
//     // if (iTimeStep < 1) getParameterArlequin_FEM(phi_,dphi_dx);

//     getParameterArlequin_FEM(phi_,dphi_dx);

//     //Computes matrixes and vectors
//     getMatrixAndVectorsSameMeshArlqStab_FEM(index,phi_,dphi_dx,ddphi_dx,arlequinStabD,arlequinStabVectorD,
//                                             arlequinStab1,arlequinStabVector1);

//     index++;

// };

// for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
// delete [] dphi_dx;

// for (int i = 0; i < 2; ++i) {
//     for (int j = 0; j < 2; ++j) delete [] ddphi_dx[i][j];
//     delete [] ddphi_dx[i];
// };
// delete [] ddphi_dx;

// for (int i = 0; i < 2; ++i) delete [] ainv_[i];
// delete [] ainv_;
// for (int i = 0; i < 2; ++i) delete [] Jac[i];
// delete [] Jac;

// return;

// };


template<>
void Element<2>::getParameterArlequinElem(double &tarlq, int &index, int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                          std::vector<IParameters *> &iparamC){

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for IGA coarse mesh computation
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9], xsiC[2], phiC_[9];
    for (int i = 0; i < 9; i++) wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]]->getINC();

    // data for computation of FEM fine mesh computations
    double xsi[2], phi_[6], djac_;

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double ***ddphi_dx;
    ddphi_dx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddphi_dx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
            ddphi_dx[i][j] = new double[6];
    };

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double *[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

    // Defines the integration points adimentional coordinates
    xsi[0] = nQuad.PointListFem(index, 0);
    xsi[1] = nQuad.PointListFem(index, 1);

    // Computes the shape functions
    shapeQuad.evaluateFem(xsi, phi_);

    // Computes the jacobian matrix
    getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

    // Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

    // Computes spatial second derivatives
    getSecondSpatialDerivatives_FEM(ainv_, ddphi_dx);

    // Integration point in the coarse mesh element
    for (int k = 0; k < 2; k++)
        xsiC[k] = intPointCorrespXsi_FEM[index][k];

    // Computes the coqrse shape functions
    shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);


    //computing tARLQ

    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters-> getGamma();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &dTime_ = parameters->getTimeStep();    
    double &k1 = parameters->getArlequinK1();

    // Interpolated variables
    int na = 6;

    //Velocity
    double u_[2], uPrev_[2], una_[2];
    getInterpVel(na, phi_, u_, uPrev_);
    for (int i = 0; i < 2; i++) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Mesh Velocity
    double uMesh_[2], uMeshPrev_[2], umeshna_[2];
    getInterpMeshVel(na, phi_, uMesh_, uMeshPrev_);
    for (int i = 0; i < 2; i++) umeshna_[i] = alpha_f * uMesh_[i] + (1. - alpha_f) * uMeshPrev_[i];

    //Velocity Derivatives
    double **du_dx;
    du_dx = new double *[2];
    for (int i = 0; i < 2; ++i) du_dx[i] = new double[2];

    double **duPrev_dx;
    duPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duPrev_dx[i] = new double[2];

    double duna_dx[2][2];
    getInterpVelDer(na, dphi_dx, du_dx, duPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duPrev_dx[i][j];
        };
    };

    //Mesh velocity derivatives
    double **duMesh_dx;
    duMesh_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMesh_dx[i] = new double[2];

    double **duMeshPrev_dx;
    duMeshPrev_dx = new double *[2];
    for (int i = 0; i < 2; ++i) duMeshPrev_dx[i] = new double[2];

    double dumeshna_dx[2][2];
    getInterpMeshVelDer(na, dphi_dx, duMesh_dx, duMeshPrev_dx);
    for (int i = 0; i <2; i++){
        for (int j = 0; j < 2; j++){
            dumeshna_dx[i][j] = alpha_f * duMesh_dx[i][j] + (1. - alpha_f) * duMeshPrev_dx[i][j];
        };
    };

    double lambda1[12][12] = {};
    double lambda0[12][18] = {};
    double convec[12][12] = {};
    double iner[12][12] = {};
    double visc[12][12] = {};
    double press[12][12] = {};
    double gradlambda[12][12] = {};

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {

           //inercia matrix
            iner[2*i][2*j] += (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * alpha_m;
            iner[2*i+1][2*j+1] += (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * alpha_m;

            //convecction matrix
            double convec1 = (dphi_dx[0][i] * ((una_[0] - umeshna_[0]) * ddphi_dx[0][0][j] + (duna_dx[0][0] - dumeshna_dx[0][0]) * dphi_dx[0][j] + 
                                               (una_[1] - umeshna_[1]) * ddphi_dx[1][0][j] + (duna_dx[1][0] - dumeshna_dx[1][0]) * dphi_dx[1][j]) +  
                              dphi_dx[1][i] * ((una_[0] - umeshna_[0]) * ddphi_dx[0][1][j] + (duna_dx[0][1] - dumeshna_dx[0][1]) * dphi_dx[0][j] + 
                                               (una_[1] - umeshna_[1]) * ddphi_dx[1][1][j] + (duna_dx[1][1] - dumeshna_dx[1][1]) * dphi_dx[1][j]));

            convec[2*i][2*j] += convec1 * alpha_f * gamma * dTime_;
            convec[2*i+1][2*j+1] += convec1 * alpha_f * gamma * dTime_;

            //pressure term
            press[2*i][j] += (dphi_dx[0][i] * ddphi_dx[0][0][j] + dphi_dx[1][i] * ddphi_dx[1][0][j])/dens_;
            press[2*i+1][j] += (dphi_dx[0][i] * ddphi_dx[0][1][j] + dphi_dx[1][i] * ddphi_dx[1][1][j])/dens_;

            gradlambda[2*i][2*j] += (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j])/dens_*k1;
            gradlambda[2*i+1][2*j+1] += (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j])/dens_*k1;
        
            // viscosity
            visc[2*i][2*j]     += (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * ddphi_dx[0][0][j] + ddphi_dx[1][1][j]) * visc_ * alpha_f * gamma * dTime_;
            visc[2*i][2*j+1]   += (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * ddphi_dx[0][1][j] * visc_ * alpha_f * gamma * dTime_;
            visc[2*i+1][2*j]   += (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * ddphi_dx[1][0][j] * visc_ * alpha_f * gamma * dTime_;
            visc[2*i+1][2*j+1] += (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * ddphi_dx[1][1][j] + ddphi_dx[0][0][j]) * visc_ * alpha_f * gamma * dTime_;

            // lagrange multipliers
            lambda1[2*i][2*j] += phi_[i] * phi_[j] * k1;
            lambda1[2*i+1][2*j+1] += phi_[i] * phi_[j] * k1;
        };

        for (int j = 0; j < 9; j++)
        {
            // lagrange multipliers
            lambda0[2*i][2*j] += phi_[i] * phiC_[j] * k1;
            lambda0[2*i+1][2*j+1] += phi_[i] * phiC_[j] * k1;
        };

    };

    double normLambda0 = 0.;
    double normLambda1 = 0.;
    double normConvec = 0.;
    double normIner = 0.;
    double normVisc = 0.;
    double normPress = 0.;
    double normGradLambda = 0.;

    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 12; j++){
            normLambda1 += lambda1[i][j] * lambda1[i][j];
            normConvec += convec[i][j] * convec[i][j];
            normIner += iner[i][j] * iner[i][j];
            normVisc += visc[i][j] * visc[i][j];
            normPress += press[i][j] * press[i][j];
            normGradLambda += gradlambda[i][j] * gradlambda[i][j];
        }
        for (int j = 0; j < 18; j++){
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
    if(tE1>1.)tE1=0.0000000001;
    if(tE2>1.)tE2=0.0000000001;

    double tA = 1. / sqrt(1. / (tA1 * tA1) + 1. / (tA2 * tA2));

    double tB = 1. / sqrt(1. / (tB1 * tB1) + 1. / (tB2 * tB2));

    double tC = 1. / sqrt(1. / (tC1 * tC1) + 1. / (tC2 * tC2));

    double tD = 1. / sqrt(1. / (tD1 * tD1) + 1. / (tD2 * tD2));

    double tE = 1. / sqrt(1. / (tE1 * tE1) + 1. / (tE2 * tE2));

    tarlq = 1. / sqrt(1. / (tA * tA) +
                        1. / (tB * tB) +
                        1. / (tC * tC) +
                        1. / (tE * tE) +
                        1. / (tD * tD));

   

    //deallocating memory
    for (int i = 0; i < 2; ++i) delete[] du_dx[i];
    delete[] du_dx;

    for (int i = 0; i < 2; ++i) delete[] duPrev_dx[i];
    delete[] duPrev_dx;

    for (int i = 0; i < 2; ++i) delete[] duMesh_dx[i];
    delete[] duMesh_dx;

    for (int i = 0; i < 2; ++i) delete[] duMeshPrev_dx[i];
    delete[] duMeshPrev_dx;



    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
            delete[] ddphi_dx[i][j];
        delete[] ddphi_dx[i];
    };
    delete[] ddphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;                                        

};

template <>
void Element<2>::getLagrangeMultipliersSameMesh_ISO(double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // data for computation of IGA basis functions
    double wpc[9], xsi[2], phi_[9];
    for (int i = 0; i < 9; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    double djac_, weight_;

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i)
        dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i)
        ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double *[2];
    for (int i = 0; i < 2; ++i)
        quadJacMat[i] = new double[2];

    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        // Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListIso(index, 0);
        xsi[1] = nQuad.PointListIso(index, 1);
        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        // Computes the jacobian matrix
        getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

        // Computes matrixes and vectors
        int na = 9;
        getMatrixAndVectorsSameMesh(na, djac_, weight_, phi_, dphi_dx, lagrMultMatrix,
                                    rhsVector1, rhsVector2);

        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;
    
    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;
    
    for (int i = 0; i < 2; ++i) delete[] quadJacMat[i];
    delete[] quadJacMat;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersSameMesh_tSUPG_tPSPG_ISO(double **jacobianNRMatrix, double *rhsVector)
{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // data for computation of IGA basis functions
    double wpc[9], xsi[2], phi_[9];
    for (int i = 0; i < 9; i++)
        wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    double tSUPG_, tPSPG_, tLSIC_;
    double djac_, weight_;

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i)
        dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i)
        ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double *[2];
    for (int i = 0; i < 2; ++i)
        quadJacMat[i] = new double[2];

    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        // Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListIso(index, 0);
        xsi[1] = nQuad.PointListIso(index, 1);
        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        // Computes the jacobian matrix
        getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

        // Computes stabilization parameters
        getNewParameterSUPG_ISO(tSUPG_, tPSPG_, tLSIC_, quadJacMat, phi_, dphi_dx);

        // Computes matrixes and vectors
        int na = 9;
        getMatrixAndVectorsSameMesh_tSUPG_tPSPG(na,djac_, weight_, tSUPG_, tPSPG_, phi_, dphi_dx,
                                                jacobianNRMatrix, rhsVector);

        index++;
    };

    for (int i = 0; i < 2; ++i)
        delete[] dphi_dx[i];
    delete[] dphi_dx;
    for (int i = 0; i < 2; ++i)
        delete[] ainv_[i];
    delete[] ainv_;
    for (int i = 0; i < 2; ++i)
        delete[] quadJacMat[i];
    delete[] quadJacMat;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersSameMeshArlqStab_ISO(double **arlequinStabD, double *arlequinStabVectorD,
                                                            double **arlequinStab1, double *arlequinStabVector1)
{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // data for computation of IGA basis functions
    double wpc[9], xsi[2], phi_[9];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    double tARLQ_,djac_, weight_;

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double *[2];
    for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];

    double ***ddphi_dx;
    ddphi_dx = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddphi_dx[i] = new double *[2];
        for (int j = 0; j < 2; j++)
            ddphi_dx[i][j] = new double[9];
    };

    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        // Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListIso(index, 0);
        xsi[1] = nQuad.PointListIso(index, 1);
        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

        // Computes the jacobian matrix
        getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

        // Computes spatial second derivatives
        getSecondSpatialDerivatives_ISO(xsi, ainv_, ddphi_dx);

        // get Arlequin stabilization parameter
        getParameterArlequin_ISO(tARLQ_, phi_, dphi_dx);

        // Computes matrixes and vectors
        int na = 9;
        getMatrixAndVectorsSameMeshArlqStab(na, djac_, weight_, tARLQ_, index, phi_, dphi_dx, ddphi_dx,
                                            arlequinStabD, arlequinStabVectorD,
                                            arlequinStab1, arlequinStabVector1);

        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] quadJacMat[i];
    delete[] quadJacMat;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
            delete[] ddphi_dx[i][j];
        delete[] ddphi_dx[i];
    };
    delete[] ddphi_dx;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersDifferentMesh_FEM_FEM(std::vector<Nodes *> &nodesCoarse_, int *connecC, int &ielem,
                                                             double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                             double **arlequinStab, double *arlequinStabVector)
{

    int dim = 2;
    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for FEM coarse mesh computation (Velocity field)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;

    double xsiC[dim], phiC_[6];

    double **dphiC_dx;
    dphiC_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[6];

    double **ainvC_;
    ainvC_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **JacC;
    JacC = new double *[dim];
    for (int i = 0; i < dim; ++i) JacC[i] = new double[dim];

    // Data for FEM fine mesh computation (Lagrange field)
    double xsi[dim], phi_[6];

    double **dphi_dx;
    dphi_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **Jac;
    Jac = new double *[dim];
    for (int i = 0; i < dim; ++i) Jac[i] = new double[dim];

    double tARLQ_, djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k++) xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi, phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_, xsi, Jac, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi, ainv_, dphi_dx);

            // Computes Arlequin Stabilization term
            getParameterArlequin_FEM(tARLQ_, phi_, dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++)
                xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsiC, phiC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC, JacC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC, ainvC_, dphiC_dx);

            // Computes Matrixes and vectors
            int na = 6;
            getMatrixAndVectorsDifferentMeshTemp(djac_,weight_,tARLQ_,na,phi_,dphi_dx,na,phiC_,dphiC_dx,
                                                 lagrMultMatrix,rhsVector1,rhsVector2,
                                                 arlequinStab,arlequinStabVector);
        };
        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < 2; ++i) delete[] dphiC_dx[i];
    delete[] dphiC_dx;
   
    for (int i = 0; i < 2; ++i) delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < 2; ++i) delete[] JacC[i];
    delete[] JacC;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersDifferentMesh_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                             std::vector<IParameters *> &iparamC, int &ielem,
                                                             double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int dim = 2;
    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9], xsiC[dim], phiC_[9];
    for (int i = 0; i < 9; i++) wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]]->getINC();

    double **dphiC_dx;
    dphiC_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[9];

    double **ainvC_;
    ainvC_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **quadJacMatC;
    quadJacMatC = new double *[dim];
    for (int i = 0; i < dim; ++i) quadJacMatC[i] = new double[dim];

    // Data for FEM fine mesh computation (Lagrange field)
    double xsi[dim], phi_[6];

    double **dphi_dx;
    dphi_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **Jac;
    Jac = new double *[dim];
    for (int i = 0; i < dim; ++i) Jac[i] = new double[dim];

    double djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k++) xsi[k] = sQuad.PointListFem(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi,phi_);

            // Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_,xsi,Jac,ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC_,(*iparametersC),NpatchC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC,quadJacMatC,ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC,ainvC_,dphiC_dx);

            // Computes Matrix and vectors
            int na = 6;
            int naC = 9;
            getMatrixAndVectorsDifferentMesh(djac_, weight_, na, phi_, dphi_dx, naC, phiC_, dphiC_dx,
                                             lagrMultMatrix, rhsVector1, rhsVector2);
        };
        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < 2; ++i) delete[] dphiC_dx[i];
    delete[] dphiC_dx;
    
    for (int i = 0; i < 2; ++i) delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < 2; ++i) delete[] quadJacMatC[i];
    delete[] quadJacMatC;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                         std::vector<IParameters *> &iparamC, int &ielem,
                                                                         double **jacobianNRMatrix, double *rhsVector)
{

    int dim = 2;
    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9], xsiC[dim], phiC_[9];
    for (int i = 0; i < 9; i++) wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]]->getINC();

    double **dphiC_dx;
    dphiC_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[9];

    double **ainvC_;
    ainvC_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **quadJacMatC;
    quadJacMatC = new double *[dim];
    for (int i = 0; i < dim; ++i) quadJacMatC[i] = new double[dim];

    // Data for FEM fine mesh computation (Lagrange field)
    double xsi[dim], phi_[6];

    double **dphi_dx;
    dphi_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **Jac;
    Jac = new double *[dim];
    for (int i = 0; i < dim; ++i) Jac[i] = new double[dim];

    double tSUPG_, tPSPG_, tLSIC_, djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k++) xsi[k] = sQuad.PointListFem(index, k);

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
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC, quadJacMatC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC, ainvC_, dphiC_dx);

            // Computes Matrix and vectors
            int na = 6;
            int naC = 9;
            getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG(djac_, weight_, tSUPG_, tPSPG_, na, phi_, naC, 
                                                         phiC_, dphiC_dx,jacobianNRMatrix, rhsVector);
        };
        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < 2; ++i) delete[] dphiC_dx[i];
    delete[] dphiC_dx;
    
    for (int i = 0; i < 2; ++i) delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < 2; ++i) delete[] quadJacMatC[i];
    delete[] quadJacMatC;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersDifferentMeshArlqStab_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                                     std::vector<IParameters *> &iparamC, int &ielem,
                                                                     double **arlequinStabD, double *arlequinStabVectorD,
                                                                     double **arlequinStab0, double *arlequinStabVector0)
{

    int dim = 2;
    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9], xsiC[dim], phiC_[9];
    for (int i = 0; i < 9; i++) wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]]->getINC();

    double **dphiC_dx;
    dphiC_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[9];

    double ***ddphiC_dx;
    ddphiC_dx = new double **[dim];
    for (int i = 0; i < dim; ++i)
    {
        ddphiC_dx[i] = new double *[2];
        for (int j = 0; j < dim; j++)
            ddphiC_dx[i][j] = new double[9];
    };

    double **ainvC_;
    ainvC_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **quadJacMatC;
    quadJacMatC = new double *[dim];
    for (int i = 0; i < dim; ++i) quadJacMatC[i] = new double[dim];

    // Data for FEM fine mesh computation (Lagrange field)
    double xsi[dim], phi_[6];

    double **dphi_dx;
    dphi_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[6];

    double ***ddphi_dx;
    ddphi_dx = new double **[dim];
    for (int i = 0; i < dim; ++i)
    {
        ddphi_dx[i] = new double *[2];
        for (int j = 0; j < dim; j++)
            ddphi_dx[i][j] = new double[6];
    };

    double **ainv_;
    ainv_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **Jac;
    Jac = new double *[dim];
    for (int i = 0; i < dim; ++i) Jac[i] = new double[dim];

    double tARLQ_,djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginFem(); it != sQuad.endFem(); it++)
    {

        if ((intPointCorrespElem_FEM[index] == ielem))
        {

            //Fine mesh computation
            //Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k++) xsi[k] = sQuad.PointListFem(index, k);

            //Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            //Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi,phi_);

            //Computes the jacobian matrix
            getJacobianMatrix_FEM(djac_,xsi,Jac,ainv_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

            //Computes second spatial derivatives
            getSecondSpatialDerivatives_FEM(ainv_,ddphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

            //Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC_,(*iparametersC),NpatchC_);

            //Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC,quadJacMatC,ainvC_);

            //Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC,ainvC_,dphiC_dx);

            //Computes second spatial derivatives
            getSecondSpatialDerivatives_COARSE_ISO(xsiC,ainvC_,ddphiC_dx);

            //get Arlequin stabilization parameter
            getParameterArlequin2_COARSE_ISO(djac_,weight_,tARLQ_,phi_,phiC_,dphi_dx,dphiC_dx, 
                                             ddphi_dx,ddphiC_dx);

            // Computes Matrix and vectors
            int na = 6;
            int naC = 9;
            getMatrixAndVectorsDifferentMeshArlqStab(djac_,weight_,tARLQ_,index,na,dphi_dx,
                                                     naC,phiC_,dphiC_dx,ddphiC_dx,
                                                     arlequinStabD,arlequinStabVectorD,
                                                     arlequinStab0,arlequinStabVector0);
        };
        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;

    for (int i = 0; i < 2; ++i) delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            ddphiC_dx[i][j];
            ddphi_dx[i][j];
        }
        delete[] ddphiC_dx[i];
        delete[] ddphi_dx[i];
    }
    delete[] ddphiC_dx;
    delete[] ddphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainvC_[i];
    delete[] ainvC_;
    
    for (int i = 0; i < 2; ++i) delete[] quadJacMatC[i];
    delete[] quadJacMatC;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersDifferentMesh_ISO_FEM(std::vector<Nodes *> &nodesCoarse_, int *connecC, int &ielem,
                                                             double **lagrMultMatrix, double *rhsVector1, double *rhsVector2)
{

    int dim = 2;
    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;

    double xsiC[dim], phiC_[6];

    double **dphiC_dx;
    dphiC_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[6];

    double **ainvC_;
    ainvC_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **JacC;
    JacC = new double *[dim];
    for (int i = 0; i < dim; ++i) JacC[i] = new double[dim];

    // Data for IGA fine mesh computation (Lagrange field)
    double wpc[9], xsi[dim], phi_[9];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    double **dphi_dx;
    dphi_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **quadJacMat;
    quadJacMat = new double *[dim];
    for (int i = 0; i < dim; ++i) quadJacMat[i] = new double[dim];

    double djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {

        if ((intPointCorrespElem_ISO[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k++) xsi[k] = sQuad.PointListIso(index, k);

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
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_ISO[index][k];

            // Computes the velocity shape function
            shapeQuad.evaluateFem(xsiC, phiC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC, JacC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC, ainvC_, dphiC_dx);

            // Computes Matrix and vectors
            int na = 9;
            int naC = 6;
            getMatrixAndVectorsDifferentMesh(djac_,weight_,na,phi_,dphi_dx,naC,phiC_,dphiC_dx,
                                             lagrMultMatrix,rhsVector1,rhsVector2);
        };
        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;
    
    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;
   
    for (int i = 0; i < 2; ++i) delete[] quadJacMat[i];
    delete[] quadJacMat;

    for (int i = 0; i < 2; ++i) delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < 2; ++i) delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < 2; ++i) delete[] JacC[i];
    delete[] JacC;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_ISO_FEM(std::vector<Nodes *> &nodesCoarse_, int *connecC, int &ielem,
                                                                         double **jacobianNRMatrix, double *rhsVector)
{

    int dim = 2;
    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;

    double xsiC[dim], phiC_[6];

    double **dphiC_dx;
    dphiC_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[6];

    double **ainvC_;
    ainvC_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **JacC;
    JacC = new double *[dim];
    for (int i = 0; i < dim; ++i) JacC[i] = new double[dim];

    // Data for IGA fine mesh computation (Lagrange field)
    double wpc[9], xsi[dim], phi_[9];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    double tSUPG_, tPSPG_, tLSIC_, djac_, weight_;

    double **dphi_dx;
    dphi_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **quadJacMat;
    quadJacMat = new double *[dim];
    for (int i = 0; i < dim; ++i) quadJacMat[i] = new double[dim];


    int index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {

        if ((intPointCorrespElem_ISO[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k++) xsi[k] = sQuad.PointListIso(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListIso(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);

            // Computes the jacobian matrix
            getJacobianMatrix_ISO(djac_,xsi,quadJacMat,ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

            // Computes stabilization parameters

            getNewParameterSUPG_ISO(tSUPG_,tPSPG_,tLSIC_,quadJacMat,phi_,dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++)xsiC[k] = intPointCorrespXsi_ISO[index][k];

            // Computes the velocity shape function
            shapeQuad.evaluateFem(xsiC,phiC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC,JacC,ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC,ainvC_,dphiC_dx);

            // Computes Matrix and vectors
            int na =9;
            int naC = 6;
            getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG(djac_, weight_, tSUPG_, tPSPG_, na, phi_,
                                                         naC, phiC_,dphiC_dx, jacobianNRMatrix, rhsVector);
        };
        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] quadJacMat[i];
    delete[] quadJacMat;

    for (int i = 0; i < 2; ++i) delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < 2; ++i) delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < 2; ++i) delete[] JacC[i];
    delete[] JacC;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersDifferentMeshArlqStab_ISO_FEM(std::vector<Nodes *> &nodesCoarse_, int *connecC, int &ielem,
                                                                     double **arlequinStabD, double *arlequinStabVectorD,
                                                                     double **arlequinStab0, double *arlequinStabVector0)
{
    int dim = 2;
    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;

    double xsiC[dim], phiC_[6];

    double **dphiC_dx;
    dphiC_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[6];

    double **ainvC_;
    ainvC_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **JacC;
    JacC = new double *[dim];
    for (int i = 0; i < dim; ++i) JacC[i] = new double[dim];

    double ***ddphiC_dx;
    ddphiC_dx = new double **[dim];
    for (int i = 0; i < dim; ++i)
    {
        ddphiC_dx[i] = new double *[2];
        for (int j = 0; j < dim; j++)
            ddphiC_dx[i][j] = new double[6];
    };

    // Data for IGA fine mesh computation (Lagrange field)
    double wpc[9], xsi[dim], phi_[9];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    double **dphi_dx;
    dphi_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **quadJacMat;
    quadJacMat = new double *[dim];
    for (int i = 0; i < dim; ++i) quadJacMat[i] = new double[dim];

    double tARLQ_, djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {

        if ((intPointCorrespElem_ISO[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k++) xsi[k] = sQuad.PointListIso(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListIso(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);

            // Computes the jacobian matrix
            getJacobianMatrix_ISO(djac_,xsi,quadJacMat,ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

            // get Arlequin stabilization parameter
            getParameterArlequin_ISO(tARLQ_,phi_,dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++)xsiC[k] = intPointCorrespXsi_ISO[index][k];

            // Computes the velocity shape function
            shapeQuad.evaluateFem(xsiC,phiC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC,JacC,ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC,ainvC_,dphiC_dx);

            // Computes second spatial derivatives
            getSecondSpatialDerivatives_FEM(ainvC_,ddphiC_dx);

            // Computes Matrix and vectors
            int na = 9;
            int naC = 6;
            getMatrixAndVectorsDifferentMeshArlqStab(djac_,weight_,tARLQ_,index,na,dphi_dx,naC,phiC_,dphiC_dx,ddphiC_dx,
                                                     arlequinStabD,arlequinStabVectorD,arlequinStab0,arlequinStabVector0);
        };
        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] quadJacMat[i];
    delete[] quadJacMat;

    for (int i = 0; i < 2; ++i) delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < 2; ++i) delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < 2; ++i) delete[] JacC[i];
    delete[] JacC;

    return;
};

template <>
void Element<2>::getLagrangeMultipliersDifferentMesh_ISO_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_, int *connecC,
                                                             std::vector<IParameters *> &iparamC, int &ielem,
                                                             double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                             double **arlequinStab, double *arlequinStabVector)
{
    int dim = 2;
    // quadrature and functions local classes
    SpecialQuad sQuad = SpecialQuad();
    QuadShapeFunction<2> shapeQuad;

    // Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9], xsiC[dim], phiC_[9];
    for (int i = 0; i < 9; i++) wpcC[i] = (*nodesC_)[connectC_[i]]->getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]]->getINC();

    double **dphiC_dx;
    dphiC_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[9];

    double **ainvC_;
    ainvC_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **quadJacMatC;
    quadJacMatC = new double *[dim];
    for (int i = 0; i < dim; ++i) quadJacMatC[i] = new double[dim];

    // Data for IGA fine mesh computation (Lagrange field)
    double wpc[9], xsi[dim], phi_[9];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    double **dphi_dx;
    dphi_dx = new double *[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **quadJacMat;
    quadJacMat = new double *[dim];
    for (int i = 0; i < dim; ++i) quadJacMat[i] = new double[dim];

    double tARLQ_ , djac_, weight_;

    int index = 0;
    for (double *it = sQuad.beginIso(); it != sQuad.endIso(); it++)
    {

        if ((intPointCorrespElem_ISO[index] == ielem))
        {

            // Fine mesh computation
            // Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k++) xsi[k] = sQuad.PointListIso(index, k);

            // Returns the quadrature integration weight
            weight_ = sQuad.WeightListIso(index);

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, (*iparameters), Npatch_);

            // Computes the jacobian matrix
            getJacobianMatrix_ISO(djac_, xsi, quadJacMat, ainv_);

            // Computes spatial derivatives
            getSpatialDerivatives_ISO(xsi, ainv_, dphi_dx);

            // Arlequin Stabilization parameter
            getParameterArlequin_ISO(tARLQ_, phi_, dphi_dx);

            // Coarse mesh computatation
            // Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_ISO[index][k];

            // Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC, phiC_, wpcC, incC_, (*iparametersC), NpatchC_);

            // Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC, quadJacMatC, ainvC_);

            // Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC, ainvC_, dphiC_dx);

            // Computes Matrix and vectors
            int na = 9;
            getMatrixAndVectorsDifferentMeshTemp(djac_, weight_, tARLQ_, na, phi_, dphi_dx, na, phiC_, dphiC_dx,
                                                 lagrMultMatrix, rhsVector1, rhsVector2,
                                                 arlequinStab, arlequinStabVector);
        };
        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] quadJacMat[i];
    delete[] quadJacMat;

    for (int i = 0; i < 2; ++i) delete[] dphiC_dx[i];
    delete[] dphiC_dx;

    for (int i = 0; i < 2; ++i) delete[] ainvC_[i];
    delete[] ainvC_;

    for (int i = 0; i < 2; ++i) delete[] quadJacMatC[i];
    delete[] quadJacMatC;

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
template <>
void Element<2>::setBoundaryConditions_ISO(double **jacobianNRMatrix, double *rhsVector)
{

    for (int i = 0; i < 9; i++)
    {
        // direction x
        int constrain = (*nodes_)[connect_[i]]->getConstrains(0);
        if ((constrain == 1) || (constrain == 3))
        {
            for (int j = 0; j < 27; j++)
            {
                jacobianNRMatrix[2 * i][j] = 0.0;
                jacobianNRMatrix[j][2 * i] = 0.0;
            }
            jacobianNRMatrix[2 * i][2 * i] = 1.0;
            rhsVector[2 * i] = 0.0;
        }

        // direction y
        constrain = (*nodes_)[connect_[i]]->getConstrains(1);
        if ((constrain == 1) || (constrain == 3))
        {
            for (int j = 0; j < 27; j++)
            {
                jacobianNRMatrix[2 * i + 1][j] = 0.0;
                jacobianNRMatrix[j][2 * i + 1] = 0.0;
            }
            jacobianNRMatrix[2 * i + 1][2 * i + 1] = 1.0;
            rhsVector[2 * i + 1] = 0.0;
        }

        // if PRESSURE CONDITION ADD NODE!
    }
}
//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------

template <>
void Element<2>::getTransientNavierStokes_FEM(double **jacobianNRMatrix, double *rhsVector)
{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // variables
    double phi_[6], xsi[2];

    double tSUPG_, tPSPG_, tLSIC_;
    double djac_, weight_;

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double *[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

    int index = 0;
    for (double *it = nQuad.beginFem(); it != nQuad.endFem(); it++)
    {

        // Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListFem(index,0);
        xsi[1] = nQuad.PointListFem(index,1);
        // Returns the quadrature integration weight
        weight_ = nQuad.WeightListFem(index);

        // Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi,phi_);

        // Computes the jacobian matrixgetElem
        getJacobianMatrix_FEM(djac_,xsi,Jac,ainv_);

        // Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

        // Compute Stabilization Parameters
        getNewParameterSUPG_FEM(tSUPG_,tPSPG_,tLSIC_,Jac,phi_,dphi_dx);

        // Computes the element matrix
        int na = 6;
        getElemMatrix(na,djac_,weight_,tSUPG_,tPSPG_,tLSIC_,index,phi_,dphi_dx,jacobianNRMatrix);

        // Computes the RHS vector
        getResidualVector(na,djac_,weight_,tSUPG_,tPSPG_,tLSIC_,index,phi_,dphi_dx,rhsVector);

        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;
    
    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i) delete[] Jac[i];
    delete[] Jac;

};

template <>
void Element<2>::getTransientNavierStokes_ISO(double **jacobianNRMatrix, double *rhsVector)
{

    // quadrature and functions local classes
    NormalQuad nQuad = NormalQuad();
    QuadShapeFunction<2> shapeQuad;

    // data for IGA elements
    double wpc[9], xsi[2], phi_[9];
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]]->getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]]->getINC();

    double **dphi_dx;
    dphi_dx = new double *[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double *[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double *[2];
    for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];

    double tSUPG_, tPSPG_, tLSIC_;
    double djac_, weight_;

    int index = 0;
    for (double *it = nQuad.beginIso(); it != nQuad.endIso(); it++)
    {

        // Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListIso(index, 0);
        xsi[1] = nQuad.PointListIso(index, 1);
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
        int na = 9;
        getElemMatrix(na, djac_, weight_, tSUPG_, tPSPG_, tLSIC_, index, phi_, dphi_dx, jacobianNRMatrix);

        // Computes the RHS vector
        getResidualVector(na, djac_, weight_, tSUPG_, tPSPG_, tLSIC_, index, phi_, dphi_dx, rhsVector);

        index++;
    };

    for (int i = 0; i < 2; ++i) delete[] dphi_dx[i];
    delete[] dphi_dx;

    for (int i = 0; i < 2; ++i) delete[] ainv_[i];
    delete[] ainv_;

    for (int i = 0; i < 2; ++i)delete[] quadJacMat[i];
    delete[] quadJacMat;
};
