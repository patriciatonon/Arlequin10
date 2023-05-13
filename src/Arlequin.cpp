
#include "Arlequin.h"
//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------SET ELEMENT BOXES--------------------------------
//------------------------------------------------------------------------------
template <int DIM>
void Arlequin<DIM>::setElementBoxes_ISO()
{

    // Compute element boxes for coarse model (IGA coarse mesh)
    double &alpha_f = parametersCoarse->getAlphaF();

    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        int *connec = elementsCoarse_[jel]->getConnectivity();
        double xk[DIM], Xk[DIM];

        int LNN = 18*DIM-27;

        // Bezier transformation matrix
        double **MatrixC_;
        MatrixC_ = new double *[LNN];
        for (int i = 0; i < LNN; i++)
            MatrixC_[i] = new double[LNN];

        elementsCoarse_[jel]->getMatrixC(MatrixC_);

        // Transposing MatrixC
        double transMatrixC_[LNN][LNN];
        for (int i = 0; i < LNN; i++)
            for (int j = 0; j < LNN; j++)
                transMatrixC_[i][j] = MatrixC_[j][i];

        // control points coordinates
        double coord_[LNN][DIM];
        for (int i = 0; i < LNN; i++)
        {
            double *xx = nodesCoarse_[connec[i]]->getCoordinates();
            for (int j = 0; j < DIM; j++)
                coord_[i][j] = xx[j];

            // double *xxp = nodesCoarse_[connec[i]]->getPreviousCoordinates();
            // for (int j = 0; j < DIM; j++) coord_[i][j] = alpha_f * xx[j] + (1. - alpha_f) * xxp[j];
        };

        // Bezier control points coordinates
        double Bcoord_[LNN][DIM] = {};
        for (int i = 0; i < LNN; i++)
            for (int j = 0; j < LNN; j++)
                for (int k = 0; k < DIM; k++)
                    Bcoord_[i][k] += transMatrixC_[i][j] * coord_[j][k];

        for (int i = 0; i < DIM; i++)
        {
            xk[i] = Bcoord_[0][i];
            Xk[i] = Bcoord_[LNN - 1][i];
        }; 

        elementsCoarse_[jel]->setIntersectionParameters(xk, Xk);
    };

    return;
};

template <>
void Arlequin<2>::setElementBoxes_FEM()
{
   
    // Compute element boxes for coarse model (IGA coarse mesh)
    int DIM = 2;
    double &alpha_f = parametersCoarse->getAlphaF();

    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        int *connec = elementsCoarse_[jel] -> getConnectivity();
        double xk[DIM], Xk[DIM],x1[DIM],x2[DIM],x3[DIM];

        double *xx1 = nodesCoarse_[connec[0]] -> getCoordinates();
        double *xx2 = nodesCoarse_[connec[1]] -> getCoordinates();
        double *xx3 = nodesCoarse_[connec[2]] -> getCoordinates();

        for (int i = 0; i < DIM; i++){
            x1[i] = xx1[i];
            x2[i] = xx2[i];
            x3[i] = xx3[i];
        };
        
        for (int i = 0; i < DIM; i++){
            xk[i] = std::min(x1[i],std::min(x2[i],x3[i]));
            Xk[i] = std::max(x1[i],std::max(x2[i],x3[i]));
        };
     
        elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk);

    };

    return;
};

template <>
void Arlequin<3>::setElementBoxes_FEM()
{
    int DIM = 3;
    // Compute element boxes for coarse model (IGA coarse mesh)
    double &alpha_f = parametersCoarse->getAlphaF();

    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        int *connec = elementsCoarse_[jel] -> getConnectivity();
        double xk[DIM], Xk[DIM],x1[DIM],x2[DIM],x3[DIM],x4[DIM];

        double *xx1 = nodesCoarse_[connec[0]] -> getCoordinates();
        double *xx2 = nodesCoarse_[connec[1]] -> getCoordinates();
        double *xx3 = nodesCoarse_[connec[2]] -> getCoordinates();
        double *xx4 = nodesCoarse_[connec[3]] -> getCoordinates();

        for (int i = 0; i < DIM; i++){
            x1[i] = xx1[i];
            x2[i] = xx2[i];
            x3[i] = xx3[i];
            x4[i] = xx4[i];
        };
        
        for (int i = 0; i < DIM; i++){
            xk[i] = std::min(x1[i],std::min(x2[i],std::min(x3[i],x4[i])));
            Xk[i] = std::max(x1[i],std::max(x2[i],std::max(x3[i],x4[i])));
        };
     
        elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk);

    };

    return;
};

//------------------------------------------------------------------------------
//--------------------------SEARCH CORRESPONDENCE-------------------------------
//------------------------------------------------------------------------------
template <>
void Arlequin<2>::searchPointCorrespondence_ISO(double *x, std::vector<Nodes *> nodes,
                                                std::vector<Element *> elements,
                                                std::vector<IsoParameters *> isopar,
                                                int numElem, double *xsiC, int &elemC, int elSearch)
{

    const int DIM = 2;
    int LNNC = 18*DIM-27;

    QuadShapeFunction<DIM> shapeQuad;

    double x_[DIM], deltaX[DIM], deltaXsi[DIM], xsi[DIM], xsiCC[DIM + 1];

    double **ainv;
    ainv = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv[i] = new double[DIM];

    // double &alpha_f = isopar -> getAlphaF();

    for (int i = 0; i < DIM; i++)
    {
        xsi[i] = 0.0; // central element cooordinates
        x_[i] = 0.0;
        xsiC[i] = 1.e50;
    };

    for (int i = 0; i < DIM + 1; i++)
        xsiCC[i] = 1.e10;

    int *connec = elements[elSearch]->getConnectivity();

    // Computing basis functions
    double phi_[LNNC], wpc[LNNC];
    for (int i = 0; i < LNNC; i++)
        wpc[i] = nodes[connec[i]]->getWeightPC();
    int *INC_ = nodes[connec[LNNC - 1]]->getINC();
    int patch = elements[elSearch]->getPatch();
    shapeQuad.evaluateIso(xsi, phi_, wpc, INC_, isopar, patch);

    for (int i = 0; i < LNNC; i++)
    {
        double *xint = nodes[connec[i]]->getCoordinates();
        for (int j = 0; j < DIM; j++)
            x_[j] += xint[j] * phi_[i];
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4))
    {

        iterations++;

        for (int i = 0; i < DIM; i++)
        {
            deltaX[i] = x[i] - x_[i];
            deltaXsi[i] = 0.0;
        };

        elements[elSearch]->getQuadJacobianMatrix_ISO(xsi, ainv);

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                deltaXsi[i] += ainv[i][j] * deltaX[j];
            };
        };

        error = 0.;
        for (int i = 0; i < DIM; i++)
        {
            xsi[i] += deltaXsi[i];
            x_[i] = 0.0;
            error += deltaXsi[i] * deltaXsi[i];
        };

        error = sqrt(error);

        shapeQuad.evaluateIso(xsi, phi_, wpc, INC_, isopar, patch);

        for (int i = 0; i < LNNC; i++)
        {
            double *xint = nodes[connec[i]]->getCoordinates();
            for (int j = 0; j < DIM; j++)
                x_[j] += xint[j] * phi_[i];
        };
    };

    double t1 = -1 - 1.e-2;
    double t2 = 1. + 1.e-2;

    if ((xsi[0] >= t1) && (xsi[1] >= t1) &&
        (xsi[0] <= t2) && (xsi[1] <= t2))
    {

        xsiC[0] = xsi[0];
        xsiC[1] = xsi[1];
        elemC = elSearch;
    }
    else
    {

        for (int jel = 0; jel < numElem; jel++)
        {

            int *connec = elements[jel]->getConnectivity();

            // get boxes information
            std::pair<double *, double *> XK;
            XK = elements[jel]->getXIntersectionParameter();

            // Chech if the node is inside the element box
            if ((x[0] < XK.first[0]) || (x[0] > XK.second[0]) ||
                (x[1] < XK.first[1]) || (x[1] > XK.second[1]))
                continue;

            // Compute the basis functions
            double phi_[LNNC], wpc[LNNC];
            for (int i = 0; i < LNNC; i++)
                wpc[i] = nodes[connec[i]]->getWeightPC();
            int *INC_ = nodes[connec[LNNC - 1]]->getINC();
            int patch = elements[jel]->getPatch();

            for (int i = 0; i < DIM; i++)
            {
                xsi[i] = 0.0; // central element cooordinates
                x_[i] = 0.0;
            };

            shapeQuad.evaluateIso(xsi, phi_, wpc, INC_, isopar, patch);

            for (int i = 0; i < LNNC; i++)
            {
                double *xint = nodes[connec[i]]->getCoordinates();
                for (int j = 0; j < DIM; j++)
                    x_[j] += xint[j] * phi_[i];
            };

            double error = 1.e6;

            int iterations = 0;

            while ((error > 1.e-8) && (iterations < 4))
            {

                iterations++;

                for (int i = 0; i < DIM; i++)
                {
                    deltaX[i] = x[i] - x_[i];
                    deltaXsi[i] = 0.0;
                };

                elements[jel]->getQuadJacobianMatrix_ISO(xsi, ainv);

                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        deltaXsi[i] += ainv[i][j] * deltaX[j];
                    };
                };

                error = 0.;
                for (int i = 0; i < DIM; i++)
                {
                    xsi[i] += deltaXsi[i];
                    x_[i] = 0.0;
                    error = deltaXsi[i] * deltaXsi[i];
                };

                error = sqrt(error);

                shapeQuad.evaluateIso(xsi, phi_, wpc, INC_, isopar, patch);

                for (int i = 0; i < LNNC; i++)
                {
                    double *xint = nodes[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                        x_[j] += xint[j] * phi_[i];
                };
            };

            double t1 = -1 - 1.e-2;
            double t2 = 1. + 1.e-2;

            if ((xsi[0] >= t1) && (xsi[1] >= t1) &&
                (xsi[0] <= t2) && (xsi[1] <= t2))
            {
                xsiC[0] = xsi[0];
                xsiC[1] = xsi[1];
                elemC = jel;
            };
        };
    };

    if (fabs(xsi[0]) > 2.)
        std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE "
                  << std::endl;

    for (int i = 0; i < DIM; ++i)
        delete[] ainv[i];
    delete[] ainv;
};

template <>
void Arlequin<3>::searchPointCorrespondence_ISO(double *x, std::vector<Nodes *> nodes,
                                                std::vector<Element *> elements,
                                                std::vector<IsoParameters *> isopar,
                                                int numElem, double *xsiC, int &elemC, int elSearch)
{

    const int DIM = 3;
    int LNNC = 18*DIM-27;

    QuadShapeFunction<DIM> shapeQuad;

    double x_[DIM], deltaX[DIM], deltaXsi[DIM], xsi[DIM], xsiCC[DIM + 1];

    double **ainv;
    ainv = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv[i] = new double[DIM];

    // double &alpha_f = isopar -> getAlphaF();

    for (int i = 0; i < DIM; i++)
    {
        xsi[i] = 0.0; // central element cooordinates
        x_[i] = 0.0;
        xsiC[i] = 1.e50;
    };

    for (int i = 0; i < DIM + 1; i++)
        xsiCC[i] = 1.e10;

    int *connec = elements[elSearch]->getConnectivity();

    // Computing basis functions
    double phi_[LNNC], wpc[LNNC];
    for (int i = 0; i < LNNC; i++) wpc[i] = nodes[connec[i]]->getWeightPC();
    int *INC_ = nodes[connec[LNNC - 1]]->getINC();
    int patch = elements[elSearch]->getPatch();
    shapeQuad.evaluateIso(xsi, phi_, wpc, INC_, isopar, patch);

    for (int i = 0; i < LNNC; i++)
    {
        double *xint = nodes[connec[i]]->getCoordinates();
        for (int j = 0; j < DIM; j++)
            x_[j] += xint[j] * phi_[i];
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4))
    {

        iterations++;

        for (int i = 0; i < DIM; i++)
        {
            deltaX[i] = x[i] - x_[i];
            deltaXsi[i] = 0.0;
        };

        elements[elSearch]->getQuadJacobianMatrix_ISO(xsi, ainv);

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                deltaXsi[i] += ainv[i][j] * deltaX[j];

        error = 0.;
        for (int i = 0; i < DIM; i++)
        {
            xsi[i] += deltaXsi[i];
            x_[i] = 0.0;
            error += deltaXsi[i] * deltaXsi[i];
        };

        error = sqrt(error);

        shapeQuad.evaluateIso(xsi, phi_, wpc, INC_, isopar, patch);

        for (int i = 0; i < LNNC; i++)
        {
            double *xint = nodes[connec[i]]->getCoordinates();
            for (int j = 0; j < DIM; j++)
                x_[j] += xint[j] * phi_[i];
        };
    };

    double t1 = -1 - 1.e-2;
    double t2 = 1. + 1.e-2;

    if ((xsi[0] >= t1) && (xsi[1] >= t1) && (xsi[2] >= t1) &&
        (xsi[0] <= t2) && (xsi[1] <= t2) && (xsi[2] <= t2))
    {
        xsiC[0] = xsi[0];
        xsiC[1] = xsi[1];
        xsiC[2] = xsi[2];
        elemC = elSearch;
    }
   
    else
    {

        for (int jel = 0; jel < numElem; jel++)
        {

            int *connec = elements[jel]->getConnectivity();

            // get boxes information
            std::pair<double *, double *> XK;
            XK = elements[jel]->getXIntersectionParameter();


            //Chech if the node is inside the element box
            if ((x[0] < XK.first[0] - 0.0001) || (x[0] > XK.second[0] + 0.001) ||
                (x[1] < XK.first[1] - 0.0001) || (x[1] > XK.second[1] + 0.001) ||
                (x[2] < XK.first[2] - 0.0001) || (x[2] > XK.second[2] + 0.001))
                continue;

            // Compute the basis functions
            double phi_[LNNC], wpc[LNNC];
            for (int i = 0; i < LNNC; i++)
                wpc[i] = nodes[connec[i]]->getWeightPC();
            int *INC_ = nodes[connec[LNNC - 1]]->getINC();
            int patch = elements[jel]->getPatch();

            for (int i = 0; i < DIM; i++)
            {
                xsi[i] = 0.0; // central element cooordinates
                x_[i] = 0.0;
            };

            shapeQuad.evaluateIso(xsi, phi_, wpc, INC_, isopar, patch);

            for (int i = 0; i < LNNC; i++)
            {
                double *xint = nodes[connec[i]]->getCoordinates();
                for (int j = 0; j < DIM; j++)
                    x_[j] += xint[j] * phi_[i];
            };

            double error = 1.e6;

            int iterations = 0;

            while ((error > 1.e-8) && (iterations < 4))
            {

                iterations++;

                for (int i = 0; i < DIM; i++)
                {
                    deltaX[i] = x[i] - x_[i];
                    deltaXsi[i] = 0.0;
                };

                elements[jel]->getQuadJacobianMatrix_ISO(xsi, ainv);

                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        deltaXsi[i] += ainv[i][j] * deltaX[j];

                error = 0.;
                for (int i = 0; i < DIM; i++)
                {
                    xsi[i] += deltaXsi[i];
                    x_[i] = 0.0;
                    error = deltaXsi[i] * deltaXsi[i];
                };

                error = sqrt(error);

                shapeQuad.evaluateIso(xsi, phi_, wpc, INC_, isopar, patch);

                for (int i = 0; i < LNNC; i++)
                {
                    double *xint = nodes[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                        x_[j] += xint[j] * phi_[i];
                };
            };

            double t1 = -1 - 1.e-2;
            double t2 = 1. + 1.e-2;

            if ((xsi[0] >= t1) && (xsi[1] >= t1) && (xsi[2] >= t1) &&
                (xsi[0] <= t2) && (xsi[1] <= t2) && (xsi[2] <= t2))
            {
                xsiC[0] = xsi[0];
                xsiC[1] = xsi[1];
                xsiC[2] = xsi[2];
                elemC = jel;
                break;
            };
        };
    };

    if (fabs(xsi[0]) > 2.){
        std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE "
                  << std::endl;
    };
        

    for (int i = 0; i < DIM; ++i)
        delete[] ainv[i];
    delete[] ainv;
};


template <>
void Arlequin<2>::searchPointCorrespondence_FEM(int ielem, int ip, double *x, std::vector<Nodes *> nodes,
                                                  std::vector<Element *> elements,
                                                  int numElem, double *xsiC, int &elemC, int elSearch)
{

    int const DIM = 2;

    int LNNC = 4*DIM-2;

    QuadShapeFunction<DIM> shapeQuad;

    double x_[DIM], deltaX[DIM], deltaXsi[DIM], xsi[DIM], xsiCC[DIM + 1];

    double **ainv;
    ainv = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv[i] = new double[DIM];

    // double &alpha_f = isopar -> getAlphaF();

    for (int i = 0; i < DIM; i++)
    {
        xsi[i] = 1./3.; // central element cooordinates
        x_[i] = 0.0;
        xsiC[i] = 1.e50;
    };

    for (int i = 0; i < DIM + 1; i++) xsiCC[i] = 1.e10;

    int *connec = elements[elSearch]->getConnectivity();

    // Computing basis functions
    double phi_[LNNC];
    shapeQuad.evaluateFem(xsi, phi_);

    for (int i = 0; i < LNNC; i++)
    {
        double *xint = nodes[connec[i]]->getCoordinates();
        for (int j = 0; j < DIM; j++) x_[j] += xint[j] * phi_[i];
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4))
    {

        iterations++;

        for (int i = 0; i < DIM; i++)
        {
            deltaX[i] = x[i] - x_[i];
            deltaXsi[i] = 0.0;
        };

        elements[elSearch]->getJacobianMatrixValues_FEM(xsi, ainv);

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                deltaXsi[i] += ainv[i][j] * deltaX[j];

        error = 0.;
        for (int i = 0; i < DIM; i++)
        {
            xsi[i] += deltaXsi[i];
            x_[i] = 0.0;
            error += deltaXsi[i] * deltaXsi[i];
        };

        error = sqrt(error);

        shapeQuad.evaluateFem(xsi, phi_);

        for (int i = 0; i < LNNC; i++)
        {
            double *xint = nodes[connec[i]]->getCoordinates();
            for (int j = 0; j < DIM; j++)
                x_[j] += xint[j] * phi_[i];
        };
    };

    double t1 = -1.e-2;
    double t2 = 1. - t1;

    if ((xsi[0] >= t1) && (xsi[1] >= t1) &&  
        (xsi[0] <= t2) && (xsi[1] <= t2))
    {
        xsiC[0] = xsi[0];
        xsiC[1] = xsi[1];
        elemC = elSearch;
    }
   
    else
    {

        for (int jel = 0; jel < numElem; jel++)
        {

            int *connec = elements[jel]->getConnectivity();

            // get boxes information
            std::pair<double *, double *> XK;
            XK = elements[jel]->getXIntersectionParameter();

            //Chech if the node is inside the element box
            if ((x[0] < XK.first[0] - 0.0001) || (x[0] > XK.second[0] + 0.001) ||
                (x[1] < XK.first[1] - 0.0001) || (x[1] > XK.second[1] + 0.001))
                continue;

            // Compute the basis functions
            double phi_[LNNC];
            for (int i = 0; i < DIM; i++)
            {
                xsi[i] = 1./3.; // central element cooordinates
                x_[i] = 0.0;
            };

            shapeQuad.evaluateFem(xsi, phi_);

            for (int i = 0; i < LNNC; i++)
            {
                double *xint = nodes[connec[i]]->getCoordinates();
                for (int j = 0; j < DIM; j++)
                    x_[j] += xint[j] * phi_[i];
            };

            double error = 1.e6;

            int iterations = 0;

            while ((error > 1.e-8) && (iterations < 4))
            {

                iterations++;

                for (int i = 0; i < DIM; i++)
                {
                    deltaX[i] = x[i] - x_[i];
                    deltaXsi[i] = 0.0;
                };

                elements[jel]->getJacobianMatrixValues_FEM(xsi, ainv);

                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        deltaXsi[i] += ainv[i][j] * deltaX[j];

                error = 0.;
                for (int i = 0; i < DIM; i++)
                {
                    xsi[i] += deltaXsi[i];
                    x_[i] = 0.0;
                    error = deltaXsi[i] * deltaXsi[i];
                };

                error = sqrt(error);

                shapeQuad.evaluateFem(xsi, phi_);

                for (int i = 0; i < LNNC; i++)
                {
                    double *xint = nodes[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                        x_[j] += xint[j] * phi_[i];
                };
            };

               double t1 = -1.e-2;
               double t2 = 1. - t1;

            if ((xsi[0] >= t1) && (xsi[1] >= t1) &&
                (xsi[0] <= t2) && (xsi[1] <= t2))
            {
                xsiC[0] = xsi[0];
                xsiC[1] = xsi[1];
                elemC = jel;
                break;
            };
        };
    };

    if (fabs(xsi[0]) > 2.){

        std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE "
                  << std::endl;
    };
        

    for (int i = 0; i < DIM; ++i)
        delete[] ainv[i];
    delete[] ainv;
};

template <>
void Arlequin<3>::searchPointCorrespondence_FEM(int ielem, int ip, double *x, std::vector<Nodes *> nodes,
                                                  std::vector<Element *> elements,
                                                  int numElem, double *xsiC, int &elemC, int elSearch)
{
    int const DIM = 3;
    int LNNC = 4*DIM-2;

    QuadShapeFunction<DIM> shapeQuad;

    double x_[DIM], deltaX[DIM], deltaXsi[DIM], xsi[DIM], xsiCC[DIM + 1];

    double **ainv;
    ainv = new double *[DIM];
    for (int i = 0; i < DIM; ++i)
        ainv[i] = new double[DIM];

    // double &alpha_f = isopar -> getAlphaF();

    for (int i = 0; i < DIM; i++)
    {
        xsi[i] = 1./3.; // central element cooordinates
        x_[i] = 0.0;
        xsiC[i] = 1.e50;
    };

    for (int i = 0; i < DIM + 1; i++) xsiCC[i] = 1.e10;

    int *connec = elements[elSearch]->getConnectivity();

    // Computing basis functions
    double phi_[LNNC];
    shapeQuad.evaluateFem(xsi, phi_);

    for (int i = 0; i < LNNC; i++)
    {
        double *xint = nodes[connec[i]]->getCoordinates();
        for (int j = 0; j < DIM; j++) x_[j] += xint[j] * phi_[i];
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4))
    {

        iterations++;

        for (int i = 0; i < DIM; i++)
        {
            deltaX[i] = x[i] - x_[i];
            deltaXsi[i] = 0.0;
        };

        elements[elSearch]->getJacobianMatrixValues_FEM(xsi, ainv);

        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                deltaXsi[i] += ainv[i][j] * deltaX[j];

        error = 0.;
        for (int i = 0; i < DIM; i++)
        {
            xsi[i] += deltaXsi[i];
            x_[i] = 0.0;
            error += deltaXsi[i] * deltaXsi[i];
        };

        error = sqrt(error);

        shapeQuad.evaluateFem(xsi, phi_);

        for (int i = 0; i < LNNC; i++)
        {
            double *xint = nodes[connec[i]]->getCoordinates();
            for (int j = 0; j < DIM; j++)
                x_[j] += xint[j] * phi_[i];
        };
    };

    double t1 = -1.e-2;
    double t2 = 1. - t1;

    if ((xsi[0] >= t1) && (xsi[1] >= t1) && (xsi[2] >= t1) &&
        (xsi[0] <= t2) && (xsi[1] <= t2) && (xsi[2] <= t2))
    {
        xsiC[0] = xsi[0];
        xsiC[1] = xsi[1];
        xsiC[2] = xsi[2];
        elemC = elSearch;
    }
   
    else
    {


        for (int jel = 0; jel < numElem; jel++)
        {

            int *connec = elements[jel]->getConnectivity();

            // get boxes information
            std::pair<double *, double *> XK;
            XK = elements[jel]->getXIntersectionParameter();

            //Chech if the node is inside the element box
            if ((x[0] < XK.first[0] - 0.0001) || (x[0] > XK.second[0] + 0.001) ||
                (x[1] < XK.first[1] - 0.0001) || (x[1] > XK.second[1] + 0.001) ||
                (x[2] < XK.first[2] - 0.0001) || (x[2] > XK.second[2] + 0.001))
                continue;

            // Compute the basis functions
            double phi_[LNNC];
            for (int i = 0; i < DIM; i++)
            {
                xsi[i] = 1./3.; // central element cooordinates
                x_[i] = 0.0;
            };

            shapeQuad.evaluateFem(xsi, phi_);

            for (int i = 0; i < LNNC; i++)
            {
                double *xint = nodes[connec[i]]->getCoordinates();
                for (int j = 0; j < DIM; j++)
                    x_[j] += xint[j] * phi_[i];
            };

            double error = 1.e6;

            int iterations = 0;

            while ((error > 1.e-8) && (iterations < 4))
            {

                iterations++;

                for (int i = 0; i < DIM; i++)
                {
                    deltaX[i] = x[i] - x_[i];
                    deltaXsi[i] = 0.0;
                };

                elements[jel]->getJacobianMatrixValues_FEM(xsi, ainv);

                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        deltaXsi[i] += ainv[i][j] * deltaX[j];

                error = 0.;
                for (int i = 0; i < DIM; i++)
                {
                    xsi[i] += deltaXsi[i];
                    x_[i] = 0.0;
                    error = deltaXsi[i] * deltaXsi[i];
                };

                error = sqrt(error);

                shapeQuad.evaluateFem(xsi, phi_);

                for (int i = 0; i < LNNC; i++)
                {
                    double *xint = nodes[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                        x_[j] += xint[j] * phi_[i];
                };
            };

               double t1 = -1.e-2;
               double t2 = 1. - t1;

            if ((xsi[0] >= t1) && (xsi[1] >= t1) && (xsi[2] >= t1) &&
                (xsi[0] <= t2) && (xsi[1] <= t2) && (xsi[2] <= t2))
            {
                xsiC[0] = xsi[0];
                xsiC[1] = xsi[1];
                xsiC[2] = xsi[2];
                elemC = jel;
                break;
            };
        };
    };

    if (fabs(xsi[0]) > 2.){

        std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE "
                  << std::endl;
    };
        

    for (int i = 0; i < DIM; ++i)
        delete[] ainv[i];
    delete[] ainv;
};

template <>
void Arlequin<2>::setCorrespondenceFine_FEM_ISO()
{
    int DIM = 2;

    // Node correspondence
    for (int inode = 0; inode < numNodesGlueZoneFine; inode++)
    {

        double *x = nodesFine_[nodesGlueZoneFine_[inode]]->getCoordinates();

        int elemC = 0;
        double xsiC[DIM] = {};

        searchPointCorrespondence_ISO(x, nodesCoarse_, elementsCoarse_, IsoParCoarse,
                                      elementsCoarse_.size(), xsiC, elemC,
                                      nodesFine_[nodesGlueZoneFine_[inode]]->getNodalElemCorrespondence());

        nodesFine_[nodesGlueZoneFine_[inode]]->setNodalCorrespondence(elemC, xsiC);

    };

    // integration points correspondence
    int LNN = 4*DIM-2;
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        SpecialQuadrature squad;

        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        double x1[LNN], x2[LNN];
        for (int j = 0; j < LNN; j++)
        {
            double *x = nodesFine_[connec[j]]->getCoordinates();
            x1[j] = x[0];
            x2[j] = x[1];
        };

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]]->getNumberOfIntegrationPointsSpecial_FEM();

        for (int ip = 0; ip < numberIntPoints; ip++)
        {

            int elemC = 0;
            double xsiC[DIM] = {};
            double x_[DIM] = {};

            x_[0] = squad.interpolateQuadraticVariableFem(x1, ip);
            x_[1] = squad.interpolateQuadraticVariableFem(x2, ip);

            searchPointCorrespondence_ISO(x_, nodesCoarse_, elementsCoarse_, IsoParCoarse,
                                          elementsCoarse_.size(), xsiC, elemC,
                                          elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_FEM(ip));

            elementsFine_[elementsGlueZoneFine_[i]]->setIntegrationPointCorrespondence_FEM(ip, xsiC, elemC);

        };
    };
};

template <>
void Arlequin<3>::setCorrespondenceFine_FEM_ISO()
{
    int DIM = 3;
    // Node correspondence
    for (int inode = 0; inode < numNodesGlueZoneFine; inode++)
    {

        double *x = nodesFine_[nodesGlueZoneFine_[inode]]->getCoordinates();

        int elemC = 0;
        double xsiC[DIM] = {};

        searchPointCorrespondence_ISO(x, nodesCoarse_, elementsCoarse_, IsoParCoarse,
                                      elementsCoarse_.size(), xsiC, elemC,
                                      nodesFine_[nodesGlueZoneFine_[inode]]->getNodalElemCorrespondence());

        nodesFine_[nodesGlueZoneFine_[inode]]->setNodalCorrespondence(elemC, xsiC);


    };

    // integration points correspondence
    int LNN = 4*DIM-2;
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        SpecialQuadrature squad;

        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        double x1[LNN], x2[LNN], x3[LNN];
        for (int j = 0; j < LNN; j++)
        {
            double *x = nodesFine_[connec[j]]->getCoordinates();
            x1[j] = x[0];
            x2[j] = x[1];
            x3[j] = x[2];
        };

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]]->getNumberOfIntegrationPointsSpecial_FEM();

        for (int ip = 0; ip < numberIntPoints; ip++)
        {

            int elemC = 0;
            double xsiC[DIM] = {};
            double x_[DIM] = {};

            x_[0] = squad.interpolateQuadraticVariableFem(x1, ip);
            x_[1] = squad.interpolateQuadraticVariableFem(x2, ip);
            x_[2] = squad.interpolateQuadraticVariableFem(x3, ip);

            searchPointCorrespondence_ISO(x_, nodesCoarse_, elementsCoarse_, IsoParCoarse,
                                          elementsCoarse_.size(), xsiC, elemC,
                                          elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_FEM(ip));

            elementsFine_[elementsGlueZoneFine_[i]]->setIntegrationPointCorrespondence_FEM(ip, xsiC, elemC);

        };
    };
};

template <>
void Arlequin<2>::setCorrespondenceFine_FEM_FEM()
{
    int DIM = 2;

    // Node correspondence
    for (int inode = 0; inode < numNodesGlueZoneFine; inode++)
    {
       
        double *x = nodesFine_[nodesGlueZoneFine_[inode]]->getCoordinates();

        int elemC = 0;
        double xsiC[DIM] = {};
        int const val = 1;
        searchPointCorrespondence_FEM(nodesGlueZoneFine_[inode], val, x, nodesCoarse_, elementsCoarse_,
                                      elementsCoarse_.size(), xsiC, elemC,
                                      nodesFine_[nodesGlueZoneFine_[inode]]->getNodalElemCorrespondence());

        nodesFine_[nodesGlueZoneFine_[inode]]->setNodalCorrespondence(elemC, xsiC);


    };

    // integration points correspondence
    int LNN = 4*DIM-2;
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        SpecialQuadrature squad;

        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        double x1[LNN], x2[LNN];
        for (int j = 0; j < LNN; j++)
        {
            double *x = nodesFine_[connec[j]]->getCoordinates();
            x1[j] = x[0];
            x2[j] = x[1];
        };

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]]->getNumberOfIntegrationPointsSpecial_FEM();

        for (int ip = 0; ip < numberIntPoints; ip++)
        {

            int elemC = 0;
            double xsiC[DIM] = {};
            double x_[DIM] = {};

            x_[0] = squad.interpolateQuadraticVariableFem(x1, ip);
            x_[1] = squad.interpolateQuadraticVariableFem(x2, ip);

            searchPointCorrespondence_FEM(elementsGlueZoneFine_[i], ip, x_, nodesCoarse_, elementsCoarse_,
                                          elementsCoarse_.size(), xsiC, elemC,
                                          elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_FEM(ip));
            
           
            elementsFine_[elementsGlueZoneFine_[i]]->setIntegrationPointCorrespondence_FEM(ip, xsiC, elemC);
                  
        };
    };
};

template <>
void Arlequin<3>::setCorrespondenceFine_FEM_FEM()
{
    int DIM = 3;
    // Node correspondence
    for (int inode = 0; inode < numNodesGlueZoneFine; inode++)
    {
       
        double *x = nodesFine_[nodesGlueZoneFine_[inode]]->getCoordinates();

        int elemC = 0;
        double xsiC[DIM] = {};
        int const val = 1;
        searchPointCorrespondence_FEM(nodesGlueZoneFine_[inode], val, x, nodesCoarse_, elementsCoarse_,
                                      elementsCoarse_.size(), xsiC, elemC,
                                      nodesFine_[nodesGlueZoneFine_[inode]]->getNodalElemCorrespondence());

        nodesFine_[nodesGlueZoneFine_[inode]]->setNodalCorrespondence(elemC, xsiC);


    };

    // integration points correspondence
    int LNN = 4*DIM-2;
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        SpecialQuadrature squad;

        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        double x1[LNN], x2[LNN], x3[LNN];
        for (int j = 0; j < LNN; j++)
        {
            double *x = nodesFine_[connec[j]]->getCoordinates();
            x1[j] = x[0];
            x2[j] = x[1];
            x3[j] = x[2];
        };

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]]->getNumberOfIntegrationPointsSpecial_FEM();

        for (int ip = 0; ip < numberIntPoints; ip++)
        {

            int elemC = 0;
            double xsiC[DIM] = {};
            double x_[DIM] = {};

            x_[0] = squad.interpolateQuadraticVariableFem(x1, ip);
            x_[1] = squad.interpolateQuadraticVariableFem(x2, ip);
            x_[2] = squad.interpolateQuadraticVariableFem(x3, ip);

            searchPointCorrespondence_FEM(elementsGlueZoneFine_[i], ip, x_, nodesCoarse_, elementsCoarse_,
                                          elementsCoarse_.size(), xsiC, elemC,
                                          elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_FEM(ip));
            
           
            elementsFine_[elementsGlueZoneFine_[i]]->setIntegrationPointCorrespondence_FEM(ip, xsiC, elemC);

        };
    };
};


template <>
void Arlequin<2>::setCorrespondenceFine_ISO_ISO()
{
    int DIM = 2;

    // Node correspondence
    for (int inode = 0; inode < numNodesGlueZoneFine; inode++)
    {
        double *x = nodesFine_[nodesGlueZoneFine_[inode]]->getCoordinates();
        int elemC = 0;
        double xsiC[DIM] = {};

        searchPointCorrespondence_ISO(x, nodesCoarse_, elementsCoarse_, IsoParCoarse,
                                      elementsCoarse_.size(), xsiC, elemC,
                                      nodesFine_[nodesGlueZoneFine_[inode]]->getNodalElemCorrespondence());

        nodesFine_[nodesGlueZoneFine_[inode]]->setNodalCorrespondence(elemC, xsiC);
    };

    // integration points correspondence
    int LNN = 18*DIM-27;
    
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        SpecialQuadrature squad;
        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        double x1[LNN], x2[LNN];
        for (int j = 0; j < LNN; j++)
        {
            double *x = nodesFine_[connec[j]]->getCoordinates();
            x1[j] = x[0];
            x2[j] = x[1];
        };

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]]->getNumberOfIntegrationPointsSpecial_ISO();

        double wpc[LNN];
        for (int i = 0; i < LNN; i++) wpc[i] = nodesFine_[connec[i]] -> getWeightPC();
        int *inc = nodesFine_[connec[LNN-1]] -> getINC();
        int patch = elementsFine_[elementsGlueZoneFine_[i]]->getPatch();

        for (int ip = 0; ip < numberIntPoints; ip++)
        {

            int elemC = 0;
            double xsiC[DIM] = {};
            double x_[DIM] = {};

            x_[0] = squad.interpolateQuadraticVariableIso(x1, ip , wpc, inc, IsoParFine, patch);
            x_[1] = squad.interpolateQuadraticVariableIso(x2, ip, wpc, inc, IsoParFine, patch);

            searchPointCorrespondence_ISO(x_, nodesCoarse_, elementsCoarse_, IsoParCoarse,
                                          elementsCoarse_.size(), xsiC, elemC,
                                          elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_ISO(ip));
           
            elementsFine_[elementsGlueZoneFine_[i]]->setIntegrationPointCorrespondence_ISO(ip, xsiC, elemC);

        };
    };
};

template <>
void Arlequin<3>::setCorrespondenceFine_ISO_ISO()
{
    int DIM = 3;
    // Node correspondence
    for (int inode = 0; inode < numNodesGlueZoneFine; inode++)
    {
        double *x = nodesFine_[nodesGlueZoneFine_[inode]]->getCoordinates();
        int elemC = 0;
        double xsiC[DIM] = {};

        searchPointCorrespondence_ISO(x, nodesCoarse_, elementsCoarse_, IsoParCoarse,
                                      elementsCoarse_.size(), xsiC, elemC,
                                      nodesFine_[nodesGlueZoneFine_[inode]]->getNodalElemCorrespondence());

        nodesFine_[nodesGlueZoneFine_[inode]]->setNodalCorrespondence(elemC, xsiC);
    };

    // integration points correspondence
    int LNN = 18*DIM-27;
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        SpecialQuadrature squad;

        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        double x1[LNN], x2[LNN], x3[LNN];
        for (int j = 0; j < LNN; j++)
        {
            double *x = nodesFine_[connec[j]]->getCoordinates();
            x1[j] = x[0];
            x2[j] = x[1];
            x3[j] = x[2];
        };

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]]->getNumberOfIntegrationPointsSpecial_ISO();

        double wpc[LNN];
        for (int i = 0; i < LNN; i++) wpc[i] = nodesFine_[connec[i]] -> getWeightPC();
        int *inc = nodesFine_[connec[LNN-1]] -> getINC();
        int patch = elementsFine_[elementsGlueZoneFine_[i]]->getPatch();


        for (int ip = 0; ip < numberIntPoints; ip++)
        {

            int elemC = 0;
            double xsiC[DIM] = {};
            double x_[DIM] = {};

            x_[0] = squad.interpolateQuadraticVariableIso(x1, ip , wpc, inc, IsoParFine, patch);
            x_[1] = squad.interpolateQuadraticVariableIso(x2, ip, wpc, inc, IsoParFine, patch);
            x_[2] = squad.interpolateQuadraticVariableIso(x3, ip, wpc, inc, IsoParFine, patch);

            searchPointCorrespondence_ISO(x_, nodesCoarse_, elementsCoarse_, IsoParCoarse,
                                          elementsCoarse_.size(), xsiC, elemC,
                                          elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_ISO(ip));
           
            elementsFine_[elementsGlueZoneFine_[i]]->setIntegrationPointCorrespondence_ISO(ip, xsiC, elemC);

        };
    };
};


//------------------------------------------------------------------------------
//-------------------------FIND SIGNALADED DISTANCE-----------------------------
//------------------------------------------------------------------------------
template <>
void Arlequin<2>::setSignaledDistance_FEM_ISO()
{

    int dim = 2;

    double x[dim], x1[dim], x1B[dim], x2[dim], x2B[dim], x3B[dim],
        n[dim], test[dim], dist;
    int bconnec[3];

    for (int inode = 0; inode < numNodesFine; inode++)
    {
        nodesFine_[inode]->clearInnerNormal();
        nodesFine_[inode]->setDistFunction(0.0);
    };

    for (int inode = 0; inode < numNodesCoarse; inode++)
    {
        nodesCoarse_[inode]->setDistFunction(0.0);
    };

    // Normal vector from a defined boundary in fine mesh
    for (int ibound = 0; ibound < numBoundElemFine; ibound++)
    {

        if (boundaryFine_[ibound]->getConstrain(0) == 2)
        {

            // double &alpha_f = parametersFine -> getAlphaF();
            int *connec = elementsFine_[boundaryFine_[ibound]->getElement()]->getConnectivity();

            // Recognizing the element side
            std::pair<std::vector<int>, std::vector<int>> elemBound;
            elemBound = elementsFine_[boundaryFine_[ibound]->getElement()]->getElemSideInBoundary();
            int numofBoundaries = elemBound.first.size();
            int side;
            for (int j = 0; j < numofBoundaries; j++)
            {
                if (elemBound.second[j] == ibound)
                {
                    side = elemBound.first[j];
                };
            };

            if (side == 0)
            {
                bconnec[0] = connec[1];
                bconnec[1] = connec[2];
                bconnec[2] = connec[4];
            };
            if (side == 1)
            {
                bconnec[0] = connec[2];
                bconnec[1] = connec[0];
                bconnec[2] = connec[5];
            };
            if (side == 2)
            {
                bconnec[0] = connec[0];
                bconnec[1] = connec[1];
                bconnec[2] = connec[3];
            };

            // first segment
            int no1 = bconnec[0];
            int no2 = bconnec[2];

            double *xx1 = nodesFine_[no1]->getCoordinates();
            double *xx2 = nodesFine_[no2]->getCoordinates();
            // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
            // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

            for (int k = 0; k < dim; k++)
            {
                x1[k] = xx1[k];
                x2[k] = xx2[k];
                // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
            }

            double sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                  (x1[0] - x2[0]) * (x1[0] - x2[0]));

            n[0] = ((x2[1] - x1[1]) / sLength);
            n[1] = ((x1[0] - x2[0]) / sLength);

            nodesFine_[no1]->setInnerNormal(n);
            nodesFine_[no2]->setInnerNormal(n);

            // second segment
            no1 = bconnec[2];
            no2 = bconnec[1];

            xx1 = nodesFine_[no1]->getCoordinates();
            xx2 = nodesFine_[no2]->getCoordinates();
            // xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
            // xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

            for (int k = 0; k < dim; k++)
            {
                x1[k] = xx1[k];
                x2[k] = xx2[k];
                // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
            }

            sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                           (x1[0] - x2[0]) * (x1[0] - x2[0]));

            n[0] = ((x2[1] - x1[1]) / sLength);
            n[1] = ((x1[0] - x2[0]) / sLength);

            nodesFine_[no1]->setInnerNormal(n);
            nodesFine_[no2]->setInnerNormal(n);

        }; // constrain == 2
    };     // numBoundFine

    // Gambiarra para malha FEM fine
    n[0] = 1.;
    n[1] = 1.;
    nodesFine_[6]->setInnerNormal(n);
    n[0] = -1.;
    n[1] = 1.;
    nodesFine_[8]->setInnerNormal(n);
    n[0] = -1.;
    n[1] = -1.;
    nodesFine_[17]->setInnerNormal(n);
    n[0] = 1.;
    n[1] = -1.;
    nodesFine_[15]->setInnerNormal(n);

    // Coarse mesh - closer distance for nodes or control points from defined fine boundary
    for (int ino = 0; ino < numNodesCoarse; ino++)
    {

        double *xx = nodesCoarse_[ino]->getCoordinates();
        for (int i = 0; i < dim; i++)
            x[i] = xx[i];

        // double &alphaC_f = parametersCoarse -> getAlphaF();
        // double *xx = nodesCoarse_[ino]->getCoordinates();
        // double *xxp = nodesCoarse_[ino]->getPreviousCoordinates();
        // for (int i = 0; i < dim; i++) x[i] = alphaC_f * xx[i] + (1. - alphaC_f) * xxp[i];

        dist = 10000000000000000000000000000.;

        for (int ibound = 0; ibound < numBoundElemFine; ibound++)
        {

            if (boundaryFine_[ibound]->getConstrain(0) == 2)
            {

                // double &alpha_f = parametersFine -> getAlphaF();
                int *connec = elementsFine_[boundaryFine_[ibound]->getElement()]->getConnectivity();

                // side in the boundary
                std::pair<std::vector<int>, std::vector<int>> elemBound;
                elemBound = elementsFine_[boundaryFine_[ibound]->getElement()]->getElemSideInBoundary();
                int numofBoundaries = elemBound.first.size();
                int side;
                for (int j = 0; j < numofBoundaries; j++)
                {
                    if (elemBound.second[j] == ibound)
                    {
                        side = elemBound.first[j];
                    };
                };

                if (side == 0)
                {
                    bconnec[0] = connec[1];
                    bconnec[1] = connec[2];
                    bconnec[2] = connec[4];
                };
                if (side == 1)
                {
                    bconnec[0] = connec[2];
                    bconnec[1] = connec[0];
                    bconnec[2] = connec[5];
                };
                if (side == 2)
                {
                    bconnec[0] = connec[0];
                    bconnec[1] = connec[1];
                    bconnec[2] = connec[3];
                };

                // first segment
                int no1, no2;
                no1 = bconnec[0];
                no2 = bconnec[2];

                double *xx1 = nodesFine_[no1]->getCoordinates();
                double *xx2 = nodesFine_[no2]->getCoordinates();
                // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

                for (int k = 0; k < dim; k++)
                {
                    x1[k] = xx1[k];
                    x2[k] = xx2[k];
                    // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                    // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
                }

                double aux0 = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                   (x2[0] - x1[0]) * (x2[0] - x1[0]));
                double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0]) +
                               (x[1] - x1[1]) * (x2[1] - x1[1])) /
                              aux0;
                double dist2 = -((x2[1] - x1[1]) * x[0] -
                                 (x2[0] - x1[0]) * x[1] +
                                 x2[0] * x1[1] - x2[1] * x1[0]) /
                               aux0;

                if (aux1 > aux0)
                {

                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    // find signal
                    // side normal vector
                    double *normal = nodesFine_[no2]->getInnerNormal();

                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];

                    double signaltest = 0.0;
                    for (int i = 0; i < 2; i++)
                    {
                        signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)
                        signal = 1.;

                    dist2 *= signal;
                };

                if (aux1 < 0.)
                {

                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    // find signal
                    // side normal vector
                    double *normal = nodesFine_[no1]->getInnerNormal();

                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];

                    double signaltest = 0.0;
                    for (int i = 0; i < 2; i++)
                    {
                        signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)
                        signal = 1.;

                    dist2 *= signal;
                };

                if (fabs(dist2) < fabs(dist))
                    dist = dist2;

                // second segment
                no1 = bconnec[2];
                no2 = bconnec[1];

                xx1 = nodesFine_[no1]->getCoordinates();
                xx2 = nodesFine_[no2]->getCoordinates();
                // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

                for (int k = 0; k < dim; k++)
                {
                    x1[k] = xx1[k];
                    x2[k] = xx2[k];
                    // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                    // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
                }

                aux0 = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                            (x2[0] - x1[0]) * (x2[0] - x1[0]));
                aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0]) +
                        (x[1] - x1[1]) * (x2[1] - x1[1])) /
                       aux0;
                dist2 = -((x2[1] - x1[1]) * x[0] - (x2[0] - x1[0]) * x[1] +
                          x2[0] * x1[1] - x2[1] * x1[0]) /
                        aux0;

                if (aux1 > aux0)
                {
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));

                    double *normal = nodesFine_[no2]->getInnerNormal();

                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];

                    double signaltest = 0.0;
                    for (int i = 0; i < 2; i++)
                    {
                        signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)
                        signal = 1.;

                    dist2 *= signal;
                };

                if (aux1 < 0.)
                {
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));

                    // find signal
                    // side normal vector
                    double *normal;
                    normal = nodesFine_[no1]->getInnerNormal();

                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];

                    double signaltest = 0.0;
                    for (int i = 0; i < 2; i++)
                    {
                        signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)
                        signal = 1.;

                    dist2 *= signal;
                };
                if (fabs(dist2) < fabs(dist))
                    dist = dist2;

            }; // if bf is the blend boundary
        };     // numboundaryfine

        if (fabs(nodesCoarse_[ino]->getDistFunction()) < 1.e-2)
        {

            nodesCoarse_[ino]->setDistFunction(dist);
        };

    }; // numnodescoarse

    // Fine mesh - closer distance for nodes or control points from defined fine boundary
    for (int ino = 0; ino < numNodesFine; ino++)
    {

        double &alpha_f = parametersFine->getAlphaF();

        double *xx = nodesFine_[ino]->getCoordinates();

        for (int i = 0; i < dim; i++)
            x[i] = xx[i];

        // double *xx = nodesFine_[ino]->getCoordinates();
        // double *xxp = nodesFine_[ino]->getPreviousCoordinates();
        // for (int i = 0; i < dim; i++) x[i] = alpha_f * xx[i] + (1. - alpha_f) * xxp[i];

        dist = 10000000000000000000000000000.;

        for (int ibound = 0; ibound < numBoundElemFine; ibound++)
        {

            if (boundaryFine_[ibound]->getConstrain(0) == 2)
            {

                int *connec = elementsFine_[boundaryFine_[ibound]->getElement()]->getConnectivity();

                std::pair<std::vector<int>, std::vector<int>> elemBound;
                elemBound = elementsFine_[boundaryFine_[ibound]->getElement()]->getElemSideInBoundary();
                int numofBoundaries = elemBound.first.size();
                int side;
                for (int j = 0; j < numofBoundaries; j++)
                {
                    if (elemBound.second[j] == ibound)
                    {
                        side = elemBound.first[j];
                    };
                };

                if (side == 0)
                {
                    bconnec[0] = connec[1];
                    bconnec[1] = connec[2];
                    bconnec[2] = connec[4];
                };
                if (side == 1)
                {
                    bconnec[0] = connec[2];
                    bconnec[1] = connec[0];
                    bconnec[2] = connec[5];
                };
                if (side == 2)
                {
                    bconnec[0] = connec[0];
                    bconnec[1] = connec[1];
                    bconnec[2] = connec[3];
                };

                // first segment
                int no1, no2;
                no1 = bconnec[0];
                no2 = bconnec[2];

                double *xx1 = nodesFine_[no1]->getCoordinates();
                double *xx2 = nodesFine_[no2]->getCoordinates();
                // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

                for (int k = 0; k < dim; k++)
                {
                    x1[k] = xx1[k];
                    x2[k] = xx2[k];
                    // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                    // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
                }

                double aux0 = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                   (x2[0] - x1[0]) * (x2[0] - x1[0]));
                double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0]) +
                               (x[1] - x1[1]) * (x2[1] - x1[1])) /
                              aux0;
                double dist2 = -((x2[1] - x1[1]) * x[0] -
                                 (x2[0] - x1[0]) * x[1] +
                                 x2[0] * x1[1] - x2[1] * x1[0]) /
                               aux0;

                if (aux1 > aux0)
                {

                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    // find signal
                    // side normal vector
                    double *normal = nodesFine_[no2]->getInnerNormal();

                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];

                    double signaltest = 0.0;
                    for (int i = 0; i < 2; i++)
                    {
                        signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)
                        signal = 1.;

                    dist2 *= signal;
                };

                if (aux1 < 0.)
                {

                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    // find signal
                    // side normal vector
                    double *normal = nodesFine_[no1]->getInnerNormal();

                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];

                    double signaltest = 0.0;
                    for (int i = 0; i < 2; i++)
                    {
                        signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)
                        signal = 1.;

                    dist2 *= signal;
                };

                if (fabs(dist2) < fabs(dist))
                    dist = dist2;

                // second segment
                no1 = bconnec[2];
                no2 = bconnec[1];

                xx1 = nodesFine_[no1]->getCoordinates();
                xx2 = nodesFine_[no2]->getCoordinates();
                // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

                for (int k = 0; k < dim; k++)
                {
                    x1[k] = xx1[k];
                    x2[k] = xx2[k];
                    // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                    // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
                }

                aux0 = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                            (x2[0] - x1[0]) * (x2[0] - x1[0]));
                aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0]) +
                        (x[1] - x1[1]) * (x2[1] - x1[1])) /
                       aux0;
                dist2 = -((x2[1] - x1[1]) * x[0] - (x2[0] - x1[0]) * x[1] +
                          x2[0] * x1[1] - x2[1] * x1[0]) /
                        aux0;

                if (aux1 > aux0)
                {
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));

                    double *normal = nodesFine_[no2]->getInnerNormal();

                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];

                    double signaltest = 0.0;
                    for (int i = 0; i < 2; i++)
                    {
                        signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)
                        signal = 1.;

                    dist2 *= signal;
                };

                if (aux1 < 0.)
                {
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));

                    // find signal
                    // side normal vector
                    double *normal;
                    normal = nodesFine_[no1]->getInnerNormal();

                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];

                    double signaltest = 0.0;
                    for (int i = 0; i < 2; i++)
                    {
                        signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)
                        signal = 1.;

                    dist2 *= signal;
                };
                if (fabs(dist2) < fabs(dist))
                    dist = dist2;

            }; // if bf is the blend boundary
        };     // numboundfine

        if (dist < 0)
            dist = 0;
        nodesFine_[ino]->setDistFunction(dist);

    }; // numNodes

}; // função

template <int DIM>
void Arlequin<DIM>::searchMinDist_FEM(double *pointCoord, int &ibound, double &mindist)
{

    // element connectivity
    int *connec = elementsFine_[boundaryFine_[ibound]->getElement()]->getConnectivity();

    // side of the element in the boundary
    std::pair<std::vector<int>, std::vector<int>> elemBound;
    elemBound = elementsFine_[boundaryFine_[ibound]->getElement()]->getElemSideInBoundary();
    int numofBoundaries = elemBound.first.size();
    int side;
    for (int j = 0; j < numofBoundaries; j++)
    {
        if (elemBound.second[j] == ibound)
        {
            side = elemBound.first[j];
        };
    };

    // FEM surface element
    int bconnec[6];
    if (side == 0)
    {
        bconnec[0] = connec[3];
        bconnec[1] = connec[1];
        bconnec[2] = connec[2];
        bconnec[3] = connec[9];
        bconnec[4] = connec[5];
        bconnec[5] = connec[8];
    };
    if (side == 1)
    {
        bconnec[0] = connec[3];
        bconnec[1] = connec[2];
        bconnec[2] = connec[0];
        bconnec[3] = connec[8];
        bconnec[4] = connec[6];
        bconnec[5] = connec[7];
    };
    if (side == 2)
    {
        bconnec[0] = connec[3];
        bconnec[1] = connec[0];
        bconnec[2] = connec[1];
        bconnec[3] = connec[7];
        bconnec[4] = connec[4];
        bconnec[5] = connec[9];
    };
    if (side == 3)
    {
        bconnec[0] = connec[2];
        bconnec[1] = connec[1];
        bconnec[2] = connec[0];
        bconnec[3] = connec[5];
        bconnec[4] = connec[4];
        bconnec[5] = connec[6];
    };

    // Computing mininum distance to the surface
    BoundShapeFunction<3> boundShape;
    double phi_[6], xsi[2] = {};

    double **dphi_dxsi;
    dphi_dxsi = new double *[2];
    for (int i = 0; i < 2; ++i)
        dphi_dxsi[i] = new double[6];

    double ***ddphi_dxsi;
    ddphi_dxsi = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddphi_dxsi[i] = new double *[2];
        for (int j = 0; j < 2; j++)
            ddphi_dxsi[i][j] = new double[6];
    };

    double normError = 1.0e06;
    int cont = 1;
    double tolerance = 1.e-8;

    while (normError >= tolerance && cont <= 20)
    {

        boundShape.evaluateBoundaryFem(xsi, phi_);
        boundShape.evaluateGradientBoundaryFem(xsi, dphi_dxsi);
        boundShape.evaluateHessianBoundaryFem(ddphi_dxsi);

        double surfCoord[3] = {}, t1[3] = {}, t2[3] = {}, dt1_dxsi1[3] = {},
               dt1_dxsi2[3] = {}, dt2_dxsi2[3] = {};

        for (int i = 0; i < 6; i++)
        {

            double *coord = nodesFine_[bconnec[i]]->getCoordinates();

            for (int j = 0; j < 3; j++)
            {

                surfCoord[j] += phi_[i] * coord[j];
                t1[j] += dphi_dxsi[0][i] * coord[j];
                t2[j] += dphi_dxsi[1][i] * coord[j];
                dt1_dxsi1[j] += ddphi_dxsi[0][0][i] * coord[j];
                dt1_dxsi2[j] += ddphi_dxsi[0][1][i] * coord[j];
                dt2_dxsi2[j] += ddphi_dxsi[1][1][i] * coord[j];
            };
        };

        double g[3];
        for (int i = 0; i < 3; i++)
            g[i] = pointCoord[i] - surfCoord[i];

        double r[2] = {};
        double a[2][2] = {};
        double gt[2][2] = {};
        for (int i = 0; i < 3; i++)
        {
            r[0] += g[i] * t1[i];
            r[1] += g[i] * t2[i];
            a[0][0] += t1[i] * t1[i];
            a[0][1] += t1[i] * t2[i];
            a[1][1] += t2[i] * t2[i];
            gt[0][0] += g[i] * dt1_dxsi1[i];
            gt[0][1] += g[i] * dt1_dxsi2[i];
            gt[1][1] += g[i] * dt2_dxsi2[i];
        };
        a[1][0] = a[0][1];
        gt[1][0] = gt[0][1];

        double m[2][2];
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                m[i][j] = a[i][j] - gt[i][j];
            };
        };

        // inverting M
        double detm = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        double invM[2][2];
        invM[0][0] = m[1][1] / detm;
        invM[0][1] = -m[0][1] / detm;
        invM[1][0] = -m[1][0] / detm;
        invM[1][1] = m[0][0] / detm;

        double deltaXsi[2] = {};
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                deltaXsi[i] += invM[i][j] * r[j];
            };
        };

        normError = 0.;
        for (int i = 0; i < 2; i++)
        {
            xsi[i] += deltaXsi[i];
            normError += deltaXsi[i] * deltaXsi[i];
        }

        normError = sqrt(normError);

        cont += 1;
    };

    if (xsi[0] < 0.)
        xsi[0] = 0.;
    else if (xsi[0] > 1.)
        xsi[0] = 1.;

    if (xsi[1] < 0.)
        xsi[1] = 0.;
    else if (xsi[1] > 1.)
        xsi[1] = 1.;

    // Computing the phisical distance

    double surfCoord[3] = {}, t1[3] = {}, t2[3] = {};

    boundShape.evaluateBoundaryFem(xsi, phi_);
    boundShape.evaluateGradientBoundaryFem(xsi, dphi_dxsi);

    for (int i = 0; i < 6; i++)
    {
        double *coord = nodesFine_[bconnec[i]]->getCoordinates();
        for (int j = 0; j < 3; j++)
        {
            surfCoord[j] += phi_[i] * coord[j];
            t1[j] += dphi_dxsi[0][i] * coord[j];
            t2[j] += dphi_dxsi[1][i] * coord[j];
        };
    };

    double normalVec[3]; // vectorial product
    normalVec[0] = t1[1] * t2[2] - t1[2] * t2[1];
    normalVec[1] = t1[2] * t2[0] - t1[0] * t2[2];
    normalVec[2] = t1[0] * t2[1] - t1[1] * t2[0];


    double g[3];
    for (int i = 0; i < 3; i++)
        g[i] = pointCoord[i] - surfCoord[i];


    mindist = 0.;
    for (int i = 0; i < 3; i++)
        mindist += g[i] * g[i];
    mindist = sqrt(mindist);

    double tfactor = 0.;
    for (int i = 0; i < 3; i++)
        tfactor += normalVec[i] * g[i];


    double factor = (tfactor >= 0.0) ? -1.0 : 1.0;

    mindist *= factor;

    for (int i = 0; i < 2; ++i)
        delete[] dphi_dxsi[i];
    delete[] dphi_dxsi;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] ddphi_dxsi[i][j];
        };
        delete[] ddphi_dxsi[i];
    };
    delete[] ddphi_dxsi;
};

template <>
void Arlequin<2>::searchMinDist_ISO(double *pointCoord, int &ibound, double &mindist)
{

}

template <>
void Arlequin<3>::searchMinDist_ISO(double *pointCoord, int &ibound, double &mindist)
{

    // element connectivity
    int *connec = elementsFine_[boundaryFine_[ibound]->getElement()]->getConnectivity();
    int npatch_ = elementsFine_[boundaryFine_[ibound]->getElement()]->getPatch();
   
    // side of the element in the boundary
    std::pair<std::vector<int>, std::vector<int>> elemBound;
    elemBound = elementsFine_[boundaryFine_[ibound]->getElement()]->getElemSideInBoundary();
    int numofBoundaries = elemBound.first.size();
    int side;
    for (int j = 0; j < numofBoundaries; j++)
    {
        if (elemBound.second[j] == ibound)
        {
            side = elemBound.first[j];
        };
    };

    // IGA surface element
   int bconnec[9];
   double wpc[9];
   int *inc_;
   if(side == 0){
        
        bconnec[0] = connec[0];
        bconnec[1] = connec[1];
        bconnec[2] = connec[2];
        bconnec[3] = connec[3];
        bconnec[4] = connec[4];
        bconnec[5] = connec[5];
        bconnec[6] = connec[6];
        bconnec[7] = connec[7];
        bconnec[8] = connec[8];

		wpc[0] = nodesFine_[connec[0]] -> getWeightPC();
		wpc[1] = nodesFine_[connec[1]] -> getWeightPC();
		wpc[2] = nodesFine_[connec[2]] -> getWeightPC();
		wpc[3] = nodesFine_[connec[3]] -> getWeightPC();
		wpc[4] = nodesFine_[connec[4]] -> getWeightPC();
		wpc[5] = nodesFine_[connec[5]] -> getWeightPC();
		wpc[6] = nodesFine_[connec[6]] -> getWeightPC();
		wpc[7] = nodesFine_[connec[7]] -> getWeightPC();
		wpc[8] = nodesFine_[connec[8]] -> getWeightPC();  

		inc_ = nodesFine_[connec[8]] -> getINC(); 
	};
	
	if(side  == 1){
		
        bconnec[0] = connec[2];
        bconnec[1] = connec[11];
        bconnec[2] = connec[20];
        bconnec[3] = connec[5];
        bconnec[4] = connec[14];
        bconnec[5] = connec[23];
        bconnec[6] = connec[8];
        bconnec[7] = connec[17];
        bconnec[8] = connec[26];

		wpc[0] = nodesFine_[connec[2]] -> getWeightPC();
		wpc[1] = nodesFine_[connec[11]] -> getWeightPC();
		wpc[2] = nodesFine_[connec[20]] -> getWeightPC();
		wpc[3] = nodesFine_[connec[5]] -> getWeightPC();
		wpc[4] = nodesFine_[connec[14]] -> getWeightPC();
		wpc[5] = nodesFine_[connec[23]] -> getWeightPC();
		wpc[6] = nodesFine_[connec[8]] -> getWeightPC();
		wpc[7] = nodesFine_[connec[17]] -> getWeightPC();
		wpc[8] = nodesFine_[connec[26]] -> getWeightPC();  
        
		inc_ = nodesFine_[connec[26]] -> getINC(); 
	};

	if(side == 2){

        bconnec[0] = connec[20];
        bconnec[1] = connec[19];
        bconnec[2] = connec[18];
        bconnec[3] = connec[23];
        bconnec[4] = connec[22];
        bconnec[5] = connec[21];
        bconnec[6] = connec[26];
        bconnec[7] = connec[25];
        bconnec[8] = connec[24];

		wpc[0] = nodesFine_[connec[20]] -> getWeightPC();
		wpc[1] = nodesFine_[connec[19]] -> getWeightPC();
		wpc[2] = nodesFine_[connec[18]] -> getWeightPC();
		wpc[3] = nodesFine_[connec[23]] -> getWeightPC();
		wpc[4] = nodesFine_[connec[22]] -> getWeightPC();
		wpc[5] = nodesFine_[connec[21]] -> getWeightPC();
		wpc[6] = nodesFine_[connec[26]] -> getWeightPC();
		wpc[7] = nodesFine_[connec[25]] -> getWeightPC();
		wpc[8] = nodesFine_[connec[24]] -> getWeightPC();  
        
		inc_ = nodesFine_[connec[26]] -> getINC(); 

	};
		

	 if(side == 3){  
		
        bconnec[0] = connec[18];
        bconnec[1] = connec[9];
        bconnec[2] = connec[0];
        bconnec[3] = connec[21];
        bconnec[4] = connec[12];
        bconnec[5] = connec[3];
        bconnec[6] = connec[24];
        bconnec[7] = connec[15];
        bconnec[8] = connec[6];

		wpc[0] = nodesFine_[connec[18]] -> getWeightPC();
		wpc[1] = nodesFine_[connec[9]] -> getWeightPC();
		wpc[2] = nodesFine_[connec[0]] -> getWeightPC();
		wpc[3] = nodesFine_[connec[21]] -> getWeightPC();
		wpc[4] = nodesFine_[connec[12]] -> getWeightPC();
		wpc[5] = nodesFine_[connec[3]] -> getWeightPC();
		wpc[6] = nodesFine_[connec[24]] -> getWeightPC();
		wpc[7] = nodesFine_[connec[15]] -> getWeightPC();
		wpc[8] = nodesFine_[connec[6]] -> getWeightPC();  
        
		inc_ = nodesFine_[connec[24]] -> getINC(); 
	};

	if(side == 4){   

        bconnec[0] = connec[18];
        bconnec[1] = connec[19];
        bconnec[2] = connec[20];
        bconnec[3] = connec[9];
        bconnec[4] = connec[10];
        bconnec[5] = connec[11];
        bconnec[6] = connec[0];
        bconnec[7] = connec[1];
        bconnec[8] = connec[2];

		wpc[0] = nodesFine_[connec[18]] -> getWeightPC();
		wpc[1] = nodesFine_[connec[19]] -> getWeightPC();
		wpc[2] = nodesFine_[connec[20]] -> getWeightPC();
		wpc[3] = nodesFine_[connec[9]] -> getWeightPC();
		wpc[4] = nodesFine_[connec[10]] -> getWeightPC();
		wpc[5] = nodesFine_[connec[11]] -> getWeightPC();
		wpc[6] = nodesFine_[connec[0]] -> getWeightPC();
		wpc[7] = nodesFine_[connec[1]] -> getWeightPC();
		wpc[8] = nodesFine_[connec[2]] -> getWeightPC();  
        
		inc_ = nodesFine_[connec[20]] -> getINC();  
	};

	if(side == 5){  

        bconnec[0] = connec[6];
        bconnec[1] = connec[7];
        bconnec[2] = connec[8];
        bconnec[3] = connec[15];
        bconnec[4] = connec[16];
        bconnec[5] = connec[17];
        bconnec[6] = connec[24];
        bconnec[7] = connec[25];
        bconnec[8] = connec[26];

		wpc[0] = nodesFine_[connec[6]] -> getWeightPC();
		wpc[1] = nodesFine_[connec[7]] -> getWeightPC();
		wpc[2] = nodesFine_[connec[8]] -> getWeightPC();
		wpc[3] = nodesFine_[connec[15]] -> getWeightPC();
		wpc[4] = nodesFine_[connec[16]] -> getWeightPC();
		wpc[5] = nodesFine_[connec[17]] -> getWeightPC();
		wpc[6] = nodesFine_[connec[24]] -> getWeightPC();
		wpc[7] = nodesFine_[connec[25]] -> getWeightPC();
		wpc[8] = nodesFine_[connec[26]] -> getWeightPC();  
        
		inc_ = nodesFine_[connec[26]] -> getINC();   

	};

    // Computing mininum distance to the surface
    BoundShapeFunction<3> boundShape;
    double phi_[9], xsi[2] = {};

    double **dphi_dxsi;
    dphi_dxsi = new double *[2];
    for (int i = 0; i < 2; ++i)
        dphi_dxsi[i] = new double[9];

    double ***ddphi_ddxsi;
    ddphi_ddxsi = new double **[2];
    for (int i = 0; i < 2; ++i)
    {
        ddphi_ddxsi[i] = new double *[2];
        for (int j = 0; j < 2; j++)
            ddphi_ddxsi[i][j] = new double[9];
    };

    double normError = 1.0e06;
    int cont = 1;
    double tolerance = 1.e-8;


    while (normError >= tolerance && cont <= 20)
    {

        boundShape.evaluateBoundaryIso(xsi,phi_,wpc,side,inc_,IsoParFine,npatch_);
        boundShape.evaluateGradientBoundaryIso(xsi,dphi_dxsi,wpc,side,inc_,IsoParFine,npatch_);
        boundShape.evaluateHessianBoundaryIso(xsi,ddphi_ddxsi,wpc,side,inc_,IsoParFine,npatch_);

        double surfCoord[3] = {}, t1[3] = {}, t2[3] = {}, dt1_dxsi1[3] = {},
               dt1_dxsi2[3] = {}, dt2_dxsi2[3] = {};

        for (int i = 0; i < 9; i++)
        {

            double *coord = nodesFine_[bconnec[i]]->getCoordinates();

            for (int j = 0; j < 3; j++)
            {

                surfCoord[j] += phi_[i] * coord[j]/wpc[i];
                t1[j] += dphi_dxsi[0][i] * coord[j]/wpc[i];
                t2[j] += dphi_dxsi[1][i] * coord[j]/wpc[i];
                dt1_dxsi1[j] += ddphi_ddxsi[0][0][i] * coord[j]/wpc[i];
                dt1_dxsi2[j] += ddphi_ddxsi[0][1][i] * coord[j]/wpc[i];
                dt2_dxsi2[j] += ddphi_ddxsi[1][1][i] * coord[j]/wpc[i];
            };
        };

        double g[3];
        for (int i = 0; i < 3; i++)
            g[i] = pointCoord[i] - surfCoord[i];

        double r[2] = {};
        double a[2][2] = {};
        double gt[2][2] = {};
        for (int i = 0; i < 3; i++)
        {
            r[0] += g[i] * t1[i];
            r[1] += g[i] * t2[i];
            a[0][0] += t1[i] * t1[i];
            a[0][1] += t1[i] * t2[i];
            a[1][1] += t2[i] * t2[i];
            gt[0][0] += g[i] * dt1_dxsi1[i];
            gt[0][1] += g[i] * dt1_dxsi2[i];
            gt[1][1] += g[i] * dt2_dxsi2[i];
        };
        a[1][0] = a[0][1];
        gt[1][0] = gt[0][1];

        double m[2][2];
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                m[i][j] = a[i][j] - gt[i][j];
            };
        };

        // inverting M
        double detm = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        double invM[2][2];
        invM[0][0] = m[1][1] / detm;
        invM[0][1] = -m[0][1] / detm;
        invM[1][0] = -m[1][0] / detm;
        invM[1][1] = m[0][0] / detm;

        double deltaXsi[2] = {};
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                deltaXsi[i] += invM[i][j] * r[j];
            };
        };

        normError = 0.;
        for (int i = 0; i < 2; i++)
        {
            xsi[i] += deltaXsi[i];
            normError += deltaXsi[i] * deltaXsi[i];
        }

        normError = sqrt(normError);

        cont += 1;
    };

    if (xsi[0] < -1.)
        xsi[0] = -1.;
    else if (xsi[0] > 1.)
        xsi[0] = 1.;

    if (xsi[1] < -1.)
        xsi[1] = -1.;
    else if (xsi[1] > 1.)
        xsi[1] = 1.;

    // Computing the phisical distance

    double surfCoord[3] = {}, t1[3] = {}, t2[3] = {};

    boundShape.evaluateBoundaryIso(xsi,phi_,wpc,side,inc_,IsoParFine,npatch_);
    boundShape.evaluateGradientBoundaryIso(xsi,dphi_dxsi,wpc,side,inc_,IsoParFine,npatch_);

    for (int i = 0; i < 9; i++)
    {
        double *coord = nodesFine_[bconnec[i]]->getCoordinates();
        for (int j = 0; j < 3; j++)
        {
            surfCoord[j] += phi_[i] * coord[j]/wpc[i];
            t1[j] += dphi_dxsi[0][i] * coord[j]/wpc[i];
            t2[j] += dphi_dxsi[1][i] * coord[j]/wpc[i];
        };
    };

    double normalVec[3]; // vectorial product
    normalVec[0] = t1[1] * t2[2] - t1[2] * t2[1];
    normalVec[1] = t1[2] * t2[0] - t1[0] * t2[2];
    normalVec[2] = t1[0] * t2[1] - t1[1] * t2[0];


    double g[3];
    for (int i = 0; i < 3; i++)
        g[i] = pointCoord[i] - surfCoord[i];

    mindist = 0.;
    for (int i = 0; i < 3; i++)
        mindist += g[i] * g[i];
    mindist = sqrt(mindist);

    double tfactor = 0.;
    for (int i = 0; i < 3; i++)
        tfactor += normalVec[i] * g[i];

    double factor = (tfactor >= 0.0) ? -1.0 : 1.0;

    mindist *= factor;

    for (int i = 0; i < 2; ++i)
        delete[] dphi_dxsi[i];
    delete[] dphi_dxsi;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; j++)
        {
            delete[] ddphi_ddxsi[i][j];
        };
        delete[] ddphi_ddxsi[i];
    };
    delete[] ddphi_ddxsi;
};

template <int DIM>
void Arlequin<DIM>::setSignaledDistance_ISO_ISO()
{

    for (int inode = 0; inode < numNodesFine; inode++)
    {
        nodesFine_[inode]->setDistFunction(0.0);
    };

    for (int inode = 0; inode < numNodesCoarse; inode++)
    {
        nodesCoarse_[inode]->setDistFunction(0.0);
    };

    double dist;

    // Fine mesh - closer distance for nodes or control points from defined fine boundary
    for (int ino = 0; ino < numNodesFine; ino++)
    {
        double *xx = nodesFine_[ino]->getCoordinates();

        // double x[3];
        // double &alpha_f = parametersFine -> getAlphaF();
        // double *xx = nodesFine_[ino]->getCoordinates();
        // double *xxp = nodesFine_[ino]->getPreviousCoordinates();
        // for (int i = 0; i < 3; i++) x[i] = alpha_f * xx[i] + (1. - alpha_f) * xxp[i];
        // dist = 10000000000000000000000000000.;
        // for (int ibound = 0; ibound < numBoundElemFine; ibound++)
        // {
        //     if (boundaryFine_[ibound]->getConstrain(0) == 2)
        //     {
        //         double mindist;
        //         searchMinDist_ISO(x, ibound, mindist);
        //         if (fabs(mindist) < fabs(dist))
        //             dist = mindist;
                
        //     };
        // };

        dist = (0.5 - xx[1]);

        // if (dist < 0) {
        //     dist = 0;
        // };

        nodesFine_[ino]->setDistFunction(dist);
    };

    // Coarse mesh - closer distance for nodes or control points from defined fine boundary
    for (int ino = 0; ino < numNodesCoarse; ino++)
    {
        double *xx = nodesCoarse_[ino]->getCoordinates();
        // double x[3];
        // double &alpha_f = parametersFine -> getAlphaF();
        // double *xx = nodesFine_[ino]->getCoordinates();
        // double *xxp = nodesFine_[ino]->getPreviousCoordinates();
        // for (int i = 0; i < 3; i++) x[i] = alpha_f * xx[i] + (1. - alpha_f) * xxp[i];
        // dist = 10000000000000000000000000000.;
        // for (int ibound = 0; ibound < numBoundElemFine; ibound++)
        // {
        //     if (boundaryFine_[ibound]->getConstrain(0) == 2)
        //     {
        //         double mindist;
        //         searchMinDist_ISO(x, ibound, mindist);
        //         if (fabs(mindist) < fabs(dist))
        //             dist = mindist;
        //     };
        // };

        dist = (0.5 - xx[1]);

        nodesCoarse_[ino]->setDistFunction(dist);
    };
};


template <>
void Arlequin<3>::setSignaledDistance_FEM_ISO()
{

    for (int inode = 0; inode < numNodesFine; inode++)
    {
        nodesFine_[inode]->setDistFunction(0.0);
    };

    for (int inode = 0; inode < numNodesCoarse; inode++)
    {
        nodesCoarse_[inode]->setDistFunction(0.0);
    };

    double dist;

    // Fine mesh - closer distance for nodes or control points from defined fine boundary
    for (int ino = 0; ino < numNodesFine; ino++)
    {
        double x[3];
        double *xx = nodesFine_[ino]->getCoordinates();
        for (int i = 0; i < 3; i++)
            x[i] = xx[i];
        // double x[3];
        // double &alpha_f = parametersFine -> getAlphaF();
        // double *xx = nodesFine_[ino]->getCoordinates();
        // double *xxp = nodesFine_[ino]->getPreviousCoordinates();
        // for (int i = 0; i < 3; i++) x[i] = alpha_f * xx[i] + (1. - alpha_f) * xxp[i];
        dist = 10000000000000000000000000000.;
        for (int ibound = 0; ibound < numBoundElemFine; ibound++)
        {
            if (boundaryFine_[ibound]->getConstrain(0) == 2)
            {
                double mindist;
                searchMinDist_FEM(x, ibound, mindist);
                if (fabs(mindist) < fabs(dist))
                    dist = mindist;
                
            };
        };

        if (dist < 0) {
            dist = -dist;
            //dist = 0;
        };

        nodesFine_[ino]->setDistFunction(dist);
    };

    // Coarse mesh - closer distance for nodes or control points from defined fine boundary
    for (int ino = 0; ino < numNodesCoarse; ino++)
    {
        double x[3];
        double *xx = nodesCoarse_[ino]->getCoordinates();
        for (int i = 0; i < 3; i++)
            x[i] = xx[i];
        // double x[3];
        // double &alpha_f = parametersFine -> getAlphaF();
        // double *xx = nodesFine_[ino]->getCoordinates();
        // double *xxp = nodesFine_[ino]->getPreviousCoordinates();
        // for (int i = 0; i < 3; i++) x[i] = alpha_f * xx[i] + (1. - alpha_f) * xxp[i];
        dist = 10000000000000000000000000000.;
        for (int ibound = 0; ibound < numBoundElemFine; ibound++)
        {
            if (boundaryFine_[ibound]->getConstrain(0) == 2)
            {
                double mindist;
                searchMinDist_FEM(x, ibound, mindist);
                if (fabs(mindist) < fabs(dist))
                    dist = mindist;
            };
        };

        nodesCoarse_[ino]->setDistFunction(dist);
    };
};


template <int DIM>
void Arlequin<DIM>::setSignaledDistance_FEM_FEM()
{

    for (int inode = 0; inode < numNodesFine; inode++)
    {
        nodesFine_[inode]->setDistFunction(0.0);
    };

    for (int inode = 0; inode < numNodesCoarse; inode++)
    {
        nodesCoarse_[inode]->setDistFunction(0.0);
    };

    double dist;

    // Fine mesh - closer distance for nodes from defined fine boundary
    for (int ino = 0; ino < numNodesFine; ino++)
    {
        double *xx = nodesFine_[ino]->getCoordinates();
        dist = (0.5 - xx[1]);
        nodesFine_[ino]->setDistFunction(dist);
    };

    // Coarse mesh - closer distance for nodes or control points from defined fine boundary
    for (int ino = 0; ino < numNodesCoarse; ino++)
    {
        double *xx = nodesCoarse_[ino]->getCoordinates();
        dist = (0.5 - xx[1]);
        nodesCoarse_[ino]->setDistFunction(dist);
    };



    // for (int inode = 0; inode < numNodesFine; inode++)
    // {
    //     nodesFine_[inode]->setDistFunction(0.0);
    // };

    // for (int inode = 0; inode < numNodesCoarse; inode++)
    // {
    //     nodesCoarse_[inode]->setDistFunction(0.0);
    // };

    // double dist;

    // // Fine mesh - closer distance for nodes from defined fine boundary
    // for (int ino = 0; ino < numNodesFine; ino++)
    // {
    //     double x[DIM];
    //     double *xx = nodesFine_[ino]->getCoordinates();
    //     for (int i = 0; i < DIM; i++) x[i] = xx[i];

    //     dist = 10000000000000000000000000000.;
    //     for (int ibound = 0; ibound < numBoundElemFine; ibound++)
    //     {
    //         if (boundaryFine_[ibound]->getConstrain(0) == 2)
    //         {
    //             double mindist;
    //             searchMinDist_FEM(x, ibound, mindist);
    //             if (fabs(mindist) < fabs(dist))
    //                 dist = mindist;
                
    //         };
    //     };

    //     if (dist < 0) {
    //         dist = -dist;
    //         //dist = 0;
    //     };

    //     nodesFine_[ino]->setDistFunction(dist);
    // };

    // // Coarse mesh - closer distance for nodes or control points from defined fine boundary
    // for (int ino = 0; ino < numNodesCoarse; ino++)
    // {
    //     double x[DIM];
    //     double *xx = nodesCoarse_[ino]->getCoordinates();
    //     for (int i = 0; i < DIM; i++) x[i] = xx[i];

    //     dist = 10000000000000000000000000000.;
    //     for (int ibound = 0; ibound < numBoundElemFine; ibound++)
    //     {
    //         if (boundaryFine_[ibound]->getConstrain(0) == 2)
    //         {
    //             double mindist;
    //             searchMinDist_FEM(x, ibound, mindist);
    //             if (fabs(mindist) < fabs(dist))
    //                 dist = mindist;
    //         };
    //     };

    //     nodesCoarse_[ino]->setDistFunction(dist);
    // };
};


//------------------------------------------------------------------------------
//------------------------------SETS GLUING ZONE--------------------------------
//------------------------------------------------------------------------------
template <int DIM>
void Arlequin<DIM>::setGluingZone_FEM_ISO()
{

    double glueZoneThickness = fineModel.glueZoneThickness;

    int LNN = 4*DIM-2;
    int LNNC = 18*DIM-27;

    int flag;
    int nodesCZ[numNodesFine];
    int nodesCZ2[numNodesCoarse];
    QuadShapeFunction<DIM> shapeQuad;

    for (int i = 0; i < numNodesFine; i++)
        nodesCZ[i] = 0;
    elementsGlueZoneFine_.reserve(numElemFine);
    glueZoneFine_.reserve(numElemFine);
    nodesGlueZoneFine_.reserve(numNodesFine);

    for (int i = 0; i < numNodesCoarse; i++)
        nodesCZ2[i] = 0;
    elementsGlueZoneCoarse_.reserve(numElemCoarse);
    nodesGlueZoneCoarse_.reserve(numNodesCoarse);

    // Defines a criterion to select the fine FEM elements that are in the gluing zone
    int index = 0;
    for (int jel = 0; jel < numElemFine; jel++)
    {

        int *connec = elementsFine_[jel]->getConnectivity();

        flag = 0;
        for (int ino = 0; ino < LNN; ino++)
        {
            double dist = nodesFine_[connec[ino]]->getDistFunction();
            if (dist <= glueZoneThickness + 0.00001)
            {
                flag += 1;
            };
        };

        if (flag == LNN)
        {
            
            elementsGlueZoneFine_.push_back(jel);
            elementsFine_[jel]->setGlueZone();
            GlueZone *el = new GlueZone(index++, jel);
            glueZoneFine_.push_back(el);
        };
    };

    // Defines which nodes are in the gluing zone
    numElemGlueZoneFine = elementsGlueZoneFine_.size();
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {
        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();
        for (int ino = 0; ino < LNN; ino++)
        {
            nodesCZ[connec[ino]] += 1;
            nodesFine_[connec[ino]]->setGlueZone();
        };
    };

    // Compute and save number of nodes in the gluing zone
    numNodesGlueZoneFine = 0;
    for (int i = 0; i < numNodesFine; i++)
    {
        if (nodesCZ[i] > 0)
        {
            numNodesGlueZoneFine += 1;
            nodesGlueZoneFine_.push_back(i);
        };
    };

    NCnumNodesGlueZoneFine = numNodesGlueZoneFine;

    // Defining the Lagrange multipliers in the fine mesh
    for (int i = 0; i < numNodesGlueZoneFine; i++)
    {
        double *x = nodesFine_[nodesGlueZoneFine_[i]]->getCoordinates();
        Nodes *no = new Nodes(x, i, 1.);
        nodesLagrangeFine_.push_back(no);
    };

    // Define the Lagrange Multipliers connectivity in the fine mesh
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        int *connecAux = new int[LNN];
        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        for (int ino = 0; ino < numNodesGlueZoneFine; ino++)
            for (int k = 0; k < LNN; k++)
                if (nodesGlueZoneFine_[ino] == connec[k])
                    connecAux[k] = ino;

        glueZoneFine_[i]->setConnectivity(connecAux);
    };

    // Defines a criterion to select the coarse elements that are in the gluing zone
    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        int *connec = elementsCoarse_[jel]->getConnectivity();

        flag = 0;

        double distance_[LNNC], Bdistance_[LNNC] = {};

        for (int i = 0; i < LNNC; i++)
            distance_[i] = nodesCoarse_[connec[i]]->getDistFunction();

        for (int icp = 0; icp < LNNC; icp++)
        {

            double wpc[LNNC], phi_[LNNC], xsi[DIM];

            for (int i = 0; i < LNNC; i++)
                wpc[i] = nodesCoarse_[connec[i]]->getWeightPC();
            int *inc_ = nodesCoarse_[connec[LNNC - 1]]->getINC();
            int patch_ = elementsCoarse_[jel]->getPatch();
            for (int i = 0; i < DIM; i++)
                xsi[i] = shapeQuad.ParametricCoordBezier(icp, i);

            shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParCoarse, patch_);

            for (int i = 0; i < LNNC; i++)
                Bdistance_[icp] += phi_[i] * distance_[i];
        };

        for (int ino = 0; ino < LNNC; ino++)
        {
            if ((Bdistance_[ino] <= glueZoneThickness + 0.00001) && (Bdistance_[ino] >= 0.00001))
            {
                flag += 1;
            };
        };

        if (flag != 0)
        {
            elementsGlueZoneCoarse_.push_back(jel);
            elementsCoarse_[jel]->setGlueZone();
        };
    };

    // Defines which coarse nodes are in the gluing zone
    numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    for (int i = 0; i < numElemGlueZoneCoarse; i++)
    {
        int *connec = elementsCoarse_[elementsGlueZoneCoarse_[i]]->getConnectivity();
        for (int ino = 0; ino < LNNC; ino++)
        {
            nodesCZ2[connec[ino]] += 1;
            nodesCoarse_[connec[ino]]->setGlueZone();
        };
    };

    // Compute number of coarse nodes in the gluing zone
    numNodesGlueZoneCoarse = 0;
    for (int i = 0; i < numNodesCoarse; i++)
    {
        if (nodesCZ2[i] > 0)
        {
            numNodesGlueZoneCoarse += 1;
            nodesGlueZoneCoarse_.push_back(i);
        };
    };
};

template <int DIM>
void Arlequin<DIM>::setGluingZone_FEM_FEM()
{

    double glueZoneThickness = fineModel.glueZoneThickness;

    int LNN = 4*DIM-2;
    int LNNC = 4*DIM-2;

    int flag;
    int nodesCZ[numNodesFine];
    int nodesCZ2[numNodesCoarse];
    QuadShapeFunction<DIM> shapeQuad;

    for (int i = 0; i < numNodesFine; i++)
        nodesCZ[i] = 0;
    elementsGlueZoneFine_.reserve(numElemFine);
    glueZoneFine_.reserve(numElemFine);
    nodesGlueZoneFine_.reserve(numNodesFine);

    for (int i = 0; i < numNodesCoarse; i++)
        nodesCZ2[i] = 0;
    elementsGlueZoneCoarse_.reserve(numElemCoarse);
    nodesGlueZoneCoarse_.reserve(numNodesCoarse);

    // Defines a criterion to select the fine FEM elements that are in the gluing zone
    int index = 0;
    for (int jel = 0; jel < numElemFine; jel++)
    {
        int *connec = elementsFine_[jel]->getConnectivity();

        flag = 0;
        for (int ino = 0; ino < LNN; ino++)
        {
            double dist = nodesFine_[connec[ino]]->getDistFunction();
            if (dist <= glueZoneThickness + 0.00001)
            {
                flag += 1;
            };
        };

        if (flag == LNN)
        {
            elementsGlueZoneFine_.push_back(jel);
            elementsFine_[jel]->setGlueZone();
            GlueZone *el = new GlueZone(index++, jel);
            glueZoneFine_.push_back(el);
        };
    };

    // Defines which nodes are in the gluing zone
    numElemGlueZoneFine = elementsGlueZoneFine_.size();
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {
        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();
        for (int ino = 0; ino < LNN; ino++)
        {
            nodesCZ[connec[ino]] += 1;
            nodesFine_[connec[ino]]->setGlueZone();
        };
    };

    // Compute and save number of nodes in the gluing zone
    numNodesGlueZoneFine = 0;
    for (int i = 0; i < numNodesFine; i++)
    {
        if (nodesCZ[i] > 0)
        {
            numNodesGlueZoneFine += 1;
            nodesGlueZoneFine_.push_back(i);
        };
    };

    NCnumNodesGlueZoneFine = numNodesGlueZoneFine;

    // Defining the Lagrange multipliers in the fine mesh
    for (int i = 0; i < numNodesGlueZoneFine; i++)
    {
        double *x = nodesFine_[nodesGlueZoneFine_[i]]->getCoordinates();
        Nodes *no = new Nodes(x, i, 1.);
        nodesLagrangeFine_.push_back(no);
    };

    // Define the Lagrange Multipliers connectivity in the fine mesh
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        int *connecAux = new int[LNN];
        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        for (int ino = 0; ino < numNodesGlueZoneFine; ino++)
            for (int k = 0; k < LNN; k++)
                if (nodesGlueZoneFine_[ino] == connec[k])
                    connecAux[k] = ino;

        glueZoneFine_[i]->setConnectivity(connecAux);
    };

    // Defines a criterion to select the coarse elements that are in the gluing zone
    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        int *connec = elementsCoarse_[jel]->getConnectivity();

        flag = 0;
        for (int ino = 0; ino < LNNC; ino++)
        {
            double dist = nodesCoarse_[connec[ino]]->getDistFunction();
            if ((dist <= glueZoneThickness + 0.00001) && (dist >= 0.00001))
            {
                flag += 1;
            };
        };

        if (flag != 0)
        {
            elementsGlueZoneCoarse_.push_back(jel);
            elementsCoarse_[jel]->setGlueZone();
        };
    };

    // Defines which coarse nodes are in the gluing zone
    numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    for (int i = 0; i < numElemGlueZoneCoarse; i++)
    {
        int *connec = elementsCoarse_[elementsGlueZoneCoarse_[i]]->getConnectivity();
        for (int ino = 0; ino < LNNC; ino++)
        {
            nodesCZ2[connec[ino]] += 1;
            nodesCoarse_[connec[ino]]->setGlueZone();
        };
    };

    // Compute number of coarse nodes in the gluing zone
    numNodesGlueZoneCoarse = 0;
    for (int i = 0; i < numNodesCoarse; i++)
    {
        if (nodesCZ2[i] > 0)
        {
            numNodesGlueZoneCoarse += 1;
            nodesGlueZoneCoarse_.push_back(i);
        };
    };
};

template <int DIM>
void Arlequin<DIM>::setGluingZone_ISO_ISO()
{

    double glueZoneThickness = fineModel.glueZoneThickness;

    int LNN = 18*DIM-27;
    int LNNC = 18*DIM-27;

    int flag;
    int nodesCZ[numNodesFine];
    int nodesCZ2[numNodesCoarse];
    QuadShapeFunction<DIM> shapeQuad;

    for (int i = 0; i < numNodesFine; i++)
        nodesCZ[i] = 0;
    elementsGlueZoneFine_.reserve(numElemFine);
    glueZoneFine_.reserve(numElemFine);
    nodesGlueZoneFine_.reserve(numNodesFine);

    for (int i = 0; i < numNodesCoarse; i++)
        nodesCZ2[i] = 0;
    elementsGlueZoneCoarse_.reserve(numElemCoarse);
    nodesGlueZoneCoarse_.reserve(numNodesCoarse);

    // Defines a criterion to select the fine FEM elements that are in the gluing zone
    int index = 0;
    for (int jel = 0; jel < numElemFine; jel++)
    {
        int *connec = elementsFine_[jel]->getConnectivity();
        
        flag = 0;

        double distance_[LNN], Bdistance_[LNN] = {};

        for (int i = 0; i < LNN; i++)
            distance_[i] = nodesFine_[connec[i]]->getDistFunction();

        for (int icp = 0; icp < LNN; icp++)
        {

            double wpc[LNN], phi_[LNN], xsi[DIM];

            for (int i = 0; i < LNN; i++)
                wpc[i] = nodesFine_[connec[i]]->getWeightPC();
            int *inc_ = nodesFine_[connec[LNN - 1]]->getINC();
            int patch_ = elementsFine_[jel]->getPatch();
            for (int i = 0; i < DIM; i++)
                xsi[i] = shapeQuad.ParametricCoordBezier(icp, i);

            shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParFine, patch_);

            for (int i = 0; i < LNN; i++)
                Bdistance_[icp] += phi_[i] * distance_[i];

        };

        for (int i = 0; i < LNN; i++) if (Bdistance_[i] <= glueZoneThickness + 0.001) flag++; 


        if (flag == LNN)
        {
            
            elementsGlueZoneFine_.push_back(jel);
            elementsFine_[jel]->setGlueZone();
            GlueZone *el = new GlueZone(index++, jel);
            glueZoneFine_.push_back(el);
        };
    };

    // Defines which nodes are in the gluing zone
    numElemGlueZoneFine = elementsGlueZoneFine_.size();

    for (int i = 0; i < numElemGlueZoneFine; i++)
    {
        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();
        for (int ino = 0; ino < LNN; ino++)
        {
            nodesCZ[connec[ino]] += 1;
            nodesFine_[connec[ino]]->setGlueZone();
        };
    };

    // Compute and save number of nodes in the gluing zone
    numNodesGlueZoneFine = 0;
    for (int i = 0; i < numNodesFine; i++)
    {
        if (nodesCZ[i] > 0)
        {
            numNodesGlueZoneFine += 1;
            nodesGlueZoneFine_.push_back(i);
        };
    };

    NCnumNodesGlueZoneFine = numNodesGlueZoneFine;

    // Defining the Lagrange multipliers in the fine mesh
    for (int i = 0; i < numNodesGlueZoneFine; i++)
    {
        double *x = nodesFine_[nodesGlueZoneFine_[i]]->getCoordinates();
        Nodes *no = new Nodes(x, i, 1.);
        nodesLagrangeFine_.push_back(no);
    };

    // Define the Lagrange Multipliers connectivity in the fine mesh
    for (int i = 0; i < numElemGlueZoneFine; i++)
    {

        int *connecAux = new int[LNN];
        int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

        for (int ino = 0; ino < numNodesGlueZoneFine; ino++)
            for (int k = 0; k < LNN; k++)
                if (nodesGlueZoneFine_[ino] == connec[k])
                    connecAux[k] = ino;

        glueZoneFine_[i]->setConnectivity(connecAux);
    };

    // Defines a criterion to select the coarse elements that are in the gluing zone
    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        int *connec = elementsCoarse_[jel]->getConnectivity();

        flag = 0;

        double distance_[LNNC], Bdistance_[LNNC] = {};

        for (int i = 0; i < LNNC; i++)
            distance_[i] = nodesCoarse_[connec[i]]->getDistFunction();

        for (int icp = 0; icp < LNNC; icp++)
        {

            double wpc[LNNC], phi_[LNNC], xsi[DIM];

            for (int i = 0; i < LNNC; i++)
                wpc[i] = nodesCoarse_[connec[i]]->getWeightPC();
            int *inc_ = nodesCoarse_[connec[LNNC - 1]]->getINC();
            int patch_ = elementsCoarse_[jel]->getPatch();
            for (int i = 0; i < DIM; i++)
                xsi[i] = shapeQuad.ParametricCoordBezier(icp, i);

            shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParCoarse, patch_);

            for (int i = 0; i < LNNC; i++)
                Bdistance_[icp] += phi_[i] * distance_[i];
        };

        for (int ino = 0; ino < LNNC; ino++)
        {
            if ((Bdistance_[ino] <= glueZoneThickness + 0.00001) && (Bdistance_[ino] >= 0.00001))
            {
                flag += 1;
            };
        };

        if (flag != 0)
        {
            elementsGlueZoneCoarse_.push_back(jel);
            elementsCoarse_[jel]->setGlueZone();
        };
    };

    // Defines which coarse nodes are in the gluing zone
    numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    for (int i = 0; i < numElemGlueZoneCoarse; i++)
    {
        int *connec = elementsCoarse_[elementsGlueZoneCoarse_[i]]->getConnectivity();
        for (int ino = 0; ino < LNNC; ino++)
        {
            nodesCZ2[connec[ino]] += 1;
            nodesCoarse_[connec[ino]]->setGlueZone();
        };
    };

    // Compute number of coarse nodes in the gluing zone
    numNodesGlueZoneCoarse = 0;
    for (int i = 0; i < numNodesCoarse; i++)
    {
        if (nodesCZ2[i] > 0)
        {
            numNodesGlueZoneCoarse += 1;
            nodesGlueZoneCoarse_.push_back(i);
        };
    };
};

template <int DIM>
void Arlequin<DIM>::setWeightFunction_FEM_ISO()
{

    double wFuncValue;

    double glueZoneThickness = fineModel.glueZoneThickness;
    double arlequinEpsilon = fineModel.arlequinEpsilon;

    glueZoneThickness *= 1.00; // Thickness from gluing zone

    // IGA coarse mesh
    for (int iNode = 0; iNode < numNodesCoarse; iNode++)
    {

        double r = nodesCoarse_[iNode]->getDistFunction();

        if (r < 0)
        {
            wFuncValue = 1.;
        }
        else
        {
            if (r >= glueZoneThickness)
            {
                wFuncValue = arlequinEpsilon;
            }
            else
            {
                wFuncValue = 1. - (1. - arlequinEpsilon) / glueZoneThickness * r;
                if (wFuncValue < arlequinEpsilon)
                    wFuncValue = arlequinEpsilon;
            };
        };

        nodesCoarse_[iNode]->setWeightFunction(wFuncValue);
    }; // ielem

    for (int jel = 0; jel < numElemCoarse; jel++)
    {
        elementsCoarse_[jel]->setIntegPointWeightFunction_ISO();
    };

    // FEM FINE MESH
    for (int iNode = 0; iNode < numNodesFine; iNode++)
    {

        double r = nodesFine_[iNode]->getDistFunction();

        if (r >= glueZoneThickness)
        {
            wFuncValue = 1. - arlequinEpsilon;
        }
        else
        {
            wFuncValue = (1. - arlequinEpsilon) / glueZoneThickness * r;
            if (wFuncValue > (1. - arlequinEpsilon))
                wFuncValue = 1. - arlequinEpsilon;
        };
        nodesFine_[iNode]->setWeightFunction(wFuncValue);
    }; // ielem

    for (int jel = 0; jel < numElemFine; jel++)
    {
        elementsFine_[jel]->setIntegPointWeightFunction_FEM();
    };

    return;
};

template <int DIM>
void Arlequin<DIM>::setWeightFunction_FEM_FEM()
{

    double wFuncValue;
    double glueZoneThickness = fineModel.glueZoneThickness;
    double arlequinEpsilon = fineModel.arlequinEpsilon;

    glueZoneThickness *= 1.01; // Thickness from gluing zone

    // FEM coarse mesh
    for (int iNode = 0; iNode < numNodesCoarse; iNode++)
    {

        double r = nodesCoarse_[iNode]->getDistFunction();

        if (r < 0)
        {
            wFuncValue = 1.;
        }
        else
        {
            if (r >= glueZoneThickness)
            {
                wFuncValue = arlequinEpsilon;
            }
            else
            {
                wFuncValue = 1. - (1. - arlequinEpsilon) / glueZoneThickness * r;
                if (wFuncValue < arlequinEpsilon)
                    wFuncValue = arlequinEpsilon;
            };
        };

        nodesCoarse_[iNode]->setWeightFunction(wFuncValue);
    }; // ielem

    for (int jel = 0; jel < numElemCoarse; jel++)
    {
        elementsCoarse_[jel]->setIntegPointWeightFunction_FEM();
    };

    // FEM FINE MESH
    for (int iNode = 0; iNode < numNodesFine; iNode++)
    {

        double r = nodesFine_[iNode]->getDistFunction();

        if (r >= glueZoneThickness)
        {
            wFuncValue = 1. - arlequinEpsilon;
        }
        else
        {
            wFuncValue = (1. - arlequinEpsilon) / glueZoneThickness * r;
            if (wFuncValue > (1. - arlequinEpsilon))
                wFuncValue = 1. - arlequinEpsilon;
        };
        nodesFine_[iNode]->setWeightFunction(wFuncValue);
    }; // ielem

    for (int jel = 0; jel < numElemFine; jel++)
    {
        elementsFine_[jel]->setIntegPointWeightFunction_FEM();
    };

    return;
};


template <int DIM>
void Arlequin<DIM>::setWeightFunction_ISO_ISO()
{

    double wFuncValue;

    double glueZoneThickness = fineModel.glueZoneThickness;
    double arlequinEpsilon = fineModel.arlequinEpsilon;

    glueZoneThickness *= 1.00; // Thickness from gluing zone

    // IGA coarse mesh
    for (int iNode = 0; iNode < numNodesCoarse; iNode++)
    {

        double r = nodesCoarse_[iNode]->getDistFunction();

        if (r < 0)
        {
            wFuncValue = 1.;
        }
        else
        {
            if (r >= glueZoneThickness)
            {
                wFuncValue = arlequinEpsilon;
            }
            else
            {
                wFuncValue = 1. - (1. - arlequinEpsilon) / glueZoneThickness * r;
                if (wFuncValue < arlequinEpsilon)
                    wFuncValue = arlequinEpsilon;
            };
        };

        nodesCoarse_[iNode]->setWeightFunction(wFuncValue);
    }; // ielem

    for (int jel = 0; jel < numElemCoarse; jel++)
    {
        elementsCoarse_[jel]->setIntegPointWeightFunction_ISO();
    };

    // IGA FINE MESH
    for (int iNode = 0; iNode < numNodesFine; iNode++)
    {

        double r = nodesFine_[iNode]->getDistFunction();

        if (r >= glueZoneThickness)
        {
            wFuncValue = 1. - arlequinEpsilon;
        }
        else
        {
            wFuncValue = (1. - arlequinEpsilon) / glueZoneThickness * r;
            if (wFuncValue > (1. - arlequinEpsilon))
                wFuncValue = 1. - arlequinEpsilon;
        };
        nodesFine_[iNode]->setWeightFunction(wFuncValue);
    }; // ielem

    for (int jel = 0; jel < numElemFine; jel++)
    {
        elementsFine_[jel]->setIntegPointWeightFunction_ISO();
    };

    return;
};

//------------------------------------------------------------------------------
//------------------------PRINT RESULTS IN PARAVIEW-----------------------------
//------------------------------------------------------------------------------
template <>
void Arlequin<2>::printResults_FEM_ISO(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    const int DIM = 2;

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 18 * DIM - 27;

            int numBezierNodes = coarseModel.NumBezierNodes;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numBezierNodes
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            // Bezier Extraction
            double Vel[numBezierNodes][DIM], realVel[numBezierNodes][DIM], Coord[numBezierNodes][DIM];
            double Press[numBezierNodes], realPress[numBezierNodes], Distance[numBezierNodes], EnergyW[numBezierNodes];

            for (int iElem = 0; iElem < numElemCoarse; iElem++)
            {

                int *Bconnec = elementsCoarse_[iElem]->getBezierConnectivity();
                int *connec = elementsCoarse_[iElem]->getConnectivity();

                double coord_[LNNC][DIM], vel_[LNNC][DIM], realvel_[LNNC][DIM];
                double press_[LNNC], realpress_[LNNC], distance_[LNNC], energyW_[LNNC];

                // Data in the NURBS control points
                for (int i = 0; i < LNNC; i++)
                {

                    double *x = nodesCoarse_[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                    {
                        coord_[i][j] = x[j];
                        vel_[i][j] = nodesCoarse_[connec[i]]->getVelocity(j);
                        realvel_[i][j] = nodesCoarse_[connec[i]]->getVelocityArlequin(j);
                    };

                    press_[i] = nodesCoarse_[connec[i]]->getPressure();
                    realpress_[i] = nodesCoarse_[connec[i]]->getPressureArlequin();
                    distance_[i] = nodesCoarse_[connec[i]]->getDistFunction();
                    energyW_[i] = nodesCoarse_[connec[i]]->getWeightFunction();
                };

                // interpolated values (Bézier variables)
                double Bcoord_[LNNC][DIM] = {};
                double Bvel_[LNNC][DIM] = {};
                double Brealvel_[LNNC][DIM] = {};
                double Bpress_[LNNC] = {};
                double Brealpress_[LNNC] = {};
                double Bdistance_[LNNC] = {};
                double BenergyW_[LNNC] = {};

                for (int i = 0; i < LNNC; i++)
                {

                    QuadShapeFunction<DIM> shapeQuad;
                    double phi_[LNNC], wpc[LNNC], xsi[DIM];
                    for (int k = 0; k < LNNC; k++)
                        wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                    int *inc_ = nodesCoarse_[connec[LNNC - 1]]->getINC();
                    int patch = elementsCoarse_[iElem]->getPatch();
                    for (int j = 0; j < DIM; j++)
                        xsi[j] = shapeQuad.ParametricCoordBezier(i, j);
                    shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParCoarse, patch);

                    for (int j = 0; j < LNNC; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            Bcoord_[i][k] += phi_[j] * coord_[j][k];
                            Bvel_[i][k] += phi_[j] * vel_[j][k];
                            Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bpress_[i] += phi_[j] * press_[j];
                        Brealpress_[i] += phi_[j] * realpress_[j];
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };
                };

                for (int i = 0; i < LNNC; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        Coord[Bconnec[i]][j] = Bcoord_[i][j];
                        Vel[Bconnec[i]][j] = Bvel_[i][j];
                        realVel[Bconnec[i]][j] = Brealvel_[i][j];
                    };
                    Press[Bconnec[i]] = Bpress_[i];
                    realPress[Bconnec[i]] = Brealpress_[i];
                    Distance[Bconnec[i]] = Bdistance_[i];
                    EnergyW[Bconnec[i]] = BenergyW_[i];
                };

            }; // iElem

            for (int i = 0; i < numBezierNodes; i++)
            {
                output_v << Coord[i][0] << " " << Coord[i][1] << " " << 0. << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {

                int *Bconnec_ = elementsCoarse_[iElem]->getBezierConnectivity();
                int Bconnec[LNNC];
                for (int i = 0; i < LNNC; i++)
                    Bconnec[i] = Bconnec_[i];

                output_v << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[8] << " "
                         << Bconnec[6] << " " << Bconnec[1] << " " << Bconnec[5] << " "
                         << Bconnec[7] << " " << Bconnec[3] << " " << Bconnec[4] << std::endl;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + LNNC << std::endl;
                aux += LNNC;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 70 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Vel[i][0] << " "
                             << Vel[i][1] << " "
                             << 0. << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << realVel[i][0] << " "
                             << realVel[i][1] << " "
                             << 0. << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printPressure)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << 0. << " " << 0. << " "
                             << Press[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealPressure)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << 0. << " " << 0. << " "
                             << realPress[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Distance[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << EnergyW[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            std::string f = "FINEoutput" + result + ".vtu";

            std::fstream output_vf(f.c_str(), std::ios_base::out);

            int LNN = 4 * DIM - 2;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                      << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                      << "  <UnstructuredGrid>" << std::endl
                      << "  <Piece NumberOfPoints=\"" << numNodesFine
                      << "\"  NumberOfCells=\"" << numElemFine
                      << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                      << "      <DataArray type=\"Float64\" "
                      << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numNodesFine; i++)
            {
                double *x = nodesFine_[i]->getCoordinates();
                output_vf << x[0] << " " << x[1] << " " << 0.1 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                      << "      <DataArray type=\"Int32\" "
                      << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                int *connec = elementsFine_[i]->getConnectivity();
                int con[LNN];
                for (int i = 0; i < LNN; i++)
                    con[i] = connec[i];
                output_vf << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                      << " Name=\"offsets\" format=\"ascii\">" << std::endl;

            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                      << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 22 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocity(0) << " "
                              << nodesFine_[i]->getVelocity(1) << " "
                              << 0. << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocityArlequin(0) << " "
                              << nodesFine_[i]->getVelocityArlequin(1) << " "
                              << 0. << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printPressure)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << 0. << " " << 0. << " "
                              << nodesFine_[i]->getPressure() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealPressure)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << 0. << " " << 0. << " "
                              << nodesFine_[i]->getPressureArlequin() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printLagrangeMultipliers)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getLagrangeMultiplier(0) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(1) << " "
                              << 0.0 << " " << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getDistFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                      << "Name=\"Normal\" format=\"ascii\">" << std::endl;
            for (int i = 0; i < numNodesFine; i++)
            {

                double *n = nodesFine_[i]->getInnerNormal();
                output_vf << n[0] << " "
                          << n[1] << " "
                          << 0.0 << " " << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;

            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getWeightFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            int cont = 0;
            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont++;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE

            output_vf << "  </Piece>" << std::endl
                      << "  </UnstructuredGrid>" << std::endl
                      << "</VTKFile>" << std::endl;

        }; // rank == 0
    };     // printfreq
};

template <>
void Arlequin<3>::printResults_FEM_ISO(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    const int DIM = 3;

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 18*DIM - 27;

            int numBezierNodes = coarseModel.NumBezierNodes;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numBezierNodes
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            // Bezier Extraction
            double Vel[numBezierNodes][DIM], realVel[numBezierNodes][DIM], Coord[numBezierNodes][DIM];
            double Press[numBezierNodes], realPress[numBezierNodes], Distance[numBezierNodes], EnergyW[numBezierNodes];

            for (int iElem = 0; iElem < numElemCoarse; iElem++)
            {

                int *Bconnec = elementsCoarse_[iElem]->getBezierConnectivity();

                int *connec = elementsCoarse_[iElem]->getConnectivity();

                double coord_[LNNC][DIM], vel_[LNNC][DIM], realvel_[LNNC][DIM];
                double press_[LNNC], realpress_[LNNC], distance_[LNNC], energyW_[LNNC];

                // Data in the NURBS control points
                for (int i = 0; i < LNNC; i++)
                {

                    double *x = nodesCoarse_[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                    {
                        coord_[i][j] = x[j];
                        vel_[i][j] = nodesCoarse_[connec[i]]->getVelocity(j);
                        realvel_[i][j] = nodesCoarse_[connec[i]]->getVelocityArlequin(j);
                    };

                    press_[i] = nodesCoarse_[connec[i]]->getPressure();
                    realpress_[i] = nodesCoarse_[connec[i]]->getPressureArlequin();
                    distance_[i] = nodesCoarse_[connec[i]]->getDistFunction();
                    energyW_[i] = nodesCoarse_[connec[i]]->getWeightFunction();
                };

                // interpolated values (Bézier variables)
                double Bcoord_[LNNC][DIM] = {};
                double Bvel_[LNNC][DIM] = {};
                double Brealvel_[LNNC][DIM] = {};
                double Bpress_[LNNC] = {};
                double Brealpress_[LNNC] = {};
                double Bdistance_[LNNC] = {};
                double BenergyW_[LNNC] = {};

                for (int i = 0; i < LNNC; i++)
                {

                    QuadShapeFunction<DIM> shapeQuad;
                    double phi_[LNNC], wpc[LNNC], xsi[DIM];
                    for (int k = 0; k < LNNC; k++)
                        wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                    int *inc_ = nodesCoarse_[connec[LNNC - 1]]->getINC();
                    int patch = elementsCoarse_[iElem]->getPatch();
                    for (int j = 0; j < DIM; j++)
                        xsi[j] = shapeQuad.ParametricCoordBezier(i, j);
                    shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParCoarse, patch);

                    for (int j = 0; j < LNNC; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            Bcoord_[i][k] += phi_[j] * coord_[j][k];
                            Bvel_[i][k] += phi_[j] * vel_[j][k];
                            Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bpress_[i] += phi_[j] * press_[j];
                        Brealpress_[i] += phi_[j] * realpress_[j];
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };
                };

                for (int i = 0; i < LNNC; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        Coord[Bconnec[i]][j] = Bcoord_[i][j];
                        Vel[Bconnec[i]][j] = Bvel_[i][j];
                        realVel[Bconnec[i]][j] = Brealvel_[i][j];
                    };
                    Press[Bconnec[i]] = Bpress_[i];
                    realPress[Bconnec[i]] = Brealpress_[i];
                    Distance[Bconnec[i]] = Bdistance_[i];
                    EnergyW[Bconnec[i]] = BenergyW_[i];
                };

            }; // iElem

            for (int i = 0; i < numBezierNodes; i++)
            {
                output_v << Coord[i][0] << " " << Coord[i][1] << " " << Coord[i][2] << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {

                int *Bconnec_ = elementsCoarse_[iElem]->getBezierConnectivity();
                int Bconnec[LNNC];
                for (int i = 0; i < LNNC; i++)
                    Bconnec[i] = Bconnec_[i];

                output_v << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[20] << " "
                         << Bconnec[18] << " " << Bconnec[6] << " " << Bconnec[8] << " "
                         << Bconnec[26] << " " << Bconnec[24] << " " << Bconnec[1] << " "
                         << Bconnec[11] << " " << Bconnec[19] << " " << Bconnec[9] << " "
                         << Bconnec[7] << " " << Bconnec[17] << " " << Bconnec[25] << " "
                         << Bconnec[15] << " " << Bconnec[3] << " " << Bconnec[5] << " "
                         << Bconnec[23] << " " << Bconnec[21] << std::endl;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + 20 << std::endl;
                aux += 20;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 25 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Vel[i][0] << " "
                             << Vel[i][1] << " "
                             << Vel[i][2] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << realVel[i][0] << " "
                             << realVel[i][1] << " "
                             << realVel[i][2] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printPressure)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << 0. << " " << 0. << " "
                             << Press[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealPressure)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << 0. << " " << 0. << " "
                             << realPress[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Distance[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << EnergyW[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            std::string f = "FINEoutput" + result + ".vtu";

            std::fstream output_vf(f.c_str(), std::ios_base::out);

            int LNN = 4*DIM-2;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                      << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                      << "  <UnstructuredGrid>" << std::endl
                      << "  <Piece NumberOfPoints=\"" << numNodesFine
                      << "\"  NumberOfCells=\"" << numElemFine
                      << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                      << "      <DataArray type=\"Float64\" "
                      << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numNodesFine; i++)
            {
                double *x = nodesFine_[i]->getCoordinates();
                output_vf << x[0] << " " << x[1] << " " << x[2]+0.1 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                      << "      <DataArray type=\"Int32\" "
                      << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                int *connec = elementsFine_[i]->getConnectivity();
                int con[LNN];
                for (int i = 0; i < LNN; i++)
                    con[i] = connec[i];
                output_vf << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << " "
                          << con[6] << " " << con[7] << " " << con[9] << " "
                          << con[8] << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                      << " Name=\"offsets\" format=\"ascii\">" << std::endl;

            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                      << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 24 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocity(0) << " "
                              << nodesFine_[i]->getVelocity(1) << " "
                              << nodesFine_[i]->getVelocity(2) << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocityArlequin(0) << " "
                              << nodesFine_[i]->getVelocityArlequin(1) << " "
                              << nodesFine_[i]->getVelocityArlequin(2) << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printPressure)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << 0. << " " << 0. << " "
                              << nodesFine_[i]->getPressure() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealPressure)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << 0. << " " << 0. << " "
                              << nodesFine_[i]->getPressureArlequin() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printLagrangeMultipliers)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getLagrangeMultiplier(0) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(1) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(2) << " " << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getDistFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            // output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            //           << "Name=\"Normal\" format=\"ascii\">" << std::endl;
            // for (int i = 0; i < numNodesFine; i++)
            // {

            //     double *n = nodesFine_[i]->getInnerNormal();
            //     output_vf << n[0] << " "
            //               << n[1] << " "
            //               << n[2] << " " << std::endl;
            // };
            // output_vf << "      </DataArray> " << std::endl;

            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getWeightFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            int cont = 0;
            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont++;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE

            output_vf << "  </Piece>" << std::endl
                      << "  </UnstructuredGrid>" << std::endl
                      << "</VTKFile>" << std::endl;

        }; // rank == 0
    };     // printfreq
};


template <>
void Arlequin<2>::printResults_FEM_FEM(int step)
{
    int DIM = 2;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 4*DIM-2;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numNodesCoarse
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;


            for (int i = 0; i < numNodesCoarse; i++)
            {
                double *x = nodesCoarse_[i]->getCoordinates();
                output_v << x[0] << " " << x[1] << " " << 0.0 << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {
                int *connec = elementsCoarse_[iElem]->getConnectivity();
                int con[LNNC];
                for (int i = 0; i < LNNC; i++) con[i] = connec[i];
                    output_v << con[0] << " " << con[1] << " " << con[2] << " "
                             << con[3] << " " << con[4] << " " << con[5] << std::endl;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + LNNC << std::endl;
                aux += LNNC;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 22 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getVelocity(0) << " "
                             << nodesCoarse_[i]->getVelocity(1) << " "
                             << 0.0 << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getVelocityArlequin(0) << " "
                             << nodesCoarse_[i]->getVelocityArlequin(1) << " "
                             << 0.0 << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printPressure)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << 0. << " " << 0. << " "
                             << nodesCoarse_[i]->getPressure() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealPressure)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << 0. << " " << 0. << " "
                             << nodesCoarse_[i]->getPressureArlequin() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getDistFunction() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getWeightFunction() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            std::string f = "FINEoutput" + result + ".vtu";

            std::fstream output_vf(f.c_str(), std::ios_base::out);

            int LNN = 4*DIM-2;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                      << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                      << "  <UnstructuredGrid>" << std::endl
                      << "  <Piece NumberOfPoints=\"" << numNodesFine
                      << "\"  NumberOfCells=\"" << numElemFine
                      << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                      << "      <DataArray type=\"Float64\" "
                      << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numNodesFine; i++)
            {
                double *x = nodesFine_[i]->getCoordinates();
                output_vf << x[0] << " " << x[1] << " " << 0.0 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                      << "      <DataArray type=\"Int32\" "
                      << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                int *connec = elementsFine_[i]->getConnectivity();
                int con[LNN];
                for (int i = 0; i < LNN; i++) con[i] = connec[i];
                     
                output_vf << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << std::endl;

            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                      << " Name=\"offsets\" format=\"ascii\">" << std::endl;

            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                      << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 22 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocity(0) << " "
                              << nodesFine_[i]->getVelocity(1) << " "
                              << 0.0 << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocityArlequin(0) << " "
                              << nodesFine_[i]->getVelocityArlequin(1) << " "
                              << 0.0 << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printPressure)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << 0. << " " << 0. << " "
                              << nodesFine_[i]->getPressure() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealPressure)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << 0. << " " << 0. << " "
                              << nodesFine_[i]->getPressureArlequin() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printLagrangeMultipliers)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getLagrangeMultiplier(0) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(1) << " "
                              << 0.0 << " " << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getDistFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };


            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getWeightFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            int cont = 0;
            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont++;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE

            output_vf << "  </Piece>" << std::endl
                      << "  </UnstructuredGrid>" << std::endl
                      << "</VTKFile>" << std::endl;

        }; // rank == 0
    };     // printfreq
};



template <>
void Arlequin<3>::printResults_FEM_FEM(int step)
{
    int DIM = 3;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 4*DIM-2;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numNodesCoarse
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;


            for (int i = 0; i < numNodesCoarse; i++)
            {
                double *x = nodesCoarse_[i]->getCoordinates();
                output_v << x[0] << " " << x[1] << " " << x[2] << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {
                int *connec = elementsCoarse_[iElem]->getConnectivity();
                int con[LNNC];
                for (int i = 0; i < LNNC; i++) con[i] = connec[i];
                output_v << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << " "
                          << con[6] << " " << con[7] << " " << con[9] << " "
                          << con[8] << std::endl;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + LNNC << std::endl;
                aux += LNNC;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 24 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getVelocity(0) << " "
                             << nodesCoarse_[i]->getVelocity(1) << " "
                             << nodesCoarse_[i]->getVelocity(2) << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getVelocityArlequin(0) << " "
                             << nodesCoarse_[i]->getVelocityArlequin(1) << " "
                             << nodesCoarse_[i]->getVelocityArlequin(2) << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printPressure)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << 0. << " " << 0. << " "
                             << nodesCoarse_[i]->getPressure() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealPressure)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << 0. << " " << 0. << " "
                             << nodesCoarse_[i]->getPressureArlequin() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getDistFunction() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getWeightFunction() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            std::string f = "FINEoutput" + result + ".vtu";

            std::fstream output_vf(f.c_str(), std::ios_base::out);

            int LNN = 4*DIM-2;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                      << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                      << "  <UnstructuredGrid>" << std::endl
                      << "  <Piece NumberOfPoints=\"" << numNodesFine
                      << "\"  NumberOfCells=\"" << numElemFine
                      << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                      << "      <DataArray type=\"Float64\" "
                      << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numNodesFine; i++)
            {
                double *x = nodesFine_[i]->getCoordinates();
                output_vf << x[0] << " " << x[1] << " " << x[2]+0.1<< std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                      << "      <DataArray type=\"Int32\" "
                      << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                int *connec = elementsFine_[i]->getConnectivity();
                int con[LNN];
                for (int i = 0; i < LNN; i++)
                    con[i] = connec[i];
                output_vf << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << " "
                          << con[6] << " " << con[7] << " " << con[9] << " "
                          << con[8] << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                      << " Name=\"offsets\" format=\"ascii\">" << std::endl;

            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                      << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 24 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocity(0) << " "
                              << nodesFine_[i]->getVelocity(1) << " "
                              << nodesFine_[i]->getVelocity(2) << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocityArlequin(0) << " "
                              << nodesFine_[i]->getVelocityArlequin(1) << " "
                              << nodesFine_[i]->getVelocityArlequin(2) << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printPressure)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << 0. << " " << 0. << " "
                              << nodesFine_[i]->getPressure() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealPressure)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << 0. << " " << 0. << " "
                              << nodesFine_[i]->getPressureArlequin() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printLagrangeMultipliers)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getLagrangeMultiplier(0) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(1) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(2) << " " << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getDistFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            // output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            //           << "Name=\"Normal\" format=\"ascii\">" << std::endl;
            // for (int i = 0; i < numNodesFine; i++)
            // {

            //     double *n = nodesFine_[i]->getInnerNormal();
            //     output_vf << n[0] << " "
            //               << n[1] << " "
            //               << n[2] << " " << std::endl;
            // };
            // output_vf << "      </DataArray> " << std::endl;

            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getWeightFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            int cont = 0;
            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont++;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE

            output_vf << "  </Piece>" << std::endl
                      << "  </UnstructuredGrid>" << std::endl
                      << "</VTKFile>" << std::endl;

        }; // rank == 0
    };     // printfreq
};

template <>
void Arlequin<2>::printResultsLaplace_FEM_ISO(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    const int DIM = 2;

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 18*DIM-27;

            int numBezierNodes = coarseModel.NumBezierNodes;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numBezierNodes
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            // Bezier Extraction
            double Vel[numBezierNodes][DIM], realVel[numBezierNodes][DIM], Coord[numBezierNodes][DIM];
            double Distance[numBezierNodes], EnergyW[numBezierNodes];

            for (int iElem = 0; iElem < numElemCoarse; iElem++)
            {

                int *Bconnec = elementsCoarse_[iElem]->getBezierConnectivity();

                int *connec = elementsCoarse_[iElem]->getConnectivity();

                double coord_[LNNC][DIM], vel_[LNNC][DIM], realvel_[LNNC][DIM];
                double distance_[LNNC], energyW_[LNNC];

                // Data in the NURBS control points
                for (int i = 0; i < LNNC; i++)
                {

                    double *x = nodesCoarse_[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                    {
                        coord_[i][j] = x[j];
                        vel_[i][j] = nodesCoarse_[connec[i]]->getVelocity(j);
                        realvel_[i][j] = nodesCoarse_[connec[i]]->getVelocityArlequin(j);
                    };

                    distance_[i] = nodesCoarse_[connec[i]]->getDistFunction();
                    energyW_[i] = nodesCoarse_[connec[i]]->getWeightFunction();
                };

                // interpolated values (Bézier variables)
                double Bcoord_[LNNC][DIM] = {};
                double Bvel_[LNNC][DIM] = {};
                double Brealvel_[LNNC][DIM] = {};
                double Bdistance_[LNNC] = {};
                double BenergyW_[LNNC] = {};

                for (int i = 0; i < LNNC; i++)
                {

                    QuadShapeFunction<DIM> shapeQuad;
                    double phi_[LNNC], wpc[LNNC], xsi[DIM];
                    for (int k = 0; k < LNNC; k++)
                        wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                    int *inc_ = nodesCoarse_[connec[LNNC - 1]]->getINC();
                    int patch = elementsCoarse_[iElem]->getPatch();
                    for (int j = 0; j < DIM; j++)
                        xsi[j] = shapeQuad.ParametricCoordBezier(i, j);
                    shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParCoarse, patch);

                    for (int j = 0; j < LNNC; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            Bcoord_[i][k] += phi_[j] * coord_[j][k];
                            Bvel_[i][k] += phi_[j] * vel_[j][k];
                            Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };
                };

                for (int i = 0; i < LNNC; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        Coord[Bconnec[i]][j] = Bcoord_[i][j];
                        Vel[Bconnec[i]][j] = Bvel_[i][j];
                        realVel[Bconnec[i]][j] = Brealvel_[i][j];
                    };
                    Distance[Bconnec[i]] = Bdistance_[i];
                    EnergyW[Bconnec[i]] = BenergyW_[i];
                };

            }; // iElem

            for (int i = 0; i < numBezierNodes; i++)
            {
                output_v << Coord[i][0] << " " << Coord[i][1] << " " << 0.0 << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {

                int *Bconnec_ = elementsCoarse_[iElem]->getBezierConnectivity();
                int Bconnec[LNNC];
                for (int i = 0; i < LNNC; i++)
                    Bconnec[i] = Bconnec_[i];

                output_v << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[8] << " "
                         << Bconnec[6] << " " << Bconnec[1] << " " << Bconnec[5] << " "
                         << Bconnec[7] << " " << Bconnec[3] << " " << Bconnec[4] << std::endl;
            
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + LNNC << std::endl;
                aux += LNNC;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 70 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Vel[i][0] << " "
                             << Vel[i][1] << " "
                             << 0.0 << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << realVel[i][0] << " "
                             << realVel[i][1] << " "
                             << 0.0 << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };


            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Distance[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << EnergyW[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            std::string f = "FINEoutput" + result + ".vtu";

            std::fstream output_vf(f.c_str(), std::ios_base::out);

            int LNN = 4*DIM-2;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                      << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                      << "  <UnstructuredGrid>" << std::endl
                      << "  <Piece NumberOfPoints=\"" << numNodesFine
                      << "\"  NumberOfCells=\"" << numElemFine
                      << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                      << "      <DataArray type=\"Float64\" "
                      << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numNodesFine; i++)
            {
                double *x = nodesFine_[i]->getCoordinates();
                output_vf << x[0] << " " << x[1] << " " << 0.1 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                      << "      <DataArray type=\"Int32\" "
                      << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                int *connec = elementsFine_[i]->getConnectivity();
                int con[LNN];
                for (int i = 0; i < LNN; i++)
                    con[i] = connec[i];
                output_vf << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                      << " Name=\"offsets\" format=\"ascii\">" << std::endl;

            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                      << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 22 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocity(0) << " "
                              << nodesFine_[i]->getVelocity(1) << " "
                              << 0.0 << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocityArlequin(0) << " "
                              << nodesFine_[i]->getVelocityArlequin(1) << " "
                              << 0.0 << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };


            if (fineModel.printLagrangeMultipliers)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getLagrangeMultiplier(0) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(1) << " "
                              << 0.0 << " " << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getDistFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            // output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            //           << "Name=\"Normal\" format=\"ascii\">" << std::endl;
            // for (int i = 0; i < numNodesFine; i++)
            // {

            //     double *n = nodesFine_[i]->getInnerNormal();
            //     output_vf << n[0] << " "
            //               << n[1] << " "
            //               << 0.0 << " " << std::endl;
            // };
            // output_vf << "      </DataArray> " << std::endl;

            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getWeightFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            int cont = 0;
            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont++;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE

            output_vf << "  </Piece>" << std::endl
                      << "  </UnstructuredGrid>" << std::endl
                      << "</VTKFile>" << std::endl;

        }; // rank == 0
    };     // printfreq
};

template <>
void Arlequin<3>::printResultsLaplace_FEM_ISO(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int const DIM = 3;

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 18*DIM-27;

            int numBezierNodes = coarseModel.NumBezierNodes;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numBezierNodes
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            // Bezier Extraction
            double Vel[numBezierNodes][DIM], realVel[numBezierNodes][DIM], Coord[numBezierNodes][DIM];
            double Distance[numBezierNodes], EnergyW[numBezierNodes];

            for (int iElem = 0; iElem < numElemCoarse; iElem++)
            {

                int *Bconnec = elementsCoarse_[iElem]->getBezierConnectivity();

                int *connec = elementsCoarse_[iElem]->getConnectivity();

                double coord_[LNNC][DIM], vel_[LNNC][DIM], realvel_[LNNC][DIM];
                double distance_[LNNC], energyW_[LNNC];

                // Data in the NURBS control points
                for (int i = 0; i < LNNC; i++)
                {

                    double *x = nodesCoarse_[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                    {
                        coord_[i][j] = x[j];
                        vel_[i][j] = nodesCoarse_[connec[i]]->getVelocity(j);
                        realvel_[i][j] = nodesCoarse_[connec[i]]->getVelocityArlequin(j);
                    };

                    distance_[i] = nodesCoarse_[connec[i]]->getDistFunction();
                    energyW_[i] = nodesCoarse_[connec[i]]->getWeightFunction();
                };

                // interpolated values (Bézier variables)
                double Bcoord_[LNNC][DIM] = {};
                double Bvel_[LNNC][DIM] = {};
                double Brealvel_[LNNC][DIM] = {};
                double Bdistance_[LNNC] = {};
                double BenergyW_[LNNC] = {};

                for (int i = 0; i < LNNC; i++)
                {

                    QuadShapeFunction<DIM> shapeQuad;
                    double phi_[LNNC], wpc[LNNC], xsi[DIM];
                    for (int k = 0; k < LNNC; k++)
                        wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                    int *inc_ = nodesCoarse_[connec[LNNC - 1]]->getINC();
                    int patch = elementsCoarse_[iElem]->getPatch();
                    for (int j = 0; j < DIM; j++)
                        xsi[j] = shapeQuad.ParametricCoordBezier(i, j);
                    shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParCoarse, patch);

                    for (int j = 0; j < LNNC; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            Bcoord_[i][k] += phi_[j] * coord_[j][k];
                            Bvel_[i][k] += phi_[j] * vel_[j][k];
                            Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };
                };

                for (int i = 0; i < LNNC; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        Coord[Bconnec[i]][j] = Bcoord_[i][j];
                        Vel[Bconnec[i]][j] = Bvel_[i][j];
                        realVel[Bconnec[i]][j] = Brealvel_[i][j];
                    };
                    Distance[Bconnec[i]] = Bdistance_[i];
                    EnergyW[Bconnec[i]] = BenergyW_[i];
                };

            }; // iElem

            for (int i = 0; i < numBezierNodes; i++)
            {
                output_v << Coord[i][0] << " " << Coord[i][1] << " " << Coord[i][2] << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {

                int *Bconnec_ = elementsCoarse_[iElem]->getBezierConnectivity();
                int Bconnec[LNNC];
                for (int i = 0; i < LNNC; i++)
                    Bconnec[i] = Bconnec_[i];

                output_v << Bconnec[18] << " " << Bconnec[20] << " " <<  Bconnec[26] << " " << Bconnec[24] << " " << Bconnec[0] << " " <<
                            Bconnec[2]  << " " << Bconnec[8]  << " " <<  Bconnec[6]  << " " << Bconnec[19] << " " << Bconnec[23] << " " <<
                            Bconnec[25] << " " << Bconnec[21] << " " <<  Bconnec[1]  << " " << Bconnec[5]  << " " << Bconnec[7] << " " << 
                            Bconnec[3] << " " << Bconnec[9]  << " " <<  Bconnec[11] << " " << Bconnec[15] << " " << Bconnec[17] << " " <<
                            Bconnec[12] << " " << Bconnec[14] << " " <<  Bconnec[10] << " " << Bconnec[16] << " " << Bconnec[22] << " " <<
                            Bconnec[4]  << " " << Bconnec[13] << std::endl;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + LNNC << std::endl;
                aux += LNNC;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 72 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Vel[i][0] << " "
                             << Vel[i][1] << " "
                             << Vel[i][2] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << realVel[i][0] << " "
                             << realVel[i][1] << " "
                             << realVel[i][2] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };


            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Distance[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << EnergyW[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            std::string f = "FINEoutput" + result + ".vtu";

            std::fstream output_vf(f.c_str(), std::ios_base::out);

            int LNN = 4*DIM-2;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                      << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                      << "  <UnstructuredGrid>" << std::endl
                      << "  <Piece NumberOfPoints=\"" << numNodesFine
                      << "\"  NumberOfCells=\"" << numElemFine
                      << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                      << "      <DataArray type=\"Float64\" "
                      << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numNodesFine; i++)
            {
                double *x = nodesFine_[i]->getCoordinates();
                output_vf << x[0] << " " << x[1] << " " << x[2]+0.2 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                      << "      <DataArray type=\"Int32\" "
                      << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                int *connec = elementsFine_[i]->getConnectivity();
                int con[LNN];
                for (int i = 0; i < LNN; i++)
                    con[i] = connec[i];
                output_vf << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << " "
                          << con[6] << " " << con[7] << " " << con[9] << " "
                          << con[8] << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                      << " Name=\"offsets\" format=\"ascii\">" << std::endl;

            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                      << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 24 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocity(0) << " "
                              << nodesFine_[i]->getVelocity(1) << " "
                              << nodesFine_[i]->getVelocity(2) << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocityArlequin(0) << " "
                              << nodesFine_[i]->getVelocityArlequin(1) << " "
                              << nodesFine_[i]->getVelocityArlequin(2) << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };


            if (fineModel.printLagrangeMultipliers)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getLagrangeMultiplier(0) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(1) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(2) << " " << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getDistFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            // output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            //           << "Name=\"Normal\" format=\"ascii\">" << std::endl;
            // for (int i = 0; i < numNodesFine; i++)
            // {

            //     double *n = nodesFine_[i]->getInnerNormal();
            //     output_vf << n[0] << " "
            //               << n[1] << " "
            //               << 0.0 << " " << std::endl;
            // };
            // output_vf << "      </DataArray> " << std::endl;

            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getWeightFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            int cont = 0;
            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont++;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE

            output_vf << "  </Piece>" << std::endl
                      << "  </UnstructuredGrid>" << std::endl
                      << "</VTKFile>" << std::endl;

        }; // rank == 0
    };     // printfreq
};


template <int DIM>
void Arlequin<DIM>::printResultsLaplace_FEM_FEM(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 4*DIM-2;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numNodesCoarse
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;


            for (int i = 0; i < numNodesCoarse; i++)
            {
                double *x = nodesCoarse_[i] -> getCoordinates();
                output_v << x[0] << " " << x[1] << " " << x[2] << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {
                int *connec = elementsCoarse_[iElem]->getConnectivity();
                int con[LNNC];
                for (int i = 0; i < LNNC; i++) con[i] = connec[i];

                output_v << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << " "
                          << con[6] << " " << con[7] << " " << con[9] << " "
                          << con[8] << std::endl;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + LNNC << std::endl;
                aux += LNNC;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 24 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]-> getVelocity(0) << " "
                             << nodesCoarse_[i]-> getVelocity(1) << " "
                             << nodesCoarse_[i]-> getVelocity(2) << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getVelocityArlequin(0) << " "
                              << nodesCoarse_[i]->getVelocityArlequin(1) << " "
                              << nodesCoarse_[i]->getVelocityArlequin(2) << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };


            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getDistFunction() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesCoarse; i++)
                {
                    output_v << nodesCoarse_[i]->getWeightFunction() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            std::string f = "FINEoutput" + result + ".vtu";

            std::fstream output_vf(f.c_str(), std::ios_base::out);

            int LNN = 4*DIM-2;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                      << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                      << "  <UnstructuredGrid>" << std::endl
                      << "  <Piece NumberOfPoints=\"" << numNodesFine
                      << "\"  NumberOfCells=\"" << numElemFine
                      << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                      << "      <DataArray type=\"Float64\" "
                      << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numNodesFine; i++)
            {
                double *x = nodesFine_[i]->getCoordinates();
                output_vf << x[0] << " " << x[1] << " " << x[2]+0.2 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                      << "      <DataArray type=\"Int32\" "
                      << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                int *connec = elementsFine_[i]->getConnectivity();
                int con[LNN];
                for (int i = 0; i < LNN; i++)
                    con[i] = connec[i];
                output_vf << con[0] << " " << con[1] << " " << con[2] << " "
                          << con[3] << " " << con[4] << " " << con[5] << " "
                          << con[6] << " " << con[7] << " " << con[9] << " "
                          << con[8] << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                      << " Name=\"offsets\" format=\"ascii\">" << std::endl;

            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                      << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 24 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                      << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocity(0) << " "
                              << nodesFine_[i]->getVelocity(1) << " "
                              << nodesFine_[i]->getVelocity(2) << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getVelocityArlequin(0) << " "
                              << nodesFine_[i]->getVelocityArlequin(1) << " "
                              << nodesFine_[i]->getVelocityArlequin(2) << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };


            if (fineModel.printLagrangeMultipliers)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                          << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getLagrangeMultiplier(0) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(1) << " "
                              << nodesFine_[i]->getLagrangeMultiplier(2) << " " << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getDistFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            // output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            //           << "Name=\"Normal\" format=\"ascii\">" << std::endl;
            // for (int i = 0; i < numNodesFine; i++)
            // {

            //     double *n = nodesFine_[i]->getInnerNormal();
            //     output_vf << n[0] << " "
            //               << n[1] << " "
            //               << 0.0 << " " << std::endl;
            // };
            // output_vf << "      </DataArray> " << std::endl;

            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numNodesFine; i++)
                {
                    output_vf << nodesFine_[i]->getWeightFunction() << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            int cont = 0;
            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                          << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont++;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE

            output_vf << "  </Piece>" << std::endl
                      << "  </UnstructuredGrid>" << std::endl
                      << "</VTKFile>" << std::endl;

        }; // rank == 0
    };     // printfreq
};

template <>
void Arlequin<2>::printResultsLaplace_ISO_ISO(int step)
{
    int const DIM = 2;
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 18*DIM-27;

            int numBezierNodes = coarseModel.NumBezierNodes;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numBezierNodes
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            // Bezier Extraction
            double Vel[numBezierNodes][DIM], realVel[numBezierNodes][DIM], Coord[numBezierNodes][DIM];
            double Distance[numBezierNodes], EnergyW[numBezierNodes];

            for (int iElem = 0; iElem < numElemCoarse; iElem++)
            {

                int *Bconnec = elementsCoarse_[iElem]->getBezierConnectivity();

                int *connec = elementsCoarse_[iElem]->getConnectivity();

                double coord_[LNNC][DIM], vel_[LNNC][DIM], realvel_[LNNC][DIM];
                double distance_[LNNC], energyW_[LNNC];

                // Data in the NURBS control points
                for (int i = 0; i < LNNC; i++)
                {

                    double *x = nodesCoarse_[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                    {
                        coord_[i][j] = x[j];
                        vel_[i][j] = nodesCoarse_[connec[i]]->getVelocity(j);
                        realvel_[i][j] = nodesCoarse_[connec[i]]->getVelocityArlequin(j);
                    };

                    distance_[i] = nodesCoarse_[connec[i]]->getDistFunction();
                    energyW_[i] = nodesCoarse_[connec[i]]->getWeightFunction();
                };

                // interpolated values (Bézier variables)
                double Bcoord_[LNNC][DIM] = {};
                double Bvel_[LNNC][DIM] = {};
                double Brealvel_[LNNC][DIM] = {};
                double Bdistance_[LNNC] = {};
                double BenergyW_[LNNC] = {};

                for (int i = 0; i < LNNC; i++)
                {

                    QuadShapeFunction<DIM> shapeQuad;
                    double phi_[LNNC], wpc[LNNC], xsi[DIM];
                    for (int k = 0; k < LNNC; k++)
                        wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                    int *inc_ = nodesCoarse_[connec[LNNC - 1]]->getINC();
                    int patch = elementsCoarse_[iElem]->getPatch();
                    for (int j = 0; j < DIM; j++)
                        xsi[j] = shapeQuad.ParametricCoordBezier(i, j);
                    shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParCoarse, patch);

                    for (int j = 0; j < LNNC; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            Bcoord_[i][k] += phi_[j] * coord_[j][k];
                            Bvel_[i][k] += phi_[j] * vel_[j][k];
                            Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };
                };

                for (int i = 0; i < LNNC; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        Coord[Bconnec[i]][j] = Bcoord_[i][j];
                        Vel[Bconnec[i]][j] = Bvel_[i][j];
                        realVel[Bconnec[i]][j] = Brealvel_[i][j];
                    };
                    Distance[Bconnec[i]] = Bdistance_[i];
                    EnergyW[Bconnec[i]] = BenergyW_[i];
                };

            }; // iElem

            for (int i = 0; i < numBezierNodes; i++)
            {
                output_v << Coord[i][0] << " " << Coord[i][1] << " " << 0.0 << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {

                int *Bconnec_ = elementsCoarse_[iElem]->getBezierConnectivity();
                int Bconnec[LNNC];
                for (int i = 0; i < LNNC; i++)
                    Bconnec[i] = Bconnec_[i];

                 output_v << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[8] << " "
                    << Bconnec[6] << " " << Bconnec[1] << " " << Bconnec[5] << " " 
                    << Bconnec[7] << " " << Bconnec[3] << " " << Bconnec[4]<<  std::endl;

            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + LNNC << std::endl;
                aux += LNNC;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 70 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Vel[i][0] << " "
                             << Vel[i][1] << " "
                             << 0.0 << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << realVel[i][0] << " "
                             << realVel[i][1] << " "
                             << 0.0 << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };


            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Distance[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << EnergyW[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            
            // FINE MESH
            s = "FINEoutput" + result + ".vtu";
            std::fstream output_vf(s.c_str(), std::ios_base::out);

            int LNN = 18*DIM-27;

            numBezierNodes = fineModel.NumBezierNodes;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numBezierNodes
                     << "\"  NumberOfCells=\"" << numElemFine
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            // Bezier Extraction
            double VelF[numBezierNodes][DIM], realVelF[numBezierNodes][DIM], CoordF[numBezierNodes][DIM];
            double DistanceF[numBezierNodes], EnergyWF[numBezierNodes];

            for (int iElem = 0; iElem < numElemFine; iElem++)
            {

                int *Bconnec = elementsFine_[iElem]->getBezierConnectivity();

                int *connec = elementsFine_[iElem]->getConnectivity();

                double coord_[LNN][DIM], vel_[LNN][DIM], realvel_[LNN][DIM];
                double distance_[LNN], energyW_[LNN];

                // Data in the NURBS control points
                for (int i = 0; i < LNN; i++)
                {

                    double *x = nodesFine_[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                    {
                        coord_[i][j] = x[j];
                        vel_[i][j] = nodesFine_[connec[i]]->getVelocity(j);
                        realvel_[i][j] = nodesFine_[connec[i]]->getVelocityArlequin(j);
                    };

                    distance_[i] = nodesFine_[connec[i]]->getDistFunction();
                    energyW_[i] = nodesFine_[connec[i]]->getWeightFunction();
                };

                // interpolated values (Bézier variables)
                double Bcoord_[LNN][DIM] = {};
                double Bvel_[LNN][DIM] = {};
                double Brealvel_[LNN][DIM] = {};
                double Bdistance_[LNN] = {};
                double BenergyW_[LNN] = {};

                for (int i = 0; i < LNN; i++)
                {

                    QuadShapeFunction<DIM> shapeQuad;
                    double phi_[LNN], wpc[LNN], xsi[DIM];
                    for (int k = 0; k < LNN; k++)
                        wpc[k] = nodesFine_[connec[k]]->getWeightPC();
                    int *inc_ = nodesFine_[connec[LNN - 1]]->getINC();
                    int patch = elementsFine_[iElem]->getPatch();
                    for (int j = 0; j < DIM; j++)
                        xsi[j] = shapeQuad.ParametricCoordBezier(i, j);
                    shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParFine, patch);

                    for (int j = 0; j < LNN; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            Bcoord_[i][k] += phi_[j] * coord_[j][k];
                            Bvel_[i][k] += phi_[j] * vel_[j][k];
                            Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };
                };

                for (int i = 0; i < LNN; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        CoordF[Bconnec[i]][j] = Bcoord_[i][j];
                        VelF[Bconnec[i]][j] = Bvel_[i][j];
                        realVelF[Bconnec[i]][j] = Brealvel_[i][j];
                    };
                    DistanceF[Bconnec[i]] = Bdistance_[i];
                    EnergyWF[Bconnec[i]] = BenergyW_[i];
                };

            }; // iElem

            for (int i = 0; i < numBezierNodes; i++)
            {
                output_vf << CoordF[i][0] << " " << CoordF[i][1] << " " << 0.0 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemFine; ++iElem)
            {

                int *Bconnec_ = elementsFine_[iElem]->getBezierConnectivity();
                int Bconnec[LNN];
                for (int i = 0; i < LNN; i++)
                    Bconnec[i] = Bconnec_[i];

                output_vf << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[8] << " "
                    << Bconnec[6] << " " << Bconnec[1] << " " << Bconnec[5] << " " 
                    << Bconnec[7] << " " << Bconnec[3] << " " << Bconnec[4]<<  std::endl;

            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 70 << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_vf << VelF[i][0] << " "
                             << VelF[i][1] << " "
                             << 0.0 << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_vf << realVelF[i][0] << " "
                             << realVelF[i][1] << " "
                             << 0.0 << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };


            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_vf << DistanceF[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_vf << EnergyWF[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_vf << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;


        }; // rank == 0
    };     // printfreq
};

template <>
void Arlequin<3>::printResultsLaplace_ISO_ISO(int step)
{
    int const DIM = 3;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (step % fineModel.printFreq == 0)
    {

        if (rank == 0)
        {

            std::string result;
            std::ostringstream convert;
            convert << step + 100000;
            result = convert.str();

            // COARSE MESH
            std::string s = "COARSEoutput" + result + ".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            int LNNC = 18*DIM-27;

            int numBezierNodes = coarseModel.NumBezierNodes;

            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numBezierNodes
                     << "\"  NumberOfCells=\"" << numElemCoarse
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            // Bezier Extraction
            double Vel[numBezierNodes][DIM], realVel[numBezierNodes][DIM], Coord[numBezierNodes][DIM];
            double Distance[numBezierNodes], EnergyW[numBezierNodes];

            for (int iElem = 0; iElem < numElemCoarse; iElem++)
            {

                int *Bconnec = elementsCoarse_[iElem]->getBezierConnectivity();

                int *connec = elementsCoarse_[iElem]->getConnectivity();

                double coord_[LNNC][DIM], vel_[LNNC][DIM], realvel_[LNNC][DIM];
                double distance_[LNNC], energyW_[LNNC];

                // Data in the NURBS control points
                for (int i = 0; i < LNNC; i++)
                {

                    double *x = nodesCoarse_[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                    {
                        coord_[i][j] = x[j];
                        vel_[i][j] = nodesCoarse_[connec[i]]->getVelocity(j);
                        realvel_[i][j] = nodesCoarse_[connec[i]]->getVelocityArlequin(j);
                    };

                    distance_[i] = nodesCoarse_[connec[i]]->getDistFunction();
                    energyW_[i] = nodesCoarse_[connec[i]]->getWeightFunction();
                };

                // interpolated values (Bézier variables)
                double Bcoord_[LNNC][DIM] = {};
                double Bvel_[LNNC][DIM] = {};
                double Brealvel_[LNNC][DIM] = {};
                double Bdistance_[LNNC] = {};
                double BenergyW_[LNNC] = {};

                for (int i = 0; i < LNNC; i++)
                {

                    QuadShapeFunction<DIM> shapeQuad;
                    double phi_[LNNC], wpc[LNNC], xsi[DIM];
                    for (int k = 0; k < LNNC; k++)
                        wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                    int *inc_ = nodesCoarse_[connec[LNNC - 1]]->getINC();
                    int patch = elementsCoarse_[iElem]->getPatch();
                    for (int j = 0; j < DIM; j++)
                        xsi[j] = shapeQuad.ParametricCoordBezier(i, j);
                    shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParCoarse, patch);

                    for (int j = 0; j < LNNC; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            Bcoord_[i][k] += phi_[j] * coord_[j][k];
                            Bvel_[i][k] += phi_[j] * vel_[j][k];
                            Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };
                };

                for (int i = 0; i < LNNC; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        Coord[Bconnec[i]][j] = Bcoord_[i][j];
                        Vel[Bconnec[i]][j] = Bvel_[i][j];
                        realVel[Bconnec[i]][j] = Brealvel_[i][j];
                    };
                    Distance[Bconnec[i]] = Bdistance_[i];
                    EnergyW[Bconnec[i]] = BenergyW_[i];
                };

            }; // iElem

            for (int i = 0; i < numBezierNodes; i++)
            {
                output_v << Coord[i][0] << " " << Coord[i][1] << " " << Coord[i][2] << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemCoarse; ++iElem)
            {

                int *Bconnec_ = elementsCoarse_[iElem]->getBezierConnectivity();
                int Bconnec[LNNC];
                for (int i = 0; i < LNNC; i++)
                    Bconnec[i] = Bconnec_[i];

                // output_v << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[20] << " "
                //          << Bconnec[18] << " " << Bconnec[6] << " " << Bconnec[8] << " "
                //          << Bconnec[26] << " " << Bconnec[24] << " " << Bconnec[1] << " "
                //          << Bconnec[11] << " " << Bconnec[19] << " " << Bconnec[9] << " "
                //          << Bconnec[7] << " " << Bconnec[17] << " " << Bconnec[25] << " "
                //          << Bconnec[15] << " " << Bconnec[3] << " " << Bconnec[5] << " "
                //          << Bconnec[23] << " " << Bconnec[21] << std::endl;

                output_v << Bconnec[18] << " " << Bconnec[20] << " " <<  Bconnec[26] << " " << Bconnec[24] << " " << Bconnec[0] << " " <<
                            Bconnec[2]  << " " << Bconnec[8]  << " " <<  Bconnec[6]  << " " << Bconnec[19] << " " << Bconnec[23] << " " <<
                            Bconnec[25] << " " << Bconnec[21] << " " <<  Bconnec[1]  << " " << Bconnec[5]  << " " << Bconnec[7] << " " << 
                            Bconnec[3] << " " << Bconnec[9]  << " " <<  Bconnec[11] << " " << Bconnec[15] << " " << Bconnec[17] << " " <<
                            Bconnec[12] << " " << Bconnec[14] << " " <<  Bconnec[10] << " " << Bconnec[16] << " " << Bconnec[22] << " " <<
                            Bconnec[4]  << " " << Bconnec[13] << std::endl;

            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            int aux = 0;
            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << aux + LNNC << std::endl;
                aux += LNNC;
            };
            output_v << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemCoarse; i++)
            {
                output_v << 72 << std::endl;
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (coarseModel.printVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Vel[i][0] << " "
                             << Vel[i][1] << " "
                             << Vel[i][2] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printRealVelocity)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << realVel[i][0] << " "
                             << realVel[i][1] << " "
                             << realVel[i][2] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };


            if (coarseModel.printDistFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << Distance[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printEnergyWeightFunction)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_v << EnergyW[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;

            if (coarseModel.printProcess)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    output_v << domDecompFine.first[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (coarseModel.printGlueZone)
            {
                output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemCoarse; i++)
                {
                    if (elementsGlueZoneCoarse_[cont] == i)
                    {
                        output_v << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_v << 0.0 << std::endl;
                    };
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;

            // PRINT FINE MODEL RESULTS - FEM mesh
            
            // FINE MESH
            s = "FINEoutput" + result + ".vtu";
            std::fstream output_vf(s.c_str(), std::ios_base::out);

            int LNN = 18*DIM-27;

            numBezierNodes = fineModel.NumBezierNodes;

            output_vf << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numBezierNodes
                     << "\"  NumberOfCells=\"" << numElemFine
                     << "\">" << std::endl;

            // WRITE NODAL COORDINATES
            output_vf << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            // Bezier Extraction
            double VelF[numBezierNodes][DIM], realVelF[numBezierNodes][DIM], CoordF[numBezierNodes][DIM];
            double DistanceF[numBezierNodes], EnergyWF[numBezierNodes];

            for (int iElem = 0; iElem < numElemFine; iElem++)
            {

                int *Bconnec = elementsFine_[iElem]->getBezierConnectivity();

                int *connec = elementsFine_[iElem]->getConnectivity();

                double coord_[LNN][DIM], vel_[LNN][DIM], realvel_[LNN][DIM];
                double distance_[LNN], energyW_[LNN];

                // Data in the NURBS control points
                for (int i = 0; i < LNN; i++)
                {

                    double *x = nodesFine_[connec[i]]->getCoordinates();
                    for (int j = 0; j < DIM; j++)
                    {
                        coord_[i][j] = x[j];
                        vel_[i][j] = nodesFine_[connec[i]]->getVelocity(j);
                        realvel_[i][j] = nodesFine_[connec[i]]->getVelocityArlequin(j);
                    };

                    distance_[i] = nodesFine_[connec[i]]->getDistFunction();
                    energyW_[i] = nodesFine_[connec[i]]->getWeightFunction();
                };

                // interpolated values (Bézier variables)
                double Bcoord_[LNN][DIM] = {};
                double Bvel_[LNN][DIM] = {};
                double Brealvel_[LNN][DIM] = {};
                double Bdistance_[LNN] = {};
                double BenergyW_[LNN] = {};

                for (int i = 0; i < LNN; i++)
                {

                    QuadShapeFunction<DIM> shapeQuad;
                    double phi_[LNN], wpc[LNN], xsi[DIM];
                    for (int k = 0; k < LNN; k++)
                        wpc[k] = nodesFine_[connec[k]]->getWeightPC();
                    int *inc_ = nodesFine_[connec[LNN - 1]]->getINC();
                    int patch = elementsFine_[iElem]->getPatch();
                    for (int j = 0; j < DIM; j++)
                        xsi[j] = shapeQuad.ParametricCoordBezier(i, j);
                    shapeQuad.evaluateIso(xsi, phi_, wpc, inc_, IsoParFine, patch);

                    for (int j = 0; j < LNN; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            Bcoord_[i][k] += phi_[j] * coord_[j][k];
                            Bvel_[i][k] += phi_[j] * vel_[j][k];
                            Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };
                };

                for (int i = 0; i < LNN; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        CoordF[Bconnec[i]][j] = Bcoord_[i][j];
                        VelF[Bconnec[i]][j] = Bvel_[i][j];
                        realVelF[Bconnec[i]][j] = Brealvel_[i][j];
                    };
                    DistanceF[Bconnec[i]] = Bdistance_[i];
                    EnergyWF[Bconnec[i]] = BenergyW_[i];
                };

            }; // iElem

            for (int i = 0; i < numBezierNodes; i++)
            {
                output_vf << CoordF[i][0] << " " << CoordF[i][1] << " " << CoordF[i][2]+0.1 << std::endl;
            };

            output_vf << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;

            // WRITE ELEMENT CONNECTIVITY
            output_vf << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

            for (int iElem = 0; iElem < numElemFine; ++iElem)
            {

                int *Bconnec_ = elementsFine_[iElem]->getBezierConnectivity();
                int Bconnec[LNN];
                for (int i = 0; i < LNN; i++)
                    Bconnec[i] = Bconnec_[i];

                output_vf << Bconnec[18] << " " << Bconnec[20] << " " <<  Bconnec[26] << " " << Bconnec[24] << " " << Bconnec[0] << " " <<
                            Bconnec[2]  << " " << Bconnec[8]  << " " <<  Bconnec[6]  << " " << Bconnec[19] << " " << Bconnec[23] << " " <<
                            Bconnec[25] << " " << Bconnec[21] << " " <<  Bconnec[1]  << " " << Bconnec[5]  << " " << Bconnec[7] << " " << 
                            Bconnec[3] << " " << Bconnec[9]  << " " <<  Bconnec[11] << " " << Bconnec[15] << " " << Bconnec[17] << " " <<
                            Bconnec[12] << " " << Bconnec[14] << " " <<  Bconnec[10] << " " << Bconnec[16] << " " << Bconnec[22] << " " <<
                            Bconnec[4]  << " " << Bconnec[13] << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE OFFSETS IN DATA ARRAY
            output_vf << "      <DataArray type=\"Int32\""
                     << " Name=\"offsets\" format=\"ascii\">" << std::endl;
            aux = 0;
            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << aux + LNN << std::endl;
                aux += LNN;
            };
            output_vf << "      </DataArray>" << std::endl;

            // WRITE ELEMENT TYPES
            output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;

            for (int i = 0; i < numElemFine; i++)
            {
                output_vf << 72 << std::endl;
            };
            output_vf << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            // WRITE NODAL RESULTS
            output_vf << "    <PointData>" << std::endl;

            if (fineModel.printVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_vf << VelF[i][0] << " "
                             << VelF[i][1] << " "
                             << VelF[i][2] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printRealVelocity)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_vf << realVelF[i][0] << " "
                             << realVelF[i][1] << " "
                             << realVelF[i][2] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };


            if (fineModel.printDistFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_vf << DistanceF[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printEnergyWeightFunction)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numBezierNodes; i++)
                {
                    output_vf << EnergyWF[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </PointData>" << std::endl;

            // WRITE ELEMENT RESULTS
            output_vf << "    <CellData>" << std::endl;

            if (fineModel.printProcess)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < numElemFine; i++)
                {
                    output_vf << domDecompFine.first[i] << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            if (fineModel.printGlueZone)
            {
                output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                int cont = 0;
                for (int i = 0; i < numElemFine; i++)
                {
                    if (elementsGlueZoneFine_[cont] == i)
                    {
                        output_vf << 1.0 << std::endl;
                        cont += 1;
                    }
                    else
                    {
                        output_vf << 0.0 << std::endl;
                    };
                };
                output_vf << "      </DataArray> " << std::endl;
            };

            output_vf << "    </CellData>" << std::endl;

            // FINALIZE OUTPUT FILE
            output_vf << "  </Piece>" << std::endl
                     << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;


        }; // rank == 0
    };     // printfreq
};

template <>
void Arlequin<2>::printResultsIP_FEM_ISO(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    const int DIM = 2;

    // Coarse Mesh (IGA MESH)
    if (rank == 0)
    {
        std::string result;
        std::ostringstream convert;

        convert << step + 100000;
        result = convert.str();
        std::string s = "saidaCoarseIP" + result + ".vtu";

        std::fstream output_v(s.c_str(), std::ios_base::out);

        int numberIntPoints = elementsFine_[0]->getNumberOfIntegrationPointsSpecial_FEM();

        output_v << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\">" << std::endl;

        output_v << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNNC = 18 * DIM - 27;

        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                QuadShapeFunction<DIM> shapeQuad;
                double qxsiC[DIM];
                int indCoarseElem;

                // correspondent integration point in the coarse mesh
                for (int j = 0; j < DIM; j++)
                    qxsiC[j] = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCoordinatesValue_FEM(ip, j);
                // coarse element index
                indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_FEM(ip);

                int *connec = elementsCoarse_[indCoarseElem]->getConnectivity();

                // Computes nurbs basis functions
                int patch = elementsCoarse_[indCoarseElem]->getPatch();
                int *inc = nodesCoarse_[connec[LNNC - 1]]->getINC();
                double wpc[LNNC], phi_[LNNC];
                for (int k = 0; k < LNNC; k++)
                    wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                shapeQuad.evaluateIso(qxsiC, phi_, wpc, inc, IsoParCoarse, patch);

                double coord[DIM] = {};
                for (int j = 0; j < LNNC; j++)
                {
                    double *x = nodesCoarse_[connec[j]]->getCoordinates();
                    for (int k = 0; k < DIM; k++)
                        coord[k] += x[k] * phi_[j];
                };

                output_v << coord[0] << " " << coord[1] << " " << 0.2 << std::endl;

            }; // loop integration points
        };     // loop glue fine mesh elements

        output_v << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_v << "    <Cells>" << std::endl
                 << "      <DataArray type=\"Int32\" "
                 << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_v << numN << std::endl;
        }

        output_v << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_v << "      <DataArray type=\"Int32\""
                 << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        int aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << aux + 1 << std::endl;
            aux += 1;
        };

        output_v << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << 1 << std::endl;
        };

        output_v << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_v << "  </Piece>" << std::endl;
        output_v << "  </UnstructuredGrid>" << std::endl
                 << "</VTKFile>" << std::endl;

        // PRINT FINE MODEL RESULTS (FEM mesh)
        std::string f = "saidaFineIP" + result + ".vtu";
        std::fstream output_vf(f.c_str(), std::ios_base::out);

        output_vf << "<?xml version=\"1.0\"?>" << std::endl
                  << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                  << "  <UnstructuredGrid>" << std::endl
                  << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\">" << std::endl;

        // WRITE NODAL COORDINATES
        output_vf << "    <Points>" << std::endl
                  << "      <DataArray type=\"Float64\" "
                  << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNN = 4 * DIM - 2;
        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            SpecialQuadrature squad;
            double x1[LNN], x2[LNN];
            int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

            for (int j = 0; j < LNN; j++)
            {
                double *x = nodesFine_[connec[j]]->getCoordinates();
                x1[j] = x[0];
                x2[j] = x[1];
            };

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                double x_[DIM];
                x_[0] = squad.interpolateQuadraticVariableFem(x1, ip);
                x_[1] = squad.interpolateQuadraticVariableFem(x2, ip);

                output_vf << x_[0] << " " << x_[1] << " " << 0.2 << std::endl;
            };
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_vf << "    <Cells>" << std::endl
                  << "      <DataArray type=\"Int32\" "
                  << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_vf << numN << std::endl;
        }

        output_vf << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_vf << "      <DataArray type=\"Int32\""
                  << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << aux + 1 << std::endl;
            aux += 1;
        };

        output_vf << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                  << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << 1 << std::endl;
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_vf << "  </Piece>" << std::endl;
        output_vf << "  </UnstructuredGrid>" << std::endl
                  << "</VTKFile>" << std::endl;
    };
};

template <>
void Arlequin<3>::printResultsIP_FEM_ISO(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    const int DIM = 3;

    // Coarse Mesh (IGA MESH)
    if (rank == 0)
    {
        std::string result;
        std::ostringstream convert;

        convert << step + 100000;
        result = convert.str();
        std::string s = "saidaCoarseIP" + result + ".vtu";

        std::fstream output_v(s.c_str(), std::ios_base::out);

        int numberIntPoints = elementsFine_[0]->getNumberOfIntegrationPointsSpecial_FEM();

        output_v << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\">" << std::endl;

        output_v << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNNC = 18*DIM-27;

        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                QuadShapeFunction<DIM> shapeQuad;
                double qxsiC[DIM];
                int indCoarseElem;

                // correspondent integration point in the coarse mesh
                for (int j = 0; j < DIM; j++)
                    qxsiC[j] = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCoordinatesValue_FEM(ip, j);

                // coarse element index
                indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_FEM(ip);

                int *connec = elementsCoarse_[indCoarseElem]->getConnectivity();

                // Computes nurbs basis functions
                int patch = elementsCoarse_[indCoarseElem]->getPatch();
                int *inc = nodesCoarse_[connec[LNNC - 1]]->getINC();
                double wpc[LNNC], phi_[LNNC];
                for (int k = 0; k < LNNC; k++)
                    wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                shapeQuad.evaluateIso(qxsiC, phi_, wpc, inc, IsoParCoarse, patch);

                double coord[DIM] = {};
                for (int j = 0; j < LNNC; j++)
                {
                    double *x = nodesCoarse_[connec[j]]->getCoordinates();
                    
                    for (int k = 0; k < DIM; k++)
                        coord[k] += x[k] * phi_[j];
                };
                
                output_v << coord[0] << " " << coord[1] << " " << coord[2]+0.1<< std::endl;


            }; // loop integration points
        };     // loop glue fine mesh elements

        output_v << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_v << "    <Cells>" << std::endl
                 << "      <DataArray type=\"Int32\" "
                 << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_v << numN << std::endl;
        }

        output_v << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_v << "      <DataArray type=\"Int32\""
                 << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        int aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << aux + 1 << std::endl;
            aux += 1;
        };

        output_v << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << 1 << std::endl;
        };

        output_v << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_v << "  </Piece>" << std::endl;
        output_v << "  </UnstructuredGrid>" << std::endl
                 << "</VTKFile>" << std::endl;

        // PRINT FINE MODEL RESULTS (FEM mesh)
        std::string f = "saidaFineIP" + result + ".vtu";
        std::fstream output_vf(f.c_str(), std::ios_base::out);

        output_vf << "<?xml version=\"1.0\"?>" << std::endl
                  << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                  << "  <UnstructuredGrid>" << std::endl
                  << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\">" << std::endl;

        // WRITE NODAL COORDINATES
        output_vf << "    <Points>" << std::endl
                  << "      <DataArray type=\"Float64\" "
                  << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNN = 4*DIM-2;
        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            SpecialQuadrature squad;
            double x1[LNN], x2[LNN], x3[LNN];
            int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

            for (int j = 0; j < LNN; j++)
            {
                double *x = nodesFine_[connec[j]]->getCoordinates();
                x1[j] = x[0];
                x2[j] = x[1];
                x3[j] = x[2];
            };

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                double x_[DIM];
                x_[0] = squad.interpolateQuadraticVariableFem(x1, ip);
                x_[1] = squad.interpolateQuadraticVariableFem(x2, ip);
                x_[2] = squad.interpolateQuadraticVariableFem(x3, ip);

                output_vf << x_[0] << " " << x_[1] << " " << x_[2]+0.1 << std::endl;
            };
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_vf << "    <Cells>" << std::endl
                  << "      <DataArray type=\"Int32\" "
                  << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_vf << numN << std::endl;
        }

        output_vf << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_vf << "      <DataArray type=\"Int32\""
                  << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << aux + 1 << std::endl;
            aux += 1;
        };

        output_vf << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                  << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << 1 << std::endl;
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_vf << "  </Piece>" << std::endl;
        output_vf << "  </UnstructuredGrid>" << std::endl
                  << "</VTKFile>" << std::endl;
    };
};

template <>
void Arlequin<2>::printResultsIP_FEM_FEM(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int const DIM = 2;

    // Coarse Mesh (FEM MESH)
    if (rank == 0)
    {
        std::string result;
        std::ostringstream convert;

        convert << step + 100000;
        result = convert.str();
        std::string s = "saidaCoarseIP" + result + ".vtu";

        std::fstream output_v(s.c_str(), std::ios_base::out);

        int numberIntPoints = elementsFine_[0]->getNumberOfIntegrationPointsSpecial_FEM();

        output_v << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\">" << std::endl;

        output_v << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNNC = 4*DIM-2;

        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                QuadShapeFunction<DIM> shapeQuad;
                double qxsiC[DIM];
                int indCoarseElem;

                // correspondent integration point in the coarse mesh
                for (int j = 0; j < DIM; j++)
                    qxsiC[j] = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCoordinatesValue_FEM(ip, j);

                // coarse element index
                indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_FEM(ip);

                int *connec = elementsCoarse_[indCoarseElem]->getConnectivity();

                // Computes nurbs basis functions
                double phi_[LNNC];
                shapeQuad.evaluateFem(qxsiC, phi_);

                double coord[DIM] = {};
                for (int j = 0; j < LNNC; j++)
                {
                    double *x = nodesCoarse_[connec[j]]->getCoordinates();
                    for (int k = 0; k < DIM; k++)
                        coord[k] += x[k] * phi_[j];
                };
                
                output_v << coord[0] << " " << coord[1] << " " << 0.0<< std::endl;


            }; // loop integration points
        };     // loop glue fine mesh elements

        output_v << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_v << "    <Cells>" << std::endl
                 << "      <DataArray type=\"Int32\" "
                 << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_v << numN << std::endl;
        }

        output_v << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_v << "      <DataArray type=\"Int32\""
                 << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        int aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << aux + 1 << std::endl;
            aux += 1;
        };

        output_v << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << 1 << std::endl;
        };

        output_v << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_v << "  </Piece>" << std::endl;
        output_v << "  </UnstructuredGrid>" << std::endl
                 << "</VTKFile>" << std::endl;

        // PRINT FINE MODEL RESULTS (FEM mesh)
        std::string f = "saidaFineIP" + result + ".vtu";
        std::fstream output_vf(f.c_str(), std::ios_base::out);

        output_vf << "<?xml version=\"1.0\"?>" << std::endl
                  << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                  << "  <UnstructuredGrid>" << std::endl
                  << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\">" << std::endl;

        // WRITE NODAL COORDINATES
        output_vf << "    <Points>" << std::endl
                  << "      <DataArray type=\"Float64\" "
                  << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNN = 4*DIM-2;
        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            SpecialQuadrature squad;
            double x1[LNN], x2[LNN];
            int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

            for (int j = 0; j < LNN; j++)
            {
                double *x = nodesFine_[connec[j]]->getCoordinates();
                x1[j] = x[0];
                x2[j] = x[1];
            };

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                double x_[DIM];
                x_[0] = squad.interpolateQuadraticVariableFem(x1, ip);
                x_[1] = squad.interpolateQuadraticVariableFem(x2, ip);

                output_vf << x_[0] << " " << x_[1] << " " << 0.0 << std::endl;
            };
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_vf << "    <Cells>" << std::endl
                  << "      <DataArray type=\"Int32\" "
                  << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_vf << numN << std::endl;
        }

        output_vf << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_vf << "      <DataArray type=\"Int32\""
                  << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << aux + 1 << std::endl;
            aux += 1;
        };

        output_vf << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                  << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << 1 << std::endl;
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_vf << "  </Piece>" << std::endl;
        output_vf << "  </UnstructuredGrid>" << std::endl
                  << "</VTKFile>" << std::endl;
    };
};


template <>
void Arlequin<3>::printResultsIP_FEM_FEM(int step)
{

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int const DIM = 3;

    // Coarse Mesh (FEM MESH)
    if (rank == 0)
    {
        std::string result;
        std::ostringstream convert;

        convert << step + 100000;
        result = convert.str();
        std::string s = "saidaCoarseIP" + result + ".vtu";

        std::fstream output_v(s.c_str(), std::ios_base::out);

        int numberIntPoints = elementsFine_[0]->getNumberOfIntegrationPointsSpecial_FEM();

        output_v << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\">" << std::endl;

        output_v << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNNC = 4*DIM-2;

        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                QuadShapeFunction<DIM> shapeQuad;
                double qxsiC[DIM];
                int indCoarseElem;

                // correspondent integration point in the coarse mesh
                for (int j = 0; j < DIM; j++)
                    qxsiC[j] = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCoordinatesValue_FEM(ip, j);

                // coarse element index
                indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_FEM(ip);

                int *connec = elementsCoarse_[indCoarseElem]->getConnectivity();

                // Computes nurbs basis functions
                double phi_[LNNC];
                shapeQuad.evaluateFem(qxsiC, phi_);

                double coord[DIM] = {};
                for (int j = 0; j < LNNC; j++)
                {
                    double *x = nodesCoarse_[connec[j]]->getCoordinates();
                    for (int k = 0; k < DIM; k++)
                        coord[k] += x[k] * phi_[j];
                };
                
                output_v << coord[0] << " " << coord[1] << " " << coord[2]+.1<< std::endl;


            }; // loop integration points
        };     // loop glue fine mesh elements

        output_v << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_v << "    <Cells>" << std::endl
                 << "      <DataArray type=\"Int32\" "
                 << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_v << numN << std::endl;
        }

        output_v << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_v << "      <DataArray type=\"Int32\""
                 << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        int aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << aux + 1 << std::endl;
            aux += 1;
        };

        output_v << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << 1 << std::endl;
        };

        output_v << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_v << "  </Piece>" << std::endl;
        output_v << "  </UnstructuredGrid>" << std::endl
                 << "</VTKFile>" << std::endl;

        // PRINT FINE MODEL RESULTS (FEM mesh)
        std::string f = "saidaFineIP" + result + ".vtu";
        std::fstream output_vf(f.c_str(), std::ios_base::out);

        output_vf << "<?xml version=\"1.0\"?>" << std::endl
                  << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                  << "  <UnstructuredGrid>" << std::endl
                  << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\">" << std::endl;

        // WRITE NODAL COORDINATES
        output_vf << "    <Points>" << std::endl
                  << "      <DataArray type=\"Float64\" "
                  << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNN = 4*DIM-2;
        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            SpecialQuadrature squad;
            double x1[LNN], x2[LNN], x3[LNN];
            int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

            for (int j = 0; j < LNN; j++)
            {
                double *x = nodesFine_[connec[j]]->getCoordinates();
                x1[j] = x[0];
                x2[j] = x[1];
                x3[j] = x[2];
            };

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                double x_[DIM];
                x_[0] = squad.interpolateQuadraticVariableFem(x1, ip);
                x_[1] = squad.interpolateQuadraticVariableFem(x2, ip);
                x_[2] = squad.interpolateQuadraticVariableFem(x3, ip);

                output_vf << x_[0] << " " << x_[1] << " " << x_[2]+.1<< std::endl;
            };
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_vf << "    <Cells>" << std::endl
                  << "      <DataArray type=\"Int32\" "
                  << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_vf << numN << std::endl;
        }

        output_vf << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_vf << "      <DataArray type=\"Int32\""
                  << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << aux + 1 << std::endl;
            aux += 1;
        };

        output_vf << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                  << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << 1 << std::endl;
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_vf << "  </Piece>" << std::endl;
        output_vf << "  </UnstructuredGrid>" << std::endl
                  << "</VTKFile>" << std::endl;
    };
};

template <>
void Arlequin<2>::printResultsIP_ISO_ISO(int step)
{
    int const DIM = 2;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Coarse Mesh (IGA MESH)
    if (rank == 0)
    {
        std::string result;
        std::ostringstream convert;

        convert << step + 100000;
        result = convert.str();
        std::string s = "saidaCoarseIP" + result + ".vtu";

        std::fstream output_v(s.c_str(), std::ios_base::out);

        int numberIntPoints = elementsFine_[0]->getNumberOfIntegrationPointsSpecial_ISO();

        output_v << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\">" << std::endl;

        output_v << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNNC = 18*DIM-27;

        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                QuadShapeFunction<DIM> shapeQuad;
                double qxsiC[DIM];
                int indCoarseElem;

                // correspondent integration point in the coarse mesh
                for (int j = 0; j < DIM; j++)
                    qxsiC[j] = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCoordinatesValue_ISO(ip, j);

                // coarse element index
                indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_ISO(ip);

                int *connec = elementsCoarse_[indCoarseElem]->getConnectivity();

                // Computes nurbs basis functions
                int patch = elementsCoarse_[indCoarseElem]->getPatch();
                int *inc = nodesCoarse_[connec[LNNC - 1]]->getINC();
                double wpc[LNNC], phi_[LNNC];
                for (int k = 0; k < LNNC; k++)
                    wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                shapeQuad.evaluateIso(qxsiC, phi_, wpc, inc, IsoParCoarse, patch);

                double coord[DIM] = {};
                for (int j = 0; j < LNNC; j++)
                {
                    double *x = nodesCoarse_[connec[j]]->getCoordinates();
                    
                    for (int k = 0; k < DIM; k++)
                        coord[k] += x[k] * phi_[j];
                };
                
                output_v << coord[0] << " " << coord[1] << " " << 0.0 << std::endl;


            }; // loop integration points
        };     // loop glue fine mesh elements

        output_v << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_v << "    <Cells>" << std::endl
                 << "      <DataArray type=\"Int32\" "
                 << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_v << numN << std::endl;
        }

        output_v << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_v << "      <DataArray type=\"Int32\""
                 << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        int aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << aux + 1 << std::endl;
            aux += 1;
        };

        output_v << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << 1 << std::endl;
        };

        output_v << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_v << "  </Piece>" << std::endl;
        output_v << "  </UnstructuredGrid>" << std::endl
                 << "</VTKFile>" << std::endl;

        // PRINT FINE MODEL RESULTS (FEM mesh)
        std::string f = "saidaFineIP" + result + ".vtu";
        std::fstream output_vf(f.c_str(), std::ios_base::out);

        output_vf << "<?xml version=\"1.0\"?>" << std::endl
                  << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                  << "  <UnstructuredGrid>" << std::endl
                  << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\">" << std::endl;

        // WRITE NODAL COORDINATES
        output_vf << "    <Points>" << std::endl
                  << "      <DataArray type=\"Float64\" "
                  << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNN = 18*DIM-27;
        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            Quadrature quad;
            double x1[LNN], x2[LNN];
            int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

            for (int j = 0; j < LNN; j++)
            {
                double *x = nodesFine_[connec[j]]->getCoordinates();
                x1[j] = x[0];
                x2[j] = x[1];
            };

            double wpc[LNN];
            for (int i = 0; i < LNN; i++) wpc[i] = nodesFine_[connec[i]] -> getWeightPC();
            int *inc = nodesFine_[connec[LNN-1]] -> getINC();
            int patch = elementsFine_[elementsGlueZoneFine_[i]]->getPatch();

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                double x_[DIM];
                x_[0] = quad.interpolateQuadraticVariableIso(x1, ip, wpc, inc, IsoParFine, patch);
                x_[1] = quad.interpolateQuadraticVariableIso(x2, ip, wpc, inc, IsoParFine, patch);

                output_vf << x_[0] << " " << x_[1] << " " << 0.0 << std::endl;
            };
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_vf << "    <Cells>" << std::endl
                  << "      <DataArray type=\"Int32\" "
                  << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_vf << numN << std::endl;
        }

        output_vf << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_vf << "      <DataArray type=\"Int32\""
                  << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << aux + 1 << std::endl;
            aux += 1;
        };

        output_vf << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                  << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << 1 << std::endl;
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_vf << "  </Piece>" << std::endl;
        output_vf << "  </UnstructuredGrid>" << std::endl
                  << "</VTKFile>" << std::endl;
    };
};


template <>
void Arlequin<3>::printResultsIP_ISO_ISO(int step)
{
    int const DIM = 3;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Coarse Mesh (IGA MESH)
    if (rank == 0)
    {
        std::string result;
        std::ostringstream convert;

        convert << step + 100000;
        result = convert.str();
        std::string s = "saidaCoarseIP" + result + ".vtu";

        std::fstream output_v(s.c_str(), std::ios_base::out);

        int numberIntPoints = elementsFine_[0]->getNumberOfIntegrationPointsSpecial_ISO();

        output_v << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                 << "\">" << std::endl;

        output_v << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNNC = 18*DIM-27;

        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                QuadShapeFunction<DIM> shapeQuad;
                double qxsiC[DIM];
                int indCoarseElem;

                // correspondent integration point in the coarse mesh
                for (int j = 0; j < DIM; j++)
                    qxsiC[j] = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCoordinatesValue_ISO(ip, j);

                // coarse element index
                indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]]->getIntegPointCorrespondenceElement_ISO(ip);

                int *connec = elementsCoarse_[indCoarseElem]->getConnectivity();

                // Computes nurbs basis functions
                int patch = elementsCoarse_[indCoarseElem]->getPatch();
                int *inc = nodesCoarse_[connec[LNNC - 1]]->getINC();
                double wpc[LNNC], phi_[LNNC];
                for (int k = 0; k < LNNC; k++)
                    wpc[k] = nodesCoarse_[connec[k]]->getWeightPC();
                shapeQuad.evaluateIso(qxsiC, phi_, wpc, inc, IsoParCoarse, patch);

                double coord[DIM] = {};
                for (int j = 0; j < LNNC; j++)
                {
                    double *x = nodesCoarse_[connec[j]]->getCoordinates();
                    
                    for (int k = 0; k < DIM; k++)
                        coord[k] += x[k] * phi_[j];
                };
                
                output_v << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;


            }; // loop integration points
        };     // loop glue fine mesh elements

        output_v << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_v << "    <Cells>" << std::endl
                 << "      <DataArray type=\"Int32\" "
                 << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_v << numN << std::endl;
        }

        output_v << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_v << "      <DataArray type=\"Int32\""
                 << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        int aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << aux + 1 << std::endl;
            aux += 1;
        };

        output_v << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_v << 1 << std::endl;
        };

        output_v << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_v << "  </Piece>" << std::endl;
        output_v << "  </UnstructuredGrid>" << std::endl
                 << "</VTKFile>" << std::endl;

        // PRINT FINE MODEL RESULTS (FEM mesh)
        std::string f = "saidaFineIP" + result + ".vtu";
        std::fstream output_vf(f.c_str(), std::ios_base::out);

        output_vf << "<?xml version=\"1.0\"?>" << std::endl
                  << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                  << "  <UnstructuredGrid>" << std::endl
                  << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\"  NumberOfCells=\"" << numElemGlueZoneFine * numberIntPoints
                  << "\">" << std::endl;

        // WRITE NODAL COORDINATES
        output_vf << "    <Points>" << std::endl
                  << "      <DataArray type=\"Float64\" "
                  << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        int LNN = 18*DIM-27;
        for (int i = 0; i < numElemGlueZoneFine; i++)
        {

            Quadrature quad;
            double x1[LNN], x2[LNN], x3[LNN];
            int *connec = elementsFine_[elementsGlueZoneFine_[i]]->getConnectivity();

            for (int j = 0; j < LNN; j++)
            {
                double *x = nodesFine_[connec[j]]->getCoordinates();
                x1[j] = x[0];
                x2[j] = x[1];
                x3[j] = x[2];
            };

            double wpc[LNN];
            for (int i = 0; i < LNN; i++) wpc[i] = nodesFine_[connec[i]] -> getWeightPC();
            int *inc = nodesFine_[connec[LNN-1]] -> getINC();
            int patch = elementsFine_[elementsGlueZoneFine_[i]]->getPatch();

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                double x_[DIM];
                x_[0] = quad.interpolateQuadraticVariableIso(x1, ip, wpc, inc, IsoParFine, patch);
                x_[1] = quad.interpolateQuadraticVariableIso(x2, ip, wpc, inc, IsoParFine, patch);
                x_[2] = quad.interpolateQuadraticVariableIso(x3, ip, wpc, inc, IsoParFine, patch);

                output_vf << x_[0] << " " << x_[1] << " " << x_[2] << std::endl;
            };
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Points>" << std::endl;

        // WRITE ELEMENT CONNECTIVITY
        output_vf << "    <Cells>" << std::endl
                  << "      <DataArray type=\"Int32\" "
                  << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int numN = 0; numN < numElemGlueZoneFine * numberIntPoints; ++numN)
        {
            output_vf << numN << std::endl;
        }

        output_vf << "      </DataArray>" << std::endl;

        // WRITE OFFSETS IN DATA ARRAY
        output_vf << "      <DataArray type=\"Int32\""
                  << " Name=\"offsets\" format=\"ascii\">" << std::endl;

        aux = 0;
        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << aux + 1 << std::endl;
            aux += 1;
        };

        output_vf << "      </DataArray>" << std::endl;

        // WRITE ELEMENT TYPES
        output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                  << "format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine * numberIntPoints; i++)
        {
            output_vf << 1 << std::endl;
        };

        output_vf << "      </DataArray>" << std::endl
                  << "    </Cells>" << std::endl;

        // FINALIZE OUTPUT FILE
        output_vf << "  </Piece>" << std::endl;
        output_vf << "  </UnstructuredGrid>" << std::endl
                  << "</VTKFile>" << std::endl;
    };
};

template <int DIM>
void Arlequin<DIM>::dragAndLiftCoefficients_FEM(std::ofstream &dragLift, int &iTimeStep)
{

    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    double pitchingMomentCoefficient = 0.;
    double pMom = 0.;
    double per = 0.;
    double &velocityInf = parametersFine->getVelocityInf(0);
    double &rhoInf = parametersFine->getDensity();
    double &dTime = parametersFine->getTimeStep();

    for (int jel = 0; jel < numBoundElemFine; jel++)
    {

        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        double aux_Mom = 0.;
        double aux_Per = 0.;

        for (int i = 0; i < fineModel.numberOfLines; i++)
        {

            if (boundaryFine_[jel]->getBoundaryGroup() == fineModel.dragAndLiftBoundary[i])
            {

                int iel = boundaryFine_[jel]->getElement();

                // Recognizing the element side
                std::pair<std::vector<int>, std::vector<int>> elemBound;
                elemBound = elementsFine_[iel]->getElemSideInBoundary();
                int numofBoundaries = elemBound.first.size();
                int side;
                for (int j = 0; j < numofBoundaries; j++)
                {
                    if (elemBound.second[j] == fineModel.dragAndLiftBoundary[i])
                    {
                        side = elemBound.first[j];
                    };
                };

                elementsFine_[iel]->computeDragAndLiftForces_FEM(side, pDForce, pLForce, fDForce, fLForce,
                                                                 dForce, lForce, aux_Mom, aux_Per);

                pMom += aux_Mom;
                per += aux_Per;
            };
        };

        pressureDragCoefficient += pDForce /
                                   (0.5 * rhoInf * velocityInf * velocityInf * 1);
        pressureLiftCoefficient += pLForce /
                                   (0.5 * rhoInf * velocityInf * velocityInf * 1);
        frictionDragCoefficient += fDForce /
                                   (0.5 * rhoInf * velocityInf * velocityInf * 1);
        frictionLiftCoefficient += fLForce /
                                   (0.5 * rhoInf * velocityInf * velocityInf * 1);
        dragCoefficient += dForce /
                           (0.5 * rhoInf * velocityInf * velocityInf * 1);
        liftCoefficient += lForce /
                           (0.5 * rhoInf * velocityInf * velocityInf * 1);
    };

    pitchingMomentCoefficient = pMom / (rhoInf * velocityInf * velocityInf * per);

    if (rank == 0)
    {
        const int timeWidth = 15;
        const int numWidth = 15;
        dragLift << std::setprecision(5) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) << iTimeStep * dTime;
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::setw(numWidth) << pitchingMomentCoefficient;
        dragLift << std::endl;
    };

    return;
};

// template<>
// void Arlequin<2>::dragAndLiftCoefficientsISO(std::ofstream& dragLift, int &iTimeStep){

//     double dragCoefficient = 0.;
//     double liftCoefficient = 0.;
//     double pressureDragCoefficient = 0.;
//     double pressureLiftCoefficient = 0.;
//     double frictionDragCoefficient = 0.;
//     double frictionLiftCoefficient = 0.;
//     double pitchingMomentCoefficient = 0.;
//     double pMom = 0.;
//     double per = 0.;
//     double& velocityInf = parametersFine -> getVelocityInf(0);
//     double& dTime = parametersFine -> getTimeStep();
//     double& rhoInf = parametersFine ->getDensity();

//     for (int jel = 0; jel < numBoundElemFine; jel++){

//         double dForce = 0.;
//         double lForce = 0.;
//         double pDForce = 0.;
//         double pLForce = 0.;
//         double fDForce = 0.;
//         double fLForce = 0.;
//         double aux_Mom = 0.;
//         double aux_Per = 0.;

//         for (int i=0; i< fineModel.numberOfLines; i++){

//             if (boundaryFine_[jel] -> getBoundaryGroup() == fineModel.dragAndLiftBoundary[i]){
//                 int iel = boundaryFine_[jel] -> getElement();

//                 //Recognizing the element side
//                 std::pair<std::vector<int>, std::vector<int> > elemBound;
//                 elemBound = elementsFine_[iel] -> getElemSideInBoundary();
//                 int numofBoundaries = elemBound.first.size();
//                 int side;
//                 for (int j = 0; j < numofBoundaries; j++) {
//                     if (elemBound.second[j] == fineModel.dragAndLiftBoundary[i]) {
//                         side = elemBound.first[j];
//                     };
//                 }

//                 elementsFine_[iel] -> computeDragAndLiftForces_ISO(side, pDForce, pLForce, fDForce, fLForce,
//                                                              dForce, lForce, aux_Mom, aux_Per);

//                 pMom += aux_Mom;
//                 per += aux_Per;

//             };
//         };

//         pressureDragCoefficient += pDForce /
//             (0.5 * rhoInf * velocityInf * velocityInf * 1);
//         pressureLiftCoefficient += pLForce /
//             (0.5 * rhoInf * velocityInf * velocityInf * 1);
//         frictionDragCoefficient += fDForce /
//             (0.5 * rhoInf * velocityInf * velocityInf * 1);
//         frictionLiftCoefficient += fLForce /
//             (0.5 * rhoInf * velocityInf * velocityInf * 1);
//         dragCoefficient += dForce /
//             (0.5 * rhoInf * velocityInf * velocityInf * 1);
//         liftCoefficient += lForce /
//             (0.5 * rhoInf * velocityInf * velocityInf * 1);
//     };

//     pitchingMomentCoefficient = pMom / (rhoInf * velocityInf * velocityInf * per);

//     if (rank == 0) {
//         const int timeWidth = 15;
//         const int numWidth = 15;
//         dragLift << std::setprecision(5) << std::scientific;
//         dragLift << std::left << std::setw(timeWidth) << iTimeStep * dTime;
//         dragLift << std::setw(numWidth) << pressureDragCoefficient;
//         dragLift << std::setw(numWidth) << pressureLiftCoefficient;
//         dragLift << std::setw(numWidth) << frictionDragCoefficient;
//         dragLift << std::setw(numWidth) << frictionLiftCoefficient;
//         dragLift << std::setw(numWidth) << dragCoefficient;
//         dragLift << std::setw(numWidth) << liftCoefficient;
//         dragLift << std::setw(numWidth) << pitchingMomentCoefficient;
//         dragLift << std::endl;
//     };

// return;

// };

//------------------------------------------------------------------------------
//---------------------------SETS FLUID MODELS----------------------------------
//------------------------------------------------------------------------------
template <int DIM>
void Arlequin<DIM>::setFluidModels_FEM_ISO(FluidMesh &coarse, FluidMesh &fine)
{

    coarseModel = coarse;
    fineModel = fine;

    // Gets Fine and Coarse models basic information from fluid data
    nodesCoarse_ = coarseModel.nodes_;
    nodesFine_ = fineModel.nodes_;

    elementsCoarse_ = coarseModel.elements_;
    elementsFine_ = fineModel.elements_;

    boundaryCoarse_ = coarseModel.boundary_;
    boundaryFine_ = fineModel.boundary_;

    numElemCoarse = coarseModel.elements_.size();
    numElemFine = fineModel.elements_.size();
    numNodesCoarse = coarseModel.nodes_.size();
    numNodesFine = fineModel.nodes_.size();

    numBoundElemFine = boundaryFine_.size();
    numBoundElemCoarse = boundaryCoarse_.size();

    domDecompCoarse = coarseModel.getDomainDecomposition();
    domDecompFine = fineModel.getDomainDecomposition();

    parametersCoarse = elementsCoarse_[0]->getFluidParameters();
    parametersFine = elementsFine_[0]->getFluidParameters();

    IsoParCoarse = coarseModel.IsoPar_;

    // starting the program with maximum dissipation
    double &IS = parametersFine->getSpectralRadius();
    integScheme = IS;
    double dd = 0.;
    parametersCoarse->setSpectralRadius(dd);
    parametersFine->setSpectralRadius(dd);

    // Non coincidente number controlPoints in fine and coarse mesh
    NCNumberNodesF = numNodesFine;
    NCNumberNodesC = coarseModel.NCNumberNodes;

    // Defines fine model as true
    for (int i = 0; i < numElemFine; i++)
    {
        elementsFine_[i]->setModel(true);
    };
    // Defines coarse model elements as false
    for (int i = 0; i < numElemCoarse; i++)
    {
        elementsCoarse_[i]->setModel(false);
    };

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Boxes for the integration points search in coarse mesh
    setElementBoxes_ISO();

    // Signaled distance from the a defined fine boundary mesh
    setSignaledDistance_FEM_ISO();

    // Sets the fine and coarse nodes/elements in the gluing zone
    // Defines the lagrange multipliers mesh
    setGluingZone_FEM_ISO();

    // Computes the Weight function for all nodes and integration points
    setWeightFunction_FEM_ISO();

    printResults_FEM_ISO(-10);
};


template <int DIM>
void Arlequin<DIM>::setFluidModels_ISO_ISO(FluidMesh &coarse, FluidMesh &fine)
{

    coarseModel = coarse;
    fineModel = fine;

    // Gets Fine and Coarse models basic information from fluid data
    nodesCoarse_ = coarseModel.nodes_;
    nodesFine_ = fineModel.nodes_;

    elementsCoarse_ = coarseModel.elements_;
    elementsFine_ = fineModel.elements_;

    boundaryCoarse_ = coarseModel.boundary_;
    boundaryFine_ = fineModel.boundary_;

    numElemCoarse = coarseModel.elements_.size();
    numElemFine = fineModel.elements_.size();
    numNodesCoarse = coarseModel.nodes_.size();
    numNodesFine = fineModel.nodes_.size();

    numBoundElemFine = boundaryFine_.size();
    numBoundElemCoarse = boundaryCoarse_.size();

    domDecompCoarse = coarseModel.getDomainDecomposition();
    domDecompFine = fineModel.getDomainDecomposition();

    parametersCoarse = elementsCoarse_[0]->getFluidParameters();
    parametersFine = elementsFine_[0]->getFluidParameters();

    IsoParCoarse = coarseModel.IsoPar_;
    IsoParFine = fineModel.IsoPar_;

    // starting the program with maximum dissipation
    double &IS = parametersFine->getSpectralRadius();
    integScheme = IS;
    double dd = 0.;
    parametersCoarse->setSpectralRadius(dd);
    parametersFine->setSpectralRadius(dd);

    // Non coincidente number controlPoints in fine and coarse mesh
    NCNumberNodesF = fineModel.NCNumberNodes;
    NCNumberNodesC = coarseModel.NCNumberNodes;

    // Defines fine model as true
    for (int i = 0; i < numElemFine; i++)
    {
        elementsFine_[i]->setModel(true);
    };
    // Defines coarse model elements as false
    for (int i = 0; i < numElemCoarse; i++)
    {
        elementsCoarse_[i]->setModel(false);
    };

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Boxes for the integration points search in coarse mesh
    setElementBoxes_ISO();

    // Signaled distance from the a defined fine boundary mesh
    setSignaledDistance_ISO_ISO();

    // Sets the fine and coarse nodes/elements in the gluing zone
    // Defines the lagrange multipliers mesh
    setGluingZone_ISO_ISO();

    // Computes the Weight function for all nodes and integration points
    setWeightFunction_ISO_ISO();

    // Computes the Nodal correspondence between fine nodes/integration points and coarse elements
    setCorrespondenceFine_ISO_ISO();

    printResultsLaplace_ISO_ISO(-10);
    printResultsIP_ISO_ISO(0);

};

template <int DIM>
void Arlequin<DIM>::setFluidModels_FEM_FEM(FluidMesh &coarse, FluidMesh &fine)
{

    coarseModel = coarse;
    fineModel = fine;

    // Gets Fine and Coarse models basic information from fluid data
    nodesCoarse_ = coarseModel.nodes_;
    nodesFine_ = fineModel.nodes_;

    elementsCoarse_ = coarseModel.elements_;
    elementsFine_ = fineModel.elements_;

    boundaryCoarse_ = coarseModel.boundary_;
    boundaryFine_ = fineModel.boundary_;

    numElemCoarse = coarseModel.elements_.size();
    numElemFine = fineModel.elements_.size();
    numNodesCoarse = coarseModel.nodes_.size();
    numNodesFine = fineModel.nodes_.size();

    numBoundElemFine = boundaryFine_.size();
    numBoundElemCoarse = boundaryCoarse_.size();

    domDecompCoarse = coarseModel.getDomainDecomposition();
    domDecompFine = fineModel.getDomainDecomposition();

    parametersCoarse = elementsCoarse_[0]->getFluidParameters();
    parametersFine = elementsFine_[0]->getFluidParameters();

    // starting the program with maximum dissipation
    double &IS = parametersFine->getSpectralRadius();
    integScheme = IS;
    double dd = 0.;
    parametersCoarse->setSpectralRadius(dd);
    parametersFine->setSpectralRadius(dd);

    // Non coincidente number controlPoints in fine and coarse mesh
    NCNumberNodesF = numNodesFine;
    NCNumberNodesC = numNodesCoarse;

    // Defines fine model as true
    for (int i = 0; i < numElemFine; i++)
    {
        elementsFine_[i]->setModel(true);
    };
    // Defines coarse model elements as false
    for (int i = 0; i < numElemCoarse; i++)
    {
        elementsCoarse_[i]->setModel(false);
    };

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

 
    // Boxes for the integration points search in coarse mesh
    setElementBoxes_FEM(); 

    // Signaled distance from the a defined fine boundary mesh
    setSignaledDistance_FEM_FEM(); 

    // Sets the fine and coarse nodes/elements in the gluing zone
    // Defines the lagrange multipliers mesh
    setGluingZone_FEM_FEM();

    // Computes the Weight function for all nodes and integration points
    setWeightFunction_FEM_FEM();

    printResults_FEM_FEM(-10);

    // Computes the Nodal correspondence between fine nodes/integration points and coarse elements
    setCorrespondenceFine_FEM_FEM();

    printResultsIP_FEM_FEM(0);
};


template <int DIM>
void Arlequin<DIM>::setDirichletConstrain_FEM_ISO(std::vector<int> &dofTemp,std::vector<double> &dofValue)
{

// Coarse mesh (IGA elements)
    int LNNC = 18*DIM-27;
    for (int ielem = 0; ielem < numElemCoarse; ielem++)
    {
        int *connec = elementsCoarse_[ielem]->getConnectivity();
        
        for (int i = 0; i < LNNC; i++)
        {
            // velocity constrain
            int constrain = nodesCoarse_[connec[i]]->getConstrains(0);
            if ((constrain == 1) || (constrain == 3))
            {
                int newconi = nodesCoarse_[connec[i]]->getnewcon();

                dofTemp.push_back(newconi);
                dofValue.push_back(nodesCoarse_[connec[i]]->getConstrainValue(0));
            };
   
        };
    };


    // Fine mesh (FEM elements)
    int LNN = 4*DIM-2;
    for (int ielem = 0; ielem < numElemFine; ielem++)
    {
        int *connec = elementsFine_[ielem]->getConnectivity();
        
        for (int i = 0; i < LNN; i++)
        {
            int constrain = nodesFine_[connec[i]]->getConstrains(0);
            if ((constrain == 1) || (constrain == 3))
            {   
                dofTemp.push_back(connec[i]+NCNumberNodesC);
                dofValue.push_back(nodesFine_[connec[i]]->getConstrainValue(0));
            };
        };
    };
};


template <int DIM>
void Arlequin<DIM>::setDirichletConstrain_FEM_FEM(std::vector<int> &dofTemp,std::vector<double> &dofValue){

    
    // Coarse mesh (FEM elements)
    int LNNC = 4*DIM-2;
    for (int ielem = 0; ielem < numElemCoarse; ielem++)
    {
        int *connec = elementsCoarse_[ielem]->getConnectivity();
        
        for (int i = 0; i < LNNC; i++)
        {
            // velocity constrain
            int constrain = nodesCoarse_[connec[i]]->getConstrains(0);
            if ((constrain == 1) || (constrain == 3))
            {
                dofTemp.push_back(connec[i]);
                dofValue.push_back(nodesCoarse_[connec[i]]->getConstrainValue(0));
            };
   
        };
    };


    // Fine mesh (FEM elements)
    int LNN = 4*DIM-2;
    for (int ielem = 0; ielem < numElemFine; ielem++)
    {
        int *connec = elementsFine_[ielem]->getConnectivity();
        
        for (int i = 0; i < LNN; i++)
        {
            int constrain = nodesFine_[connec[i]]->getConstrains(0);
            if ((constrain == 1) || (constrain == 3))
            {
                dofTemp.push_back(connec[i]+NCNumberNodesC);
                dofValue.push_back(nodesFine_[connec[i]]->getConstrainValue(0));
            };
        };
    };

}

template <int DIM>
void Arlequin<DIM>::setDirichletConstrain_ISO_ISO(std::vector<int> &dofTemp,std::vector<double> &dofValue){

    
    // Coarse mesh (IGA elements)
    int LNNC = 18 * DIM - 27;
    for (int ielem = 0; ielem < numElemCoarse; ielem++)
    {
        int *connec = elementsCoarse_[ielem]->getConnectivity();
        
        for (int i = 0; i < LNNC; i++)
        {
            // velocity constrain
            int constrain = nodesCoarse_[connec[i]]->getConstrains(0);
            if ((constrain == 1) || (constrain == 3))
            {
                int newconi = nodesCoarse_[connec[i]]->getnewcon();

                dofTemp.push_back(newconi);
                dofValue.push_back(nodesCoarse_[connec[i]]->getConstrainValue(0));
            };
   
        };
    };


    // Fine mesh (IGA elements)
    int LNN = 18 * DIM - 27;
    for (int ielem = 0; ielem < numElemFine; ielem++)
    {
        int *connec = elementsFine_[ielem]->getConnectivity();
        
        for (int i = 0; i < LNN; i++)
        {
            int constrain = nodesFine_[connec[i]]->getConstrains(0);
            if ((constrain == 1) || (constrain == 3))
            {   
                int newconi = nodesFine_[connec[i]]->getnewcon();
                dofTemp.push_back(newconi+NCNumberNodesC);
                dofValue.push_back(nodesFine_[connec[i]]->getConstrainValue(0));
            };
        };
    };

}

template <int DIM>
void Arlequin<DIM>::setDirichletConstrainLaplace_FEM_ISO(std::vector<int> &dofTemp)
{

    // Coarse mesh (IGA elements)
    int LNNC = 18 * DIM - 27;
    for (int ielem = 0; ielem < numElemCoarse; ielem++)
    {
        int *connec = elementsCoarse_[ielem]->getConnectivity();
        for (int i = 0; i < LNNC; i++)
        {
            int newconi = nodesCoarse_[connec[i]]->getnewcon();
            // velocity constrain
            for (int j = 0; j < DIM; j++)
            {
                int constrain = nodesCoarse_[connec[i]]->getConstrains(j);
                if (constrain == 1)
                {
                    dofTemp.push_back(newconi * DIM + j);
                };
            };
        };
    };

    // Fine mesh (FEM elements)
    int LNN = 4*DIM-2;
    for (int ielem = 0; ielem < numElemFine; ielem++)
    {
        int *connec = elementsFine_[ielem]->getConnectivity();
        for (int i = 0; i < LNN; i++)
        {
            // velocity constrain
            for (int j = 0; j < DIM; j++)
            {
                int constrain = nodesFine_[connec[i]]->getConstrains(j);
                if (constrain == 1)
                {
                    dofTemp.push_back(connec[i] * DIM + DIM * NCNumberNodesC + j);
                };
            };
        };
    };

};

template <int DIM>
void Arlequin<DIM>::setMatVecValuesCoarse_ISO()
{

    int LNNC = 18 * DIM - 27;

    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        if (domDecompCoarse.first[jel] == rank)
        {

            int *connec = elementsCoarse_[jel]->getConnectivity();

            double **elemMatrix;
            elemMatrix = new double *[LNNC * (DIM + 1)]();
            for (int i = 0; i < LNNC * (DIM + 1); ++i)
                elemMatrix[i] = new double[LNNC * (DIM + 1)]();
            double elemVector[LNNC * (DIM + 1)] = {};

            elementsCoarse_[jel]->getTransientNavierStokes_ISO(elemMatrix, elemVector);

            for (int i = 0; i < LNNC; i++)
            {

                int newconi = nodesCoarse_[connec[i]]->getnewcon();

                for (int j = 0; j < LNNC; j++)
                {

                    int newconj = nodesCoarse_[connec[j]]->getnewcon();

                    for (int k = 0; k < DIM; k++)
                    {

                        for (int l = 0; l < DIM; l++)
                        {

                            int dof_i = DIM * newconi + k;
                            int dof_j = DIM * newconj + l;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[DIM * i + k][DIM * j + l],
                                                ADD_VALUES);
                        };

                        int dof_i = DIM * newconi + k;
                        int dof_j = DIM * NCNumberNodesC + newconj;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[DIM * i + k][DIM * LNNC + j],
                                            ADD_VALUES);

                        dof_i = DIM * NCNumberNodesC + newconi;
                        dof_j = DIM * newconj + k;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[DIM * LNNC + i][DIM * j + k],
                                            ADD_VALUES);
                    };

                    int dof_i = DIM * NCNumberNodesC + newconi;
                    int dof_j = DIM * NCNumberNodesC + newconj;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[DIM * LNNC + i][DIM * LNNC + j],
                                        ADD_VALUES);
                }; // loop j

                // Rhs vector

                for (int k = 0; k < DIM; k++)
                {

                    int dof_i = DIM * newconi + k;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVector[DIM * i + k],
                                        ADD_VALUES);
                };

                int dof_i = DIM * NCNumberNodesC + newconi;
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[DIM * LNNC + i],
                                    ADD_VALUES);

            }; // loop i

            for (int i = 0; i < (DIM + 1) * LNNC; ++i)
                delete[] elemMatrix[i];
            delete[] elemMatrix;
        };
    };
};



template <int DIM>
void Arlequin<DIM>::setMatVecValuesCoarseLaplace_ISO()
{

    int LNNC = 18*DIM-27;

    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        if (domDecompCoarse.first[jel] == rank)
        {

            int *connec = elementsCoarse_[jel]->getConnectivity();

            double **elemMatrix;
            elemMatrix = new double *[LNNC]();
            for (int i = 0; i < LNNC; ++i)
                elemMatrix[i] = new double[LNNC]();
            double elemVector[LNNC] = {};

            elementsCoarse_[jel]->getLaplace_ISO(elemMatrix, elemVector);
         

            for (int i = 0; i < LNNC; i++)
            {

                int newconi = nodesCoarse_[connec[i]]->getnewcon();

                for (int j = 0; j < LNNC; j++)
                {

                    int newconj = nodesCoarse_[connec[j]]->getnewcon();

                    int dof_i = newconi;
                    int dof_j = newconj;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[i][j],
                                        ADD_VALUES);


                }; // loop j

                // Rhs vector

                int dof_i = newconi;
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[i],
                                    ADD_VALUES);

            }; // loop i

            for (int i = 0; i < LNNC; ++i)
                delete[] elemMatrix[i];
            delete[] elemMatrix;
        };
    };
};

template <int DIM>
void Arlequin<DIM>::setMatVecValuesCoarseLaplace_FEM()
{

    int LNNC = 4*DIM-2;

    for (int jel = 0; jel < numElemCoarse; jel++)
    {

        if (domDecompCoarse.first[jel] == rank)
        {

            int *connec = elementsCoarse_[jel]->getConnectivity();

            double **elemMatrix;
            elemMatrix = new double *[LNNC]();
            for (int i = 0; i < LNNC; ++i)
                elemMatrix[i] = new double[LNNC]();
            double elemVector[LNNC] = {};

            elementsCoarse_[jel]->getLaplace_FEM(elemMatrix, elemVector);
         

            for (int i = 0; i < LNNC; i++)
            {

                for (int j = 0; j < LNNC; j++)
                {

                    int dof_i = connec[i];
                    int dof_j = connec[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[i][j],
                                        ADD_VALUES);


                }; // loop j

                // Rhs vector

                int dof_i = connec[i];
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[i],
                                    ADD_VALUES);

            }; // loop i

            for (int i = 0; i < LNNC; ++i)
                delete[] elemMatrix[i];
            delete[] elemMatrix;
        };
    };
};

template <int DIM>
void Arlequin<DIM>::setMatVecValuesFine_FEM()
{

    int LNN = 4*DIM-2;
    for (int jel = 0; jel < numElemFine; jel++)
    {

        if (domDecompFine.first[jel] == rank)
        {

            int *connec = elementsFine_[jel]->getConnectivity();

            double **elemMatrix;
            elemMatrix = new double *[(DIM + 1) * LNN]();
            for (int i = 0; i < (DIM + 1) * LNN; ++i)
                elemMatrix[i] = new double[(DIM + 1) * LNN]();
            double elemVector[(DIM + 1) * LNN] = {};

            elementsFine_[jel]->getTransientNavierStokes_FEM(elemMatrix, elemVector);

            // Disperse local contributions into the global matrix
            // Matrix K and C
            for (int i = 0; i < LNN; i++)
            {
                for (int j = 0; j < LNN; j++)
                {

                    for (int k = 0; k < DIM; k++)
                    {
                        for (int l = 0; l < DIM; l++)
                        {

                            int dof_i = DIM * connec[i] + (DIM + 1) * NCNumberNodesC + k;
                            int dof_j = DIM * connec[j] + (DIM + 1) * NCNumberNodesC + l;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[DIM * i + k][DIM * j + l],
                                                ADD_VALUES);
                        };

                        int dof_i = DIM * connec[i] + (DIM + 1) * NCNumberNodesC + k;
                        int dof_j = DIM * NCNumberNodesF + connec[j] + (DIM + 1) * NCNumberNodesC;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[DIM * i + k][DIM * LNN + j],
                                            ADD_VALUES);

                        dof_i = DIM * NCNumberNodesF + connec[i] + (DIM + 1) * NCNumberNodesC;
                        dof_j = DIM * connec[j] + (DIM + 1) * NCNumberNodesC + k;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[DIM * LNN + i][DIM * j + k],
                                            ADD_VALUES);
                    };

                    int dof_i = DIM * NCNumberNodesF + connec[i] + (DIM + 1) * NCNumberNodesC;
                    int dof_j = DIM * NCNumberNodesF + connec[j] + (DIM + 1) * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[DIM * LNN + i][DIM * LNN + j],
                                        ADD_VALUES);
                }; // loop j

                // Rhs vector
                for (int k = 0; k < DIM; k++)
                {

                    int dof_i = DIM * connec[i] + (DIM + 1) * NCNumberNodesC + k;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVector[DIM * i + k],
                                        ADD_VALUES);
                };

                int dof_i = DIM * NCNumberNodesF + connec[i] + (DIM + 1) * NCNumberNodesC;
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[DIM * LNN + i],
                                    ADD_VALUES);

            }; // loop i
            for (int i = 0; i < (DIM + 1) * LNN; ++i)
                delete[] elemMatrix[i];
            delete[] elemMatrix;

        }; // domain decomposition
    };     // Elements Fine
};


template <int DIM>
void Arlequin<DIM>::setMatVecValuesFineLaplace_FEM()
{

    int LNN = 4*DIM-2;
    for (int jel = 0; jel < numElemFine; jel++)
    {

        if (domDecompFine.first[jel] == rank)
        {

            int *connec = elementsFine_[jel]->getConnectivity();

            double **elemMatrix;
            elemMatrix = new double *[LNN]();
            for (int i = 0; i < LNN; ++i)
                elemMatrix[i] = new double[LNN]();
            double elemVector[LNN] = {};

            elementsFine_[jel]->getLaplace_FEM(elemMatrix, elemVector);

            // Disperse local contributions into the global matrix
            // Matrix K and C
            for (int i = 0; i < LNN; i++)
            {
                for (int j = 0; j < LNN; j++)
                {
                    int dof_i = connec[i] + NCNumberNodesC;
                    int dof_j = connec[j] + NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[i][j],
                                        ADD_VALUES);
                }; // loop j

                // Rhs vector

                int dof_i = connec[i] + NCNumberNodesC;
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[i],
                                    ADD_VALUES);

            }; // loop i
            for (int i = 0; i < LNN; ++i)
                delete[] elemMatrix[i];
            delete[] elemMatrix;

        }; // domain decomposition
    };     // Elements Fine
};


template <int DIM>
void Arlequin<DIM>::setMatVecValuesFineLaplace_ISO()
{

    int LNN = 18*DIM-27;
    for (int jel = 0; jel < numElemFine; jel++)
    {

        if (domDecompFine.first[jel] == rank)
        {

            int *connec = elementsFine_[jel]->getConnectivity();

            double **elemMatrix;
            elemMatrix = new double *[LNN]();
            for (int i = 0; i < LNN; ++i)
                elemMatrix[i] = new double[LNN]();
            double elemVector[LNN] = {};

            elementsFine_[jel]->getLaplace_ISO(elemMatrix, elemVector);

            // Disperse local contributions into the global matrix
            // Matrix K and C
            for (int i = 0; i < LNN; i++)
            {
                int newconi = nodesFine_[connec[i]]->getnewcon();

                for (int j = 0; j < LNN; j++)
                {
                    int newconj = nodesFine_[connec[j]]->getnewcon();
                    
                    int dof_i = newconi + NCNumberNodesC;
                    int dof_j = newconj + NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[i][j],
                                        ADD_VALUES);
          
                }; // loop j

                // Rhs vector
                int dof_i = newconi + NCNumberNodesC;
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[i],
                                    ADD_VALUES);
      

            }; // loop i
            for (int i = 0; i < LNN; ++i)
                delete[] elemMatrix[i];
            delete[] elemMatrix;

        }; // domain decomposition
    };     // Elements Fine
};



template <int DIM>
void Arlequin<DIM>::setMatVecValuesLagrangeFine_FEM_ISO()
{

    double &alpha_f = parametersFine->getAlphaF();
    double &gamma = parametersFine->getGamma();
    double &dTime = parametersFine->getTimeStep();
    double integ = alpha_f * gamma * dTime;

    int LNN = 4*DIM-2;

    for (int l = 0; l < numElemGlueZoneFine; l++)
    {

        int jel = elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank)
        {

            int *connec = elementsFine_[jel]->getConnectivity();
            int *connecL = glueZoneFine_[l]->getConnectivity();

            int numberIntPoints = elementsFine_[jel]->getNumberOfIntegrationPointsSpecial_FEM();

            for (int ip = 0; ip < numberIntPoints; ip++)
            {

                int iElemCoarse = elementsFine_[jel]->getIntegPointCorrespondenceElement_FEM(ip);
                int *connecC = elementsCoarse_[iElemCoarse]->getConnectivity();
                int patch = elementsCoarse_[iElemCoarse]->getPatch();

                // LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
                double **elemMatrixLag1;
                elemMatrixLag1 = new double *[DIM*LNN]();
                for (int i = 0; i < DIM*LNN; ++i)
                    elemMatrixLag1[i] = new double[DIM*LNN]();

                double elemVectorLag1_1[DIM*LNN] = {};
                double elemVectorLag1_2[DIM*LNN] = {};

                // tSUPG and tPSPG STABILIZATION
                double **jacobianNRMatrix;
                jacobianNRMatrix = new double *[(DIM+1)*LNN]();
                for (int i = 0; i < (DIM+1)*LNN; ++i)
                    jacobianNRMatrix[i] = new double[DIM*LNN]();
                double rhsVector[(DIM+1)*LNN] = {};

                //ARLEQUIN STABILIZATION MATRIXES
                double **elemStabMatrixD;
                elemStabMatrixD = new double *[DIM*LNN]();
                for (int i = 0; i < DIM*LNN; ++i)
                    elemStabMatrixD[i] = new double[DIM*LNN]();

                double **elemStabMatrix1;
                elemStabMatrix1 = new double *[DIM*LNN]();
                for (int i = 0; i < DIM*LNN; ++i)
                    elemStabMatrix1[i] = new double[(DIM+1)*LNN]();

                double elemStabVectorD[DIM*LNN] = {};
                double elemStabVector1[DIM*LNN] = {};

                elementsFine_[jel]->getLagrangeMultipliersSameMesh_FEM(ip, elemMatrixLag1, elemVectorLag1_1, elemVectorLag1_2);
                // elementsFine_[jel] -> getLagrangeMultipliersSameMesh_tSUPG_tPSPG_FEM(ip,jacobianNRMatrix,rhsVector);
                elementsFine_[jel]->getLagrangeMultipliersSameMeshArlqStab_FEM_ISO(ip, patch, nodesCoarse_, connecC, IsoParCoarse,
                                                                                   elemStabMatrixD, elemStabVectorD,
                                                                                   elemStabMatrix1, elemStabVector1);

                //  for (int i = 0; i < DIM*LNN; i++)
                // {
                //     // std::cout << elemVectorLag1_2[i] << std::endl;
                //     for (int j = 0; j < DIM*LNN; j++)
                //     {
                //         std::cout << elemMatrixLag1[i][j] << " ";
                //     }

                //     std::cout << std::endl;
                // }

                // for (int i = 0; i < DIM*LNN; i++)
                // {
                //     std::cout << elemVectorLag1_1[i] << std::endl;
                // };


                for (int i = 0; i < LNN; i++)
                {
                    for (int j = 0; j < LNN; j++)
                    {

                        for (int k = 0; k < DIM; k++)
                        {
                            for (int l = 0; l < DIM; l++)
                            {

                                // Lagrange Multipliers
                                int dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                                int dof_j = (DIM + 1) * NCNumberNodesC + DIM * connec[j] + l;
                                double value = integ * elemMatrixLag1[DIM * i + k][DIM * j + l];
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &value, ADD_VALUES);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &elemMatrixLag1[DIM * i + k][DIM * j + l],
                                                    ADD_VALUES);

                                // tSUPG
                                dof_i = (DIM + 1) * NCNumberNodesC + DIM * connec[i] + k;
                                dof_j = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[j] + l;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &jacobianNRMatrix[DIM * i + k][DIM * j + l], ADD_VALUES);

                                //Arlequin Stabilization from diagonal matrix
                                dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                                dof_j = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[j] + l;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemStabMatrixD[DIM*i+k][DIM*j+l],
                                                    ADD_VALUES);

                                // Stabilization Arlequin terms from fine mesh matrix
                                dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                                dof_j = (DIM + 1) * NCNumberNodesC + DIM * connec[j] + l;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemStabMatrix1[DIM * i + k][DIM * j + l],
                                                    ADD_VALUES);
                            };

                            // tPSPG
                            int dof_i = (DIM + 1) * NCNumberNodesC + DIM * NCNumberNodesF + connec[i];
                            int dof_j = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[j] + k;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &jacobianNRMatrix[DIM * LNN + i][DIM * j + k], ADD_VALUES);

                            dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                            dof_j = (DIM + 1) * NCNumberNodesC +  DIM * NCNumberNodesF + connec[j];
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemStabMatrix1[DIM * i + k][DIM * LNN + j], ADD_VALUES);


                        };

                    }; // j

                    for (int k = 0; k < DIM; k++)
                    {

                        // Lagrange Multipliers
                        int dof_i = (DIM + 1) * NCNumberNodesC + DIM * connec[i] + k;
                        ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_1[DIM*i+k], ADD_VALUES);

                        dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                        ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_2[DIM*i+k], ADD_VALUES);

                        // tSUPG
                        dof_i = (DIM + 1) * NCNumberNodesC + DIM * connec[i] + k;
                        ierr = VecSetValues(b, 1, &dof_i, &rhsVector[DIM * i + k], ADD_VALUES);

                       // Arlequin stabilization term from diagonal
                        dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                        ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[DIM * i + k], ADD_VALUES);

                        // Arlequin stabilization term from fine mesh
                        dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                        ierr = VecSetValues(b, 1, &dof_i, &elemStabVector1[DIM * i + k], ADD_VALUES);
                    };

                    //tPSPG
                    int dof_i = (DIM + 1) * NCNumberNodesC + DIM * NCNumberNodesF + connec[i];
                    ierr = VecSetValues(b, 1, &dof_i, &rhsVector[DIM * LNN + i], ADD_VALUES);

                }; // i

                for (int i = 0; i < DIM*LNN; ++i)
                {
                    delete[] elemMatrixLag1[i];
                    delete[] elemStabMatrixD[i];
                    delete[] elemStabMatrix1[i];
                };
                delete[] elemMatrixLag1;
                delete[] elemStabMatrixD;
                delete[] elemStabMatrix1;

                for (int i = 0; i < (DIM+1)*LNN; ++i)
                {
                    delete[] jacobianNRMatrix[i];
                };
                delete[] jacobianNRMatrix;

            }; // integration points

        }; // decomposition

    }; // gluezonefine
};



template <int DIM>
void Arlequin<DIM>::setMatVecValuesLagrangeFineLaplace_FEM_ISO()
{

    int LNN = 4*DIM-2;

    for (int l = 0; l < numElemGlueZoneFine; l++)
    {

        int jel = elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank)
        {

            int *connec = elementsFine_[jel]->getConnectivity();
            int *connecL = glueZoneFine_[l]->getConnectivity();


            // LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
            double **elemMatrixLag1;
            elemMatrixLag1 = new double *[LNN]();
            for (int i = 0; i < LNN; ++i)
                elemMatrixLag1[i] = new double[LNN]();

            double elemVectorLag1_1[LNN] = {};
            double elemVectorLag1_2[LNN] = {};

            elementsFine_[jel]->getLagrangeMultipliersSameMeshLaplace_FEM(elemMatrixLag1, elemVectorLag1_1, elemVectorLag1_2);


            for (int i = 0; i < LNN; i++)
            {
                for (int j = 0; j < LNN; j++)
                {

                    // Lagrange Multipliers
                    int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                    int dof_j = NCNumberNodesC + connec[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrixLag1[i][j], ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag1[i][j],
                                        ADD_VALUES);

                }; // j

                // // Lagrange Multipliers
                // int dof_i = NCNumberNodesC + connec[i];
                // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_1[i], ADD_VALUES);

                // dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_2[i], ADD_VALUES);


            }; // i

            for (int i = 0; i < LNN; ++i)
            {
                delete[] elemMatrixLag1[i];
            };
            delete[] elemMatrixLag1;

        }; // decomposition

    }; // gluezonefine
};


template <int DIM>
void Arlequin<DIM>::setMatVecValuesLagrangeFineLaplace_FEM_FEM()
{

    int LNN = 4*DIM-2;

    for (int l = 0; l < numElemGlueZoneFine; l++)
    {

        int jel = elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank)
        {

            int *connec = elementsFine_[jel]->getConnectivity();
            int *connecL = glueZoneFine_[l]->getConnectivity();


            // LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
            double **elemMatrixLag1;
            elemMatrixLag1 = new double *[LNN]();
            for (int i = 0; i < LNN; ++i)
                elemMatrixLag1[i] = new double[LNN]();

            double elemVectorLag1_1[LNN] = {};
            double elemVectorLag1_2[LNN] = {};

            elementsFine_[jel]->getLagrangeMultipliersSameMeshLaplace_FEM(elemMatrixLag1, elemVectorLag1_1, elemVectorLag1_2);


            for (int i = 0; i < LNN; i++)
            {
                for (int j = 0; j < LNN; j++)
                {

                    // Lagrange Multipliers
                    int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                    int dof_j = NCNumberNodesC + connec[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrixLag1[i][j], ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag1[i][j],
                                        ADD_VALUES);

                }; // j

                // // Lagrange Multipliers
                // int dof_i = NCNumberNodesC + connec[i];
                // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_1[i], ADD_VALUES);

                // dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_2[i], ADD_VALUES);


            }; // i

            for (int i = 0; i < LNN; ++i)
            {
                delete[] elemMatrixLag1[i];
            };
            delete[] elemMatrixLag1;

        }; // decomposition

    }; // gluezonefine
};



template <int DIM>
void Arlequin<DIM>::setMatVecValuesLagrangeFineLaplace_ISO_ISO()
{

    int LNN = 18*DIM-27;

    for (int l = 0; l < numElemGlueZoneFine; l++)
    {

        int jel = elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank)
        {

            int *connec = elementsFine_[jel]->getConnectivity();
            int *connecL = glueZoneFine_[l]->getConnectivity();


            // LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
            double **elemMatrixLag1;
            elemMatrixLag1 = new double *[LNN]();
            for (int i = 0; i < LNN; ++i)
                elemMatrixLag1[i] = new double[LNN]();

            double elemVectorLag1_1[LNN] = {};
            double elemVectorLag1_2[LNN] = {};

            elementsFine_[jel]->getLagrangeMultipliersSameMeshLaplace_ISO(elemMatrixLag1, elemVectorLag1_1, elemVectorLag1_2);

            //ARRUMAR PARA NEWCONI
            for (int i = 0; i < LNN; i++)
            {
                for (int j = 0; j < LNN; j++)
                {

                    // Lagrange Multipliers
                    int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                    int dof_j = NCNumberNodesC + connec[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrixLag1[i][j], ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag1[i][j],
                                        ADD_VALUES);

    
                }; // j

                    // // Lagrange Multipliers
                    // int dof_i = NCNumberNodesC + connec[i];
                    // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_1[i], ADD_VALUES);

                    // dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                    // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_2[i], ADD_VALUES);

            }; // i

            for (int i = 0; i < LNN; ++i)
            {
                delete[] elemMatrixLag1[i];
            };
            delete[] elemMatrixLag1;

        }; // decomposition

    }; // gluezonefine
};


template <int DIM>
void Arlequin<DIM>::setMatVecValuesLagrangeCoarse_FEM_ISO()
{

    double &alpha_f = parametersFine->getAlphaF();
    double &gamma = parametersFine->getGamma();
    double &dTime = parametersFine->getTimeStep();

    double integ = alpha_f * gamma * dTime;

    int LNN = 4*DIM-2;
    int LNNC = 18*DIM-27;

    // Coarse mesh
    for (int l = 0; l < numElemGlueZoneFine; l++)
    {

        int jel = elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank)
        {

            int *connecL = glueZoneFine_[l]->getConnectivity();

            int numberIntPoints = elementsFine_[jel]->getNumberOfIntegrationPointsSpecial_FEM();

            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            // Finding coarse elements
            for (int i = 0; i < numberIntPoints; i++)
            {
                int aux = elementsFine_[jel]->getIntegPointCorrespondenceElement_FEM(i);
                ele.push_back(aux);
            };

            int numElemIntersect = 1;
            int flag = 0;
            diffElem.push_back(ele[0]);

            for (int i = 1; i < numberIntPoints; i++)
            {
                flag = 0;
                for (int j = 0; j < numElemIntersect; j++)
                {
                    if (ele[i] == diffElem[j])
                    {
                        break;
                    }
                    else
                    {
                        flag++;
                    };
                    if (flag == numElemIntersect)
                    {
                        numElemIntersect++;
                        diffElem.push_back(ele[i]);
                    };
                };
            };

            for (int ielem = 0; ielem < numElemIntersect; ielem++)
            {

                // LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
                double **elemMatrixLag0;
                elemMatrixLag0 = new double *[DIM * LNN]();
                for (int i = 0; i < DIM * LNN; ++i)
                    elemMatrixLag0[i] = new double[DIM * LNNC]();
                double elemVectorLag0_1[DIM * LNNC] = {};
                double elemVectorLag0_2[DIM * LNN] = {};

                // tSUPG and tPSPG STABILIZATION
                double **jacobianNRMatrix;
                jacobianNRMatrix = new double *[(DIM + 1) * LNNC]();
                for (int i = 0; i < (DIM + 1) * LNNC; ++i)
                    jacobianNRMatrix[i] = new double[DIM * LNN]();
                double rhsVector[(DIM + 1) * LNNC] = {};

                //ARLEQUIN STABILIZATION MATRIX AND VECTOR
                double **elemStabMatrixD;
                elemStabMatrixD = new double *[DIM * LNN]();
                for (int i = 0; i < DIM * LNN; ++i)
                    elemStabMatrixD[i] = new double[DIM * LNN]();
                double elemStabVectorD[DIM * LNN] = {};

                double **elemStabMatrix0;
                elemStabMatrix0 = new double *[DIM * LNN]();
                for (int i = 0; i < DIM * LNN; ++i)
                    elemStabMatrix0[i] = new double[(DIM + 1) * LNNC]();
                double elemStabVector0[DIM * LNN] = {};

                int iElemCoarse = diffElem[ielem];

                int *connecC = elementsCoarse_[iElemCoarse]->getConnectivity();
                int patch = elementsCoarse_[iElemCoarse]->getPatch();

                elementsFine_[jel]->getLagrangeMultipliersDifferentMesh_FEM_ISO(patch, nodesCoarse_, connecC, IsoParCoarse, iElemCoarse,
                                                                                elemMatrixLag0, elemVectorLag0_1, elemVectorLag0_2);
                // elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_FEM_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,
                //  																			  iElemCoarse, jacobianNRMatrix,rhsVector);
                elementsFine_[jel]->getLagrangeMultipliersDifferentMeshArlqStab_FEM_ISO(patch, nodesCoarse_, connecC, IsoParCoarse, iElemCoarse,
                                                                                        elemStabMatrixD, elemStabVectorD, elemStabMatrix0,
                                                                                        elemStabVector0);

            
                for (int i = 0; i < LNN; i++)
                {

                    for (int j = 0; j < LNNC; j++)
                    {

                        int newconj = nodesCoarse_[connecC[j]]->getnewcon();

                        for (int k = 0; k < DIM; k++)
                        {
                            for (int l = 0; l < DIM; l++)
                            {

                                // Lagrange multipliers matrixes
                                int dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                                int dof_j = DIM * newconj + l;
                                double value = integ * elemMatrixLag0[DIM * i + k][DIM * j + l];
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &value, ADD_VALUES);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &elemMatrixLag0[DIM * i + k][DIM * j + l],
                                                    ADD_VALUES);

                                // Arlequin Stabilization terms coarse mesh
                                dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                                dof_j = DIM * newconj + l;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemStabMatrix0[DIM * i + k][DIM * j + l],
                                                    ADD_VALUES);
                            };

                            int dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                            int dof_j = DIM * NCNumberNodesC + newconj;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemStabMatrix0[DIM * i + k][DIM * LNNC + j],
                                                ADD_VALUES);
                        };

                    }; // j

                    for (int k = 0; k < DIM; k++)
                    {

                        // Lagrange multipliers vector
                        int dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                        ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[DIM * i + k], ADD_VALUES);

                        // Arlequin Stabilization terms coarse mesh
                        dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                        ierr = VecSetValues(b, 1, &dof_i, &elemStabVector0[DIM * i + k], ADD_VALUES);
                    };

                }; // j

                for (int i = 0; i < LNNC; i++)
                {

                    int newconi = nodesCoarse_[connecC[i]]->getnewcon();

                    for (int j = 0; j < LNN; j++)
                    {

                        for (int k = 0; k < DIM; k++)
                        {

                            for (int l = 0; l < DIM; l++)
                            {

                                // tSUPG stabilization
                                int dof_i = DIM * newconi + k;
                                int dof_j = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[j] + l;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &jacobianNRMatrix[DIM * i + k][DIM * j + l],
                                                    ADD_VALUES);
                            };

                            // tPSPG stabilization
                            int dof_i = DIM * NCNumberNodesC + newconi;
                            int dof_j = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[j] + k;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &jacobianNRMatrix[DIM * LNNC + i][DIM * j + k],
                                                ADD_VALUES);
                        };

                    }; // j

                    for (int k = 0; k < DIM; k++)
                    {

                        // Lagrange multipliers vector
                        int dof_i = DIM * newconi + k;
                        ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[DIM * i + k], ADD_VALUES);

                        // tSUPG stabilization
                        dof_i = DIM * newconi + k;
                        ierr = VecSetValues(b, 1, &dof_i, &rhsVector[DIM * i + k], ADD_VALUES);
                    };

                    // tPSPG stabilization
                    int dof_i = DIM * NCNumberNodesC + newconi;
                    ierr = VecSetValues(b, 1, &dof_i, &rhsVector[DIM * LNNC + i], ADD_VALUES);

                }; // i

                //Lagrange stabilization terms (Lagrange multipliers defined in the fine mesh)
                for (int i = 0; i < LNN; i++)
                {
                    for (int j = 0; j < LNN; j++)
                    {

                        for (int k = 0; k < DIM; k++)
                        {
                            for (int l = 0; l < DIM; l++)
                            {

                                int dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                                int dof_j = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[j] + l;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemStabMatrixD[DIM * i + k][DIM * j + l],
                                                    ADD_VALUES);
                            };
                        };
                    };

                    for (int k = 0; k < DIM; k++)
                    {

                        int dof_i = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * connecL[i] + k;
                        ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[DIM * i + k], ADD_VALUES);
                    };
                };

                for (int i = 0; i < DIM * LNN; ++i)
                {
                    delete[] elemMatrixLag0[i];
                    delete[] elemStabMatrixD[i];
                    delete[] elemStabMatrix0[i];
                }
                delete[] elemMatrixLag0;
                delete[] elemStabMatrixD;
                delete[] elemStabMatrix0;

                for (int i = 0; i < (DIM + 1) * LNNC; ++i)
                    delete[] jacobianNRMatrix[i];
                delete[] jacobianNRMatrix;

            }; // intersect

        }; // decomposition

    }; // gluezonefine
};


template <int DIM>
void Arlequin<DIM>::setMatVecValuesLagrangeCoarseLaplace_FEM_ISO()
{
    int LNN = 4*DIM-2;
    int LNNC = 18*DIM-27;

    // Coarse mesh
    for (int l = 0; l < numElemGlueZoneFine; l++)
    {
        int jel = elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank)
        {

            int *connecL = glueZoneFine_[l]->getConnectivity();

            int numberIntPoints = elementsFine_[jel]->getNumberOfIntegrationPointsSpecial_FEM();

            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            // Finding coarse elements
            for (int i = 0; i < numberIntPoints; i++)
            {
                int aux = elementsFine_[jel]->getIntegPointCorrespondenceElement_FEM(i);
                ele.push_back(aux);
            };

            int numElemIntersect = 1;
            int flag = 0;
            diffElem.push_back(ele[0]);

            for (int i = 1; i < numberIntPoints; i++)
            {
                flag = 0;
                for (int j = 0; j < numElemIntersect; j++)
                {
                    if (ele[i] == diffElem[j])
                    {
                        break;
                    }
                    else
                    {
                        flag++;
                    };
                    if (flag == numElemIntersect)
                    {
                        numElemIntersect++;
                        diffElem.push_back(ele[i]);
                    };
                };
            };

            for (int ielem = 0; ielem < numElemIntersect; ielem++)
            {

                // LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
                double **elemMatrixLag0;
                elemMatrixLag0 = new double *[LNN]();
                for (int i = 0; i < LNN; ++i)
                    elemMatrixLag0[i] = new double[LNNC]();
                double elemVectorLag0_1[LNNC] = {};
                double elemVectorLag0_2[LNN] = {};

                int iElemCoarse = diffElem[ielem];

                int *connecC = elementsCoarse_[iElemCoarse]->getConnectivity();
                int patch = elementsCoarse_[iElemCoarse]->getPatch();

                elementsFine_[jel]->getLagrangeMultipliersDifferentMeshLaplace_FEM_ISO(patch, nodesCoarse_, connecC, IsoParCoarse, iElemCoarse,
                                                                                       elemMatrixLag0, elemVectorLag0_1, elemVectorLag0_2);

            
                for (int i = 0; i < LNN; i++)
                {
                    for (int j = 0; j < LNNC; j++)
                    {
                        int newconj = nodesCoarse_[connecC[j]]->getnewcon();

                        // Lagrange multipliers matrixes
                        int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                        int dof_j = newconj;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrixLag0[i][j], ADD_VALUES);
                        ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                            &elemMatrixLag0[i][j],
                                            ADD_VALUES);

                    }; // j


                    // Lagrange multipliers vector
                    int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[i], ADD_VALUES);

                }; // j

                for (int i = 0; i < LNNC; i++)
                {

                    int newconi = nodesCoarse_[connecC[i]]->getnewcon();

                    // Lagrange multipliers vector
                    int dof_i = newconi;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[i], ADD_VALUES);

                }; // i


                for (int i = 0; i < LNN; ++i)
                {
                    delete[] elemMatrixLag0[i];
                }
                delete[] elemMatrixLag0;

            }; // intersect

        }; // decomposition

    }; // gluezonefine
};


template <int DIM>
void Arlequin<DIM>::setMatVecValuesLagrangeCoarseLaplace_FEM_FEM()
{
    int LNN = 4*DIM-2;
    int LNNC = 4*DIM-2;

    // Coarse mesh
    for (int l = 0; l < numElemGlueZoneFine; l++)
    {
        int jel = elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank)
        {

            int *connecL = glueZoneFine_[l]->getConnectivity();

            int numberIntPoints = elementsFine_[jel]->getNumberOfIntegrationPointsSpecial_FEM();

            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            // Finding coarse elements
            for (int i = 0; i < numberIntPoints; i++)
            {
                int aux = elementsFine_[jel]->getIntegPointCorrespondenceElement_FEM(i);
                ele.push_back(aux);
            };

            int numElemIntersect = 1;
            int flag = 0;
            diffElem.push_back(ele[0]);

            for (int i = 1; i < numberIntPoints; i++)
            {
                flag = 0;
                for (int j = 0; j < numElemIntersect; j++)
                {
                    if (ele[i] == diffElem[j])
                    {
                        break;
                    }
                    else
                    {
                        flag++;
                    };
                    if (flag == numElemIntersect)
                    {
                        numElemIntersect++;
                        diffElem.push_back(ele[i]);
                    };
                };
            };

            for (int ielem = 0; ielem < numElemIntersect; ielem++)
            {

                // LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
                double **elemMatrixLag0;
                elemMatrixLag0 = new double *[LNN]();
                for (int i = 0; i < LNN; ++i)
                    elemMatrixLag0[i] = new double[LNNC]();
                double elemVectorLag0_1[LNNC] = {};
                double elemVectorLag0_2[LNN] = {};

                int iElemCoarse = diffElem[ielem];

                int *connecC = elementsCoarse_[iElemCoarse]->getConnectivity();


                elementsFine_[jel]->getLagrangeMultipliersDifferentMeshLaplace_FEM_FEM(nodesCoarse_, connecC, iElemCoarse,
                                                                                       elemMatrixLag0, elemVectorLag0_1, elemVectorLag0_2);              


                for (int i = 0; i < LNN; i++)
                {
                    for (int j = 0; j < LNNC; j++)
                    {

                        // Lagrange multipliers matrixes
                        int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                        int dof_j = connecC[j];
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrixLag0[i][j], ADD_VALUES);
                        ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                            &elemMatrixLag0[i][j],
                                            ADD_VALUES);
                    }; // j


                    // // Lagrange multipliers vector
                    // int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                    // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[i], ADD_VALUES);

                }; // j

                // for (int i = 0; i < LNNC; i++)
                // {
                //     // // Lagrange multipliers vector
                //     // int dof_i = connecC[i];
                //     // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[i], ADD_VALUES);
                // }; // i


                for (int i = 0; i < LNN; ++i)
                {
                    delete[] elemMatrixLag0[i];
                }
                delete[] elemMatrixLag0;

            }; // intersect

        }; // decomposition

    }; // gluezonefine
};

template <int DIM>
void Arlequin<DIM>::setMatVecValuesLagrangeCoarseLaplace_ISO_ISO()
{
    int LNN = 18*DIM-27;
    int LNNC = 18*DIM-27;

    // Coarse mesh
    for (int l = 0; l < numElemGlueZoneFine; l++)
    {
        int jel = elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank)
        {

            int *connecL = glueZoneFine_[l]->getConnectivity();

            int numberIntPoints = elementsFine_[jel]->getNumberOfIntegrationPointsSpecial_ISO();

            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            // Finding coarse elements
            for (int i = 0; i < numberIntPoints; i++)
            {
                int aux = elementsFine_[jel]->getIntegPointCorrespondenceElement_ISO(i);
                ele.push_back(aux);
            };

            int numElemIntersect = 1;
            int flag = 0;
            diffElem.push_back(ele[0]);

            for (int i = 1; i < numberIntPoints; i++)
            {
                flag = 0;
                for (int j = 0; j < numElemIntersect; j++)
                {
                    if (ele[i] == diffElem[j])
                    {
                        break;
                    }
                    else
                    {
                        flag++;
                    };
                    if (flag == numElemIntersect)
                    {
                        numElemIntersect++;
                        diffElem.push_back(ele[i]);
                    };
                };
            };

            for (int ielem = 0; ielem < numElemIntersect; ielem++)
            {

                // LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
                double **elemMatrixLag0;
                elemMatrixLag0 = new double *[LNN]();
                for (int i = 0; i < LNN; ++i)
                    elemMatrixLag0[i] = new double[LNNC]();
                double elemVectorLag0_1[LNNC] = {};
                double elemVectorLag0_2[LNN] = {};

                int iElemCoarse = diffElem[ielem];

                int *connecC = elementsCoarse_[iElemCoarse]->getConnectivity();
                int patch = elementsCoarse_[iElemCoarse]->getPatch();

                elementsFine_[jel]->getLagrangeMultipliersDifferentMeshLaplace_ISO_ISO(patch, nodesCoarse_, connecC, IsoParCoarse, iElemCoarse,
                                                                                       elemMatrixLag0, elemVectorLag0_1, elemVectorLag0_2);

            
                for (int i = 0; i < LNN; i++)
                {

                    for (int j = 0; j < LNNC; j++)
                    {

                        int newconj = nodesCoarse_[connecC[j]]->getnewcon();

                        // Lagrange multipliers matrixes
                        int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                        int dof_j = newconj;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrixLag0[i][j], ADD_VALUES);
                        ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                            &elemMatrixLag0[i][j],
                                            ADD_VALUES);

                    }; // j

                    // // Lagrange multipliers vector
                    // int dof_i = NCNumberNodesC + NCNumberNodesF + connecL[i];
                    // ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[i], ADD_VALUES);

                }; // j

                // for (int i = 0; i < LNNC; i++)
                // {

                //     int newconi = nodesCoarse_[connecC[i]]->getnewcon();

                //     // Lagrange multipliers vector
                //     int dof_i = newconi;
                //     ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[i], ADD_VALUES);

                // }; // i


                for (int i = 0; i < LNN; ++i)
                {
                    delete[] elemMatrixLag0[i];
                }
                delete[] elemMatrixLag0;

            }; // intersect

        }; // decomposition

    }; // gluezonefine
};

template <int DIM>
int Arlequin<DIM>::solveArlequinProblem_FEM_ISO(int iterNumber, double tolerance)
{
    
    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    if (rank == 0)
    {
        dragLift << "Time   Pressure Drag   Pressure Lift "
                 << "Friction Drag  Friction Lift Drag    Lift "
                 << std::endl;
    };

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Computes the Nodal correspondence between fine nodes/integration points and coarse elements
    setCorrespondenceFine_FEM_ISO();

    printResultsIP_FEM_ISO(-10);

    // Set the degre of freedom index with boundary condition
    std::vector<int> dofTemp;
    // setDirichletConstrain_FEM_ISO(dofTemp);
    PetscMalloc1(dofTemp.size(), &dof);
    for (size_t i = 0; i < dofTemp.size(); i++)
    {
        dof[i] = dofTemp[i];
    };

    int sysSize = (DIM + 1) * NCNumberNodesC + (DIM + 1) * NCNumberNodesF + DIM * NCnumNodesGlueZoneFine;

    double sumtime = 0;

    int numTimeSteps = fineModel.numTimeSteps;
    double &dTime = parametersFine->getTimeStep();
    int iTimeStep;

    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++)
    {

        if (rank == 0)
        {
            std::cout << "------------------------- TIME STEP = "
                      << iTimeStep << " -------------------------"
                      << std::endl;
        }

        if (iTimeStep == 10)
        {
            parametersCoarse->setSpectralRadius(integScheme);
            parametersFine->setSpectralRadius(integScheme);
        };

        double &gamma = parametersCoarse->getGamma();

        for (int i = 0; i < numNodesCoarse; i++)
        {

            double accel[DIM], u[DIM];

            for (int j = 0; j < DIM; j++)
                u[j] = nodesCoarse_[i]->getVelocity(j);
            nodesCoarse_[i]->setPreviousVelocity(u);

            for (int j = 0; j < DIM; j++)
                accel[j] = nodesCoarse_[i]->getAcceleration(j);
            nodesCoarse_[i]->setPreviousAcceleration(accel);
            for (int j = 0; j < DIM; j++)
                accel[j] *= (gamma - 1.) / gamma;
            nodesCoarse_[i]->setAcceleration(accel);
        };

        for (int i = 0; i < numNodesFine; i++)
        {

            double accel[DIM], u[DIM];

            for (int j = 0; j < DIM; j++)
                u[j] = nodesFine_[i]->getVelocity(j);
            nodesFine_[i]->setPreviousVelocity(u);
            for (int j = 0; j < DIM; j++)
                accel[j] = nodesFine_[i]->getAcceleration(j);

            nodesFine_[i]->setPreviousAcceleration(accel);
            for (int j = 0; j < DIM; j++)
                accel[j] *= (gamma - 1.) / gamma;
            nodesFine_[i]->setAcceleration(accel);
        };

        double duNorm = 100.;

        for (int inewton = 0; inewton < iterNumber; inewton++)
        {

            std::clock_t t1 = std::clock();

            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize,sysSize,100000,NULL,100000,NULL,&A);CHKERRQ(ierr);

            for (int i=0; i<sysSize; i++){
                double valu = 1.e-20;
                ierr = MatSetValues(A,1,&i,1,&i,&valu,ADD_VALUES);
            }

            // Divides the matrix between the processes
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);

            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,sysSize);CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All);CHKERRQ(ierr);

            //Matrix and vectors - COARSE MESH - IGA mesh
            setMatVecValuesCoarse_ISO();
            //Matrix and vectors - Lagrange multiplieres - COARSE MESH
            setMatVecValuesLagrangeCoarse_FEM_ISO();

            // //Matrix and vectors -FINE MESH  - FEM mesh
            setMatVecValuesFine_FEM();
            // // //Matrix and vectors - Lagrange multiplieres - FINE MESH
            setMatVecValuesLagrangeFine_FEM_ISO();


            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

            //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

            MatZeroRowsColumns(A,dofTemp.size(),dof,1.0,u,b);

            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

            //Solve using GMRES
            // ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,500); CHKERRQ(ierr);
            // ierr = KSPGetPC(ksp,&pc);
            // ierr = PCSetType(pc,PCNONE);
            // ierr = KSPSetType(ksp,KSPDGMRES); CHKERRQ(ierr);
            // ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);

            // //Solve using Mumps
            #if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp,KSPPREONLY);
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCLU);
            #endif
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);

            //Solve Linear System
            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            ierr = KSPGetTotalIterations(ksp, &iterations);

                //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);

            //Updates nodal values

            double u_;
            double p_;
            double lag_;
            double normU = 0.;
            double normP = 0.;
            double normL = 0.;
            Ione = 1;

            for (int i = 0; i < numNodesCoarse; ++i){
                int newconi = nodesCoarse_[i] -> getnewcon();
                for (int j = 0; j < DIM; j++){
                    Ii = DIM*newconi+j;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    u_ = val;
                    nodesCoarse_[i] -> incrementAcceleration(j,u_);
                    nodesCoarse_[i] -> incrementVelocity(j,u_*gamma*dTime);
                    normU += val*val;
                };
                Ii = DIM*NCNumberNodesC + newconi;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                normP += val*val;
                nodesCoarse_[i] -> incrementPressure(p_);
            };

            for (int i = 0; i < numNodesFine; ++i){
                for (int j = 0; j < DIM; j++){
                    Ii = DIM*i + (DIM+1)*NCNumberNodesC + j;
                    ierr = VecGetValues(All, Ione, &Ii, &val);
                    u_ = val;
                    nodesFine_[i] -> incrementAcceleration(j,u_);
                    nodesFine_[i] -> incrementVelocity(j,u_*gamma*dTime);
                    normU += val*val;
                };
                Ii = DIM*NCNumberNodesF + i + (DIM+1)*NCNumberNodesC;
                ierr = VecGetValues(All,Ione,&Ii,&val);
                p_ = val;
                nodesFine_[i] -> incrementPressure(p_);
                normP += val*val;
            };

            for (int i = 0; i < NCnumNodesGlueZoneFine; ++i){
                for (int j = 0; j < DIM; j++){
                    Ii = (DIM+1)*NCNumberNodesF + (DIM+1)*NCNumberNodesC + DIM*i + j;
                    ierr = VecGetValues(All, Ione, &Ii, &val);
                    lag_ = val;
                    normL += val * val;
                    nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(j,lag_);
                };
            };

            std::clock_t t2 = std::clock();
            sumtime += 1000.*(t2-t1)/CLOCKS_PER_SEC/1000.;

            if(rank == 0){

                std::cout << "Iteration = " << inewton
                            << " (" << iterations << ")"
                            << "   Du Norm = " << std::scientific << sqrt(normU)
                            << " " << sqrt(normP)  << " " << sqrt(normL)
                            << "  Time (s) = " << std::fixed
                            << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
            };

            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);
            // if(iTimeStep<5){normU=0.;}
            // if(iTimeStep<100){if(inewton>1)normU=0.;}
            if (sqrt(normU) <= tolerance) {
                break;
            };

        }; // Newton-Raphson

        if(rank == 0){
            std::cout << "Accumulated time = " << sumtime << std::endl;
        };

        //Computes the real velocity
        QuadShapeFunction<DIM> shapeQuad;

        for (int i = 0; i < numNodesFine; i++){
            for (int k = 0; k < DIM; k++) nodesFine_[i] -> setVelocityArlequin(k,nodesFine_[i] -> getVelocity(k));
            nodesFine_[i] -> setPressureArlequin(nodesFine_[i] ->getPressure());
        };

        int LNNC = 18*DIM - 27;

        for (int i = 0; i<numNodesGlueZoneFine; i++){

            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
            double* xsiC = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();

            int *connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();

            int patchC = elementsCoarse_[elCoarse] -> getPatch();
            int *incC = nodesCoarse_[connecCoarse[LNNC-1]] -> getINC();
            double wpcC[LNNC],phiC_[LNNC];
            for (int k = 0; k < LNNC; k++) wpcC[k] = nodesCoarse_[connecCoarse[k]] -> getWeightPC();
            shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC,IsoParCoarse,patchC);

            double u[DIM] = {};
            double p = 0;
            for (int j=0; j<LNNC; j++){
                for (int k = 0 ; k < DIM; k++){
                    u[k] += nodesCoarse_[connecCoarse[j]] -> getVelocity(k) * phiC_[j];
                };
                p += nodesCoarse_[connecCoarse[j]] -> getPressure() * phiC_[j];
            };

            double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> getWeightFunction();

            for (int k = 0; k < DIM; k++){
                double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> getVelocity(k) * wFunc + u[k] * (1. - wFunc);
                nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(k,u_int);
            };

            double p_int = nodesFine_[nodesGlueZoneFine_[i]] -> getPressure() * wFunc + p * (1. - wFunc);
            nodesFine_[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);

            for (int i=0; i<numNodesCoarse; i++){
                for (int k = 0; k < DIM; k++) nodesCoarse_[i] -> setVelocityArlequin(k,nodesCoarse_[i] -> getVelocity(k));
                nodesCoarse_[i] -> setPressureArlequin(nodesCoarse_[i] -> getPressure());
            };
        };

        //Compute and print drag and lift coefficients
        if (fineModel.computeDragAndLift){
            dragAndLiftCoefficients_FEM(dragLift,iTimeStep);
        };

        // Printing results
        printResults_FEM_ISO(iTimeStep);
    };

    PetscFree(dof);
    return 0;
};




template <int DIM>
int Arlequin<DIM>::solveArlequinProblemLaplace_FEM_ISO(int iterNumber, double tolerance)
{
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Computes the Nodal correspondence between fine nodes/integration points and coarse elements
    setCorrespondenceFine_FEM_ISO();

    printResultsIP_FEM_ISO(-10);

    std::vector<int> dofTemp;
    std::vector<double>dofValue;
    setDirichletConstrain_FEM_ISO(dofTemp,dofValue);
    PetscMalloc1(dofTemp.size(), &dof);
    for (size_t i = 0; i < dofTemp.size(); i++)
    {
        dof[i] = dofTemp[i];
    };

    std::clock_t t1 = std::clock();

    int sysSize = NCNumberNodesC + NCNumberNodesF + NCnumNodesGlueZoneFine;

    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        sysSize,sysSize,10000,NULL,10000,NULL,&A);CHKERRQ(ierr);

    for (int i=0; i<sysSize; i++){
        double valu = 1.e-20;
        ierr = MatSetValues(A,1,&i,1,&i,&valu,ADD_VALUES);
    }

    // Divides the matrix between the processes
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
    ierr = VecSetSizes(b,PETSC_DECIDE,sysSize);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&b2);CHKERRQ(ierr);
    ierr = VecSetSizes(b2,PETSC_DECIDE,sysSize);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b2);CHKERRQ(ierr);
       
    ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&All);CHKERRQ(ierr);

    for (int i = 0; i < dofTemp.size(); i++){

        int ind = dofTemp[i];
        ierr = VecSetValues(b2, 1, &ind, &dofValue[i],
                            INSERT_VALUES);
    };
   
    ierr = VecAssemblyBegin(b2);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b2);CHKERRQ(ierr);

    //Matrix and vectors - COARSE MESH - IGA mesh
    setMatVecValuesCoarseLaplace_ISO();
    // //Matrix and vectors - Lagrange multiplieres - COARSE MESH
    setMatVecValuesLagrangeCoarseLaplace_FEM_ISO();

    // //Matrix and vectors -FINE MESH  - FEM mesh
    setMatVecValuesFineLaplace_FEM();
    // //Matrix and vectors - Lagrange multiplieres - FINE MESH
    setMatVecValuesLagrangeFineLaplace_FEM_ISO();
   
    

    //Assemble matrices and vectors
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

    MatZeroRows(A,dofTemp.size(),dof,1.0,b2,b);

    //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    // MatZeroRowsColumns(A,dofTemp.size(),dof,1.0,u,b);

    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

    //Solve using GMRES
    // ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,500); CHKERRQ(ierr);
    // ierr = KSPGetPC(ksp,&pc);
    // ierr = PCSetType(pc,PCNONE);
    // ierr = KSPSetType(ksp,KSPDGMRES); CHKERRQ(ierr);
    // ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);

    // //Solve using Mumps
    #if defined(PETSC_HAVE_MUMPS)
    ierr = KSPSetType(ksp,KSPPREONLY);
    ierr = KSPGetPC(ksp,&pc);
    ierr = PCSetType(pc, PCLU);
    #endif
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);

    //Solve Linear System
    ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
    ierr = KSPGetTotalIterations(ksp, &iterations);

    // //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // //VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);

    // //Updates nodal values

    double u_;
    double lag_;
    double normU = 0.;
    double normL = 0.;
    Ione = 1;

    for (int i = 0; i < numNodesCoarse; ++i){
        int newconi = nodesCoarse_[i] -> getnewcon();
        Ii = newconi;
        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
        u_ = val;
        nodesCoarse_[i] -> setVelocityComponent(0,u_);

        normU += val*val;
    };

    for (int i = 0; i < numNodesFine; ++i){
        Ii = i + NCNumberNodesC;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        u_ = val;
        nodesFine_[i] -> setVelocityComponent(0,u_);
        normU += val*val;
    };

    for (int i = 0; i < NCnumNodesGlueZoneFine; ++i){
        Ii = NCNumberNodesF + NCNumberNodesC + i;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        lag_ = val;
        normL += val * val;
        nodesFine_[nodesGlueZoneFine_[i]] -> setLagrangeMultiplier(0,lag_);
    };

    std::clock_t t2 = std::clock();

    if(rank == 0){

        std::cout   << "   Du Norm = " << std::scientific << sqrt(normU)
                    << " " << sqrt(normL)
                    << "  Time (s) = " << std::fixed
                    << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
    };

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&u); CHKERRQ(ierr);
    ierr = VecDestroy(&All); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);


    //Computes the real velocity
    QuadShapeFunction<DIM> shapeQuad;

    for (int i = 0; i < numNodesFine; i++){
        for (int k = 0; k < DIM; k++) nodesFine_[i] -> setVelocityArlequin(k,nodesFine_[i] -> getVelocity(k));
    };

    int LNNC = 18*DIM-27;

    for (int i = 0; i<numNodesGlueZoneFine; i++){

        int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
        double* xsiC = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();

        int *connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();

        int patchC = elementsCoarse_[elCoarse] -> getPatch();
        int *incC = nodesCoarse_[connecCoarse[LNNC-1]] -> getINC();
        double wpcC[LNNC],phiC_[LNNC];
        for (int k = 0; k < LNNC; k++) wpcC[k] = nodesCoarse_[connecCoarse[k]] -> getWeightPC();
        shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC,IsoParCoarse,patchC);

        double u[DIM] = {};
        for (int j=0; j<LNNC; j++){
            for (int k = 0 ; k < DIM; k++){
                u[k] += nodesCoarse_[connecCoarse[j]] -> getVelocity(k) * phiC_[j];
            };
        };

        double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> getWeightFunction();

        for (int k = 0; k < DIM; k++){
            double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> getVelocity(k) * wFunc + u[k] * (1. - wFunc);
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(k,u_int);
        };

    };

    // Printing results
    printResultsLaplace_FEM_ISO(0);

    // PetscFree(dof);
    return 0;
};




template <int DIM>
int Arlequin<DIM>::solveArlequinProblemLaplace_FEM_FEM(int iterNumber, double tolerance)
{
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::vector<int> dofTemp;
    std::vector<double>dofValue;
    setDirichletConstrain_FEM_FEM(dofTemp,dofValue);
    PetscMalloc1(dofTemp.size(), &dof);
    for (size_t i = 0; i < dofTemp.size(); i++)
    {
        dof[i] = dofTemp[i];
    };


    std::clock_t t1 = std::clock();

    int sysSize = NCNumberNodesC + NCNumberNodesF + NCnumNodesGlueZoneFine;

    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        sysSize,sysSize,10000,NULL,10000,NULL,&A);CHKERRQ(ierr);

    for (int i=0; i<sysSize; i++){
        double valu = 1.e-20;
        ierr = MatSetValues(A,1,&i,1,&i,&valu,ADD_VALUES);
    }

    // Divides the matrix between the processes
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
    ierr = VecSetSizes(b,PETSC_DECIDE,sysSize);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&b2);CHKERRQ(ierr);
    ierr = VecSetSizes(b2,PETSC_DECIDE,sysSize);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b2);CHKERRQ(ierr);
       
    ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&All);CHKERRQ(ierr);

    for (int i = 0; i < dofTemp.size(); i++){

        int ind = dofTemp[i];
        ierr = VecSetValues(b2, 1, &ind, &dofValue[i],
                            INSERT_VALUES);
    };
   
    ierr = VecAssemblyBegin(b2);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b2);CHKERRQ(ierr);
  
    //Matrix and vectors - COARSE MESH - IGA mesh
    setMatVecValuesCoarseLaplace_FEM();
    //Matrix and vectors - Lagrange multiplieres - COARSE MESH
    setMatVecValuesLagrangeCoarseLaplace_FEM_FEM();

    // // //Matrix and vectors -FINE MESH  - FEM mesh
    setMatVecValuesFineLaplace_FEM();
    // // // //Matrix and vectors - Lagrange multiplieres - FINE MESH
    setMatVecValuesLagrangeFineLaplace_FEM_FEM();

    //Assemble matrices and vectors
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

    MatZeroRows(A,dofTemp.size(),dof,1.0,b2,b);

    // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

    //Solve using GMRES
    // ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,500); CHKERRQ(ierr);
    // ierr = KSPGetPC(ksp,&pc);
    // ierr = PCSetType(pc,PCNONE);
    // ierr = KSPSetType(ksp,KSPDGMRES); CHKERRQ(ierr);
    // ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);

    // //Solve using Mumps
    #if defined(PETSC_HAVE_MUMPS)
    ierr = KSPSetType(ksp,KSPPREONLY);
    ierr = KSPGetPC(ksp,&pc);
    ierr = PCSetType(pc, PCLU);
    #endif
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);

    //Solve Linear System
    ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
    ierr = KSPGetTotalIterations(ksp, &iterations);

    // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);

    // //Updates nodal values

    double u_;
    double lag_;
    double normU = 0.;
    double normL = 0.;
    Ione = 1;

    for (int i = 0; i < numNodesCoarse; ++i){
        Ii = i;
        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
        u_ = val;
        nodesCoarse_[i] -> setVelocityComponent(0,u_);
        normU += val*val;
    };

    for (int i = 0; i < numNodesFine; ++i){
        Ii = i + NCNumberNodesC;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        u_ = val;
        nodesFine_[i] -> setVelocityComponent(0,u_);
        normU += val*val;
    };

    for (int i = 0; i < NCnumNodesGlueZoneFine; ++i){
        Ii = NCNumberNodesF + NCNumberNodesC + i;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        lag_ = val;
        normL += val * val;
        nodesFine_[nodesGlueZoneFine_[i]] -> setLagrangeMultiplier(0,lag_);
    };

    std::clock_t t2 = std::clock();

    if(rank == 0){

        std::cout   << "   Du Norm = " << std::scientific << sqrt(normU)
                    << " " << sqrt(normL)
                    << "  Time (s) = " << std::fixed
                    << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
    };

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&u); CHKERRQ(ierr);
    ierr = VecDestroy(&All); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);


    //Computes the real velocity
    QuadShapeFunction<DIM> shapeQuad;

    for (int i = 0; i < numNodesFine; i++){
        for (int k = 0; k < DIM; k++) nodesFine_[i] -> setVelocityArlequin(k,nodesFine_[i] -> getVelocity(k));
    };

    int LNNC = 4*DIM-2;

    for (int i = 0; i<numNodesGlueZoneFine; i++){

        int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
        double* xsiC = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();

        int *connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();

        double phiC_[LNNC];
        shapeQuad.evaluateFem(xsiC,phiC_);

        double u[DIM] = {};
        for (int j=0; j<LNNC; j++)
            for (int k = 0 ; k < DIM; k++)
                u[k] += nodesCoarse_[connecCoarse[j]] -> getVelocity(k) * phiC_[j];
       

        double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> getWeightFunction();

        for (int k = 0; k < DIM; k++){
            double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> getVelocity(k) * wFunc + u[k] * (1. - wFunc);
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(k,u_int);
        };

    };

    // Printing results
    printResults_FEM_FEM(0);

    PetscFree(dof);
    return 0;
};



template <int DIM>
int Arlequin<DIM>::solveArlequinProblemLaplace_ISO_ISO(int iterNumber, double tolerance)
{
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::vector<int> dofTemp;
    std::vector<double>dofValue;
    setDirichletConstrain_ISO_ISO(dofTemp,dofValue);
    PetscMalloc1(dofTemp.size(), &dof);
    for (size_t i = 0; i < dofTemp.size(); i++)
    {
        dof[i] = dofTemp[i];
    };   

    std::clock_t t1 = std::clock();

    int sysSize = NCNumberNodesC + NCNumberNodesF + NCnumNodesGlueZoneFine;

    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        sysSize,sysSize,10000,NULL,10000,NULL,&A);CHKERRQ(ierr);

    for (int i=0; i<sysSize; i++){
        double valu = 1.e-20;
        ierr = MatSetValues(A,1,&i,1,&i,&valu,ADD_VALUES);
    }

    // Divides the matrix between the processes
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
    ierr = VecSetSizes(b,PETSC_DECIDE,sysSize);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&b2);CHKERRQ(ierr);
    ierr = VecSetSizes(b2,PETSC_DECIDE,sysSize);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b2);CHKERRQ(ierr);

    ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&All);CHKERRQ(ierr);

    for (int i = 0; i < dofTemp.size(); i++){

        int ind = dofTemp[i];
        ierr = VecSetValues(b2, 1, &ind, &dofValue[i],
                            INSERT_VALUES);
    };
   
    ierr = VecAssemblyBegin(b2);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b2);CHKERRQ(ierr);

    //Matrix and vectors - COARSE MESH - IGA mesh
    setMatVecValuesCoarseLaplace_ISO();
    // //Matrix and vectors - Lagrange multiplieres - COARSE MESH
    setMatVecValuesLagrangeCoarseLaplace_ISO_ISO();

    //Matrix and vectors -FINE MESH  - FEM mesh
    setMatVecValuesFineLaplace_ISO();
    //Matrix and vectors - Lagrange multiplieres - FINE MESH
    setMatVecValuesLagrangeFineLaplace_ISO_ISO();

    //Assemble matrices and vectors
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

    MatZeroRows(A,dofTemp.size(),dof,1.0,b2,b);

    // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

    //Solve using GMRES
    // ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,500); CHKERRQ(ierr);
    // ierr = KSPGetPC(ksp,&pc);
    // ierr = PCSetType(pc,PCNONE);
    // ierr = KSPSetType(ksp,KSPDGMRES); CHKERRQ(ierr);
    // ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);

    // //Solve using Mumps
    #if defined(PETSC_HAVE_MUMPS)
    ierr = KSPSetType(ksp,KSPPREONLY);
    ierr = KSPGetPC(ksp,&pc);
    ierr = PCSetType(pc, PCLU);
    #endif
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);

    //Solve Linear System
    ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
    ierr = KSPGetTotalIterations(ksp, &iterations);

    // //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // //VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);

    // //Updates nodal values

    double u_;
    double lag_;
    double normU = 0.;
    double normL = 0.;
    Ione = 1;

    for (int i = 0; i < numNodesCoarse; ++i){
        int newconi = nodesCoarse_[i] -> getnewcon();
        Ii = newconi;
        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
        u_ = val;
        nodesCoarse_[i] -> setVelocityComponent(0,u_);
        normU += val*val;
    };

    for (int i = 0; i < numNodesFine; ++i){
        int newconi = nodesFine_[i] -> getnewcon();
        Ii = newconi + NCNumberNodesC;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        u_ = val;
        nodesFine_[i] -> setVelocityComponent(0,u_);
        normU += val*val;
    };

    for (int i = 0; i < NCnumNodesGlueZoneFine; ++i){
        Ii = NCNumberNodesF + NCNumberNodesC + i;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        lag_ = val;
        normL += val * val;
        nodesFine_[nodesGlueZoneFine_[i]] -> setLagrangeMultiplier(0,lag_);
    };

    std::clock_t t2 = std::clock();

    if(rank == 0){

        std::cout   << "   Du Norm = " << std::scientific << sqrt(normU)
                    << " " << sqrt(normL)
                    << "  Time (s) = " << std::fixed
                    << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
    };

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&u); CHKERRQ(ierr);
    ierr = VecDestroy(&All); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);


    // //Computes the real velocity
    QuadShapeFunction<DIM> shapeQuad;

    for (int i = 0; i < numNodesFine; i++){
        for (int k = 0; k < DIM; k++) nodesFine_[i] -> setVelocityArlequin(k,nodesFine_[i] -> getVelocity(k));
    };

    int LNNC = 18*DIM-27;

    for (int i = 0; i<numNodesGlueZoneFine; i++){

        int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
        double* xsiC = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();

        int *connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();

        int patchC = elementsCoarse_[elCoarse] -> getPatch();
        int *incC = nodesCoarse_[connecCoarse[LNNC-1]] -> getINC();
        double wpcC[LNNC],phiC_[LNNC];
        for (int k = 0; k < LNNC; k++) wpcC[k] = nodesCoarse_[connecCoarse[k]] -> getWeightPC();
        shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC,IsoParCoarse,patchC);

        double u[DIM] = {};
        for (int j=0; j<LNNC; j++){
            for (int k = 0 ; k < DIM; k++){
                u[k] += nodesCoarse_[connecCoarse[j]] -> getVelocity(k) * phiC_[j];
            };
        };

        double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> getWeightFunction();

        for (int k = 0; k < DIM; k++){
            double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> getVelocity(k) * wFunc + u[k] * (1. - wFunc);
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(k,u_int);
        };

    };

    // Printing results
    printResultsLaplace_ISO_ISO(0);

    return 0;
};




template class Arlequin<2>;
template class Arlequin<3>;