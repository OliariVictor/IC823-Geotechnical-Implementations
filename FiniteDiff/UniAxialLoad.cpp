//
// Created by victor on 29/09/2020.
//

// Uniaxial load on poles with constant section
// Winkler-like model

//// N(z) - N(z+dz) - kvbw(z)dz = 0
//// N = -sigmaA = -EeA
//// e = dw/dz
//// EA d^2w/dz^2 - kvbw(z) = 0

//// EA/h^2(w_{i-1}-2w_{i}+w_{i+1})-kvbw_i = 0 on nodes \in [0,n];
//// Insertion of virtual nodes, node = -1 and node = n+1;
//// w_{-1} =
//// w_{n+1} =
//// N_i = -EA(w_{i+1}-w{i-1})

//// Boundary condition
//// N(z=0) = P;
//// N(z = L) = kvbw(L)

#include <stdlib.h>
#include <vector>
#include "pzvec.h"
#include "pzmatrix.h"


struct problemData
{
    //Material Properties
    double kv = 50000;    // Reaction coefficient on pole's side (N/m^3)
    double kp = 50000;    // Reaction coefficient on pole's end  (N/m^3)
    double E = 22000000;      // Young Modulus N/m^2

    double P = 1e5;       // Load(N)

    //Geometric properties;
    double d = -1;   // Diameter
    double b = -1.;  // section perimeter
    double A = -1.;  // Section Area
    double L = 12.;  // Pole length (m)
    int n = 20;           // Number of Elements
    double h = L/n;
};

class GeoElement{
private:
    double fx;
    double fw;
    int fIndex;
public:
    GeoElement() : fx(-1.), fw(-1.), fIndex(-1) {}
    ~GeoElement(){};

    void SetX(double x){
        fx = x;
    }
    double GetX(){
        return fx;
    }
    void SetW(double w){
        fw = w;
    }
    double GetW(){
        return fw;
    }
    void SetIndex(int ind){
        fIndex = ind;
    }
    int GetIndex(int ind){
        return fIndex;
    }
};


class Mesh{
private:
/** @brief List of pointers to finite elements */
    std::vector<GeoElement> *fElementVec;
public:
    Mesh(){
        fElementVec = new std::vector<GeoElement>;
    }

    ~Mesh(){};

    void BuildMesh(problemData &pData){
        int elNum = pData.n+1;
        fElementVec->resize(elNum);
        for(int elInd = 0; elInd < elNum; elInd++){
            GeoElement *gel = new GeoElement();
            double x = elInd*pData.h;
            gel->SetX(x);
            gel->SetIndex(elInd);
        }
    }

    std::vector<GeoElement> *GetElementVec(){
        return fElementVec;
    }

    void LoadSolution(TPZFMatrix<double> *F){
        for(int i = 0; i < fElementVec->size(); i++)
            fElementVec->at(i).SetW(F->GetVal(i,0));
    }

    void PrintSolution(){
        double w = -1000;
        for(int i = 0; i < fElementVec->size(); i++) {
            w = fElementVec->at(i).GetW();
            std::cout << "node = " << i << " ; w = " << w << std::endl;
        }
    }
};

class Solution{
protected:
    TPZFMatrix<double> fK;
    TPZVec<double> fF;
    Mesh *fgmesh;
    int sysSize;
public:
    Solution(Mesh *gmesh){
        std::vector<GeoElement> *elVec = gmesh->GetElementVec();
        fgmesh = gmesh;
        sysSize = elVec->size();
        fK.Resize(sysSize,sysSize);
        fF.resize(sysSize);
        for(int i = 0 ; i < sysSize; i++){
            fF[i] = 0;
            for(int j = 0 ; j < sysSize; j++){
                fK(i,j)=0;
            }
        }
    }
    void Assemble_MonoMaterial(problemData &pData){
        double cisCoeff =  pData.kv*pData.b* pData.h*pData.h/(pData.E*pData.A);
        double pCoeff =  pData.kp*pData.h*pData.A/(pData.E*pData.A);
        for(int i = 1; i < sysSize -1; i++) {
            fK(i, i) = cisCoeff + 2;
            fK(i, i-1) = -1;
            fK(i, i+1) = -1;
        }
        //First node
        fK(0, 0) = cisCoeff+1.;
        fK(0, 1) = -1.;
        //End node
        fK(sysSize-1, sysSize-2) = -1.;
        fK(sysSize-1, sysSize-1) = cisCoeff+pCoeff+1;

        fF[0] = pData.P*pData.h/(pData.E*pData.A);
        fK.Print(std::cout);
    }

    void Assemble_DuoMaterial(problemData &pData1, problemData &pData2, int nodeLimit){
        double cisCoeff1 =  pData1.kv*pData1.b* pData1.h*pData1.h/(pData1.E*pData1.A);
        double cisCoeff2 =  pData2.kv*pData2.b* pData2.h*pData2.h/(pData2.E*pData2.A);
        double pCoeff =  pData2.kp*pData2.h*pData2.A/(pData2.E*pData2.A);

        //First node
        fK(0, 0) = cisCoeff1+1.;
        fK(0, 1) = -1.;

        //Internal nodes
        for(int i = 1; i < nodeLimit; i++) {
            fK(i, i) = cisCoeff1 + 2;
            fK(i, i-1) = -1;
            fK(i, i+1) = -1;
        }

        //Node limit
        double EAh1 = pData1.E*pData1.A/pData1.h;
        double EAh2 = pData2.E*pData2.A/pData2.h;
        double kvbh1 = pData1.kv*pData1.b*pData1.h/2;
        double kvbh2 = pData2.kv*pData2.b*pData2.h/2;

        fK(nodeLimit, nodeLimit-1) = -EAh1;
        fK(nodeLimit, nodeLimit) = kvbh1+EAh1+kvbh2+EAh2;
        fK(nodeLimit, nodeLimit+1) = -EAh2;

        //Internal nodes
        for(int i = nodeLimit+1; i < sysSize-1; i++) {
            fK(i, i) = cisCoeff2 + 2;
            fK(i, i-1) = -1;
            fK(i, i+1) = -1;
        }

        fK(sysSize-1, sysSize-2) = -1.;
        fK(sysSize-1, sysSize-1) = cisCoeff2+pCoeff+1;

        fF[0] = pData1.P*pData1.h/(pData1.E*pData1.A);
        fK.Print(std::cout);
        /*for(int i = 0; i< sysSize ; i++)
            std::cout << std::endl << "ind: " << i << ": " << fF[i] << std::endl;*/
    }

    void SolveSystem(Mesh *gmesh){
        TPZFMatrix<double> *Fvec = new TPZFMatrix<double>;
        Fvec->Resize(sysSize,1);
        for (int i =0; i < sysSize; i++) Fvec->PutVal(i,0,fF[i]); fF.Print(std::cout);
        fK.Solve_LU(Fvec);
        gmesh->LoadSolution(Fvec);
    }
};

void SectionProperties(problemData &pData);

int main(int argc, char *argv[]) {
    problemData pData1,pData2;
    pData1.d = 0.25;
    pData2.d = 0.125;
    SectionProperties(pData1);
    SectionProperties(pData2);
    Mesh *gmesh = new Mesh;
    gmesh->BuildMesh(pData1);
    Solution *sol = new Solution(gmesh);
    sol->Assemble_DuoMaterial(pData1,pData2,10);
    sol->SolveSystem(gmesh);
    gmesh->PrintSolution();
}

void SectionProperties(problemData &pData){
    pData.b = 3.14*pData.d;
    pData.A = 3.14*pData.d*pData.d/4;
}