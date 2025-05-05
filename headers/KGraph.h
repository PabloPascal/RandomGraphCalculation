#ifndef KGRAPH_H
#define KGRAPH_H


#include <vector>
#include <fstream>


#include "base_def.h"
#include "GraphUtil.h"


namespace math {


    u_int FindIndexEdgeForVertex(KGraph* G, u_int u, u_int v);


    class KGraph {

        friend class kGraphFactory;
        friend u_int FindIndexEdgeForVertex(KGraph* G, u_int u, u_int v);


        u_int* KAO;     //������ ���� �������� ������ (��������� ������ ��� � 1). KAO[0] = 0. deg(i) = KAO[i] - KAO[i - 1], i = 1, ..., n
        u_int* FO;      //������, ��� �������� �� FO[KAO[i-1]] �� FO[KAO[i]] - �������, ������� � �������� i
        double* FORel;  //������ ������������ ����������� ���� � �����. ����������� FORel[j] ������������� ����� FO[j] (�� ��������� ����� 0.75) 
        bool* targets;  // ������ ������� ������

        u_int vertNumb; //���������� ������ � �����
        u_int edgNumb;  //���������� ���� � �����

        u_int Boolka(u_int i, u_int cutv1, u_int cutv2);

    private://Graph.cpp

        void memclear();

        void memcopy(const KGraph& graph);

        void memallocNewGraph(u_int _vertNumb, u_int _edgNumb, u_int*& newKAO, u_int*& newFO, double*& newFORel, bool*& Targets) const;

        void memupdate(u_int newVertNumb, u_int newEdgeNumb, u_int*& _KAO, u_int*& _FO, double*& _FORel, bool*& Targets);

        void memrebase(u_int vertLim, u_int edgLim, u_int*& _KAO, u_int*& _FO, double*& _FORel, bool*& Targets);

    public:

        KGraph();

        KGraph(const KGraph& kGraph);

        KGraph(u_int _vertNumb, u_int _edgNumb, u_int*& _KAO, u_int*& _FO, bool*& _targets, double*& _FORel);

        KGraph(u_int _vertNumb, u_int _edgNumb);

#if __cplusplus >= 201103L
        void memmove(KGraph&& other);

        KGraph(KGraph&& kGraph);

        KGraph& operator=(KGraph&& kGraph);
#endif

        ~KGraph();

        friend std::ifstream& operator >>(std::ifstream& input, KGraph& kGraph);

        KGraph& operator=(const KGraph& kGraph);

        KGraph& operator=(const KGraph* kGraph);

        
        //support func
        u_int getVertNumb();
        u_int getEdgeNumb();
        u_int* getKAO();
        u_int* getFO();
        bool* getTarget();

        void print();

    public: //kGraphOperation.cpp

        KGraph* MergeVertex(u_int u, u_int v);

        bool Kconnective();

        bool isEdge(u_int u, u_int v) const;

        bool isEdgeWithPos(u_int u, u_int v, u_int& pos) const;

        KGraph* deleteEdge(u_int u, u_int v);

        int searchEdge(u_int u, u_int v);

        void KParallelSeriesTransformation(double& p);
        void Transformation(u_int u, u_int v, u_int w, double& p);


        double baseProbabilities();

        u_int getPolusNumb();

        //KGraph* MergeVertex(KGraph* G, u_int cutv1, u_int cutv2);
    };


}
#endif // !KGRAPH_H



