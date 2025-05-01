#ifndef KGRAPH_H
#define KGRAPH_H


#include <vector>
#include <fstream>


#include "base_def.h"
#include "GraphUtil.h"


namespace math {

    class KGraph {

        friend class kGraphFactory;

        u_int* KAO;     //массив сумм степеней вершин (нумерация вершин идёт с 1). KAO[0] = 0. deg(i) = KAO[i] - KAO[i - 1], i = 1, ..., n
        u_int* FO;      //массив, где элементы от FO[KAO[i-1]] до FO[KAO[i]] - вершины, смежные с вершиной i
        double* FORel;  //массив вероятностей присутствия рёбер в графе. Вероятность FORel[j] соответствует ребру FO[j] (по умолчанию равна 0.75) 
        bool* targets;  // массив целевых вершин

        u_int vertNumb; //количество вершин в графе
        u_int edgNumb;  //количество рёбер в графе

        u_int Boolka(u_int i, u_int cutv1, u_int cutv2);

    private://Graph.cpp

        void memclear();

        void memcopy(const KGraph& graph);

        void memallocNewGraph(u_int _vertNumb, u_int _edgNumb, u_int*& newKAO, u_int*& newFO, double*& newFORel, bool*& Targets) const;

        void memupdate(u_int vertLim, u_int edgLim, u_int*& _KAO, u_int*& _FO, double*& _FORel, bool*& Targets);

        void memrebase(u_int vertLim, u_int edgLim, u_int*& _KAO, u_int*& _FO, double*& _FORel, bool*& Targets);

    public:

        KGraph();

        KGraph(const KGraph& kGraph);

        KGraph(u_int _vertNumb, u_int _edgNumb, u_int*& _KAO, u_int*& _FO, bool*& targets, double*& _FORel);

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

        KGraph* deleteEdgeK(u_int u, u_int v);

        int searchEdge(u_int u, u_int v);

        void KParallelSeriesTransformation(double& p);

        KGraph* KParallelSeriesTransformation(KGraph*& G, double& p);

        void Transformation(u_int u, u_int v, u_int w, double& p);

        KGraph* Transformation(KGraph*& G, u_int u, u_int v, u_int w, double& p);

        u_int FindIndexEdgeForVertex(u_int u, u_int v);

        u_int FindIndexEdgeForVertex(KGraph* G, u_int u, u_int v);

        double baseProbabilities();

        u_int getPolusNumb();

        //KGraph* MergeVertex(KGraph* G, u_int cutv1, u_int cutv2);
    };
}
#endif // !KGRAPH_H



