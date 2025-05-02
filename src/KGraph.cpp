
#include "../headers/KGraph.h"
#include <iostream>

namespace math {

    void KGraph::memclear() {
        delete[] FO;
        delete[] KAO;
        delete[] FORel;
        delete[] targets;

        KAO = nullptr;
        FO = nullptr;
        FORel = nullptr;
        targets = nullptr;
        vertNumb = 0;
        edgNumb = 0;
    }

    KGraph::KGraph() {
        KAO = nullptr;
        FO = nullptr;
        FORel = nullptr;
        targets = nullptr;
        vertNumb = 0;
        edgNumb = 0;
    }

    KGraph::KGraph(const KGraph& kGraph)
    {
        memcopy(kGraph);

    }

    KGraph::KGraph(u_int _vertNumb, u_int _edgNumb, u_int*& _KAO, u_int*& _FO, bool*& _targets, double*& _FORel)
    {
        vertNumb = _vertNumb;
        edgNumb = _edgNumb;
        KAO = _KAO;
        FO = _FO;
        targets = _targets;
        FORel = _FORel;
    }

    KGraph::KGraph(u_int _vertNumb, u_int _edgNumb) {
        vertNumb = _vertNumb;
        edgNumb = _edgNumb;
        memallocNewGraph(vertNumb, edgNumb, KAO, FO, FORel, targets);
    }


#if __cplusplus >= 201103L

    void KGraph::memmove(KGraph&& other) {
        memclear();

        vertNumb = other.vertNumb;
        edgNumb = other.edgNumb;

        KAO = other.KAO;
        FO = other.FO;
        FORel = other.FORel;
        targets = other.targets;

        other.KAO = nullptr;
        other.FO = nullptr;
        other.FORel = nullptr;
        other.target = nullptr;
        other.vertNumb = 0;
        other.edgNumb = 0;
    }



    KGraph::KGraph(KGraph&& kGraph)
    {

    }
    KGraph& KGraph::operator=(KGraph&& kGraph)
    {

    }
#endif


    KGraph::~KGraph()
    {
        memclear();
    }

    KGraph& KGraph::operator=(const KGraph& kGraph)
    {
        if (this == &kGraph) {
            return *this;
        }

        if (edgNumb != 0 && vertNumb != 0) {
            memclear();
        }

        memcopy(kGraph);

        return *this;
    }


    KGraph& KGraph::operator=(const KGraph* kGraph) {

        if (this == kGraph) {
            return *this;
        }

        if (edgNumb != 0 && vertNumb != 0) {
            memclear();
        }

        memcopy(*kGraph);

        return *this;

    }


    void KGraph::memcopy(const KGraph& kGraph) {
        vertNumb = kGraph.vertNumb;
        edgNumb = kGraph.edgNumb;

        KAO = new u_int[vertNumb + 1];
        FO = new u_int[2 * edgNumb];
        FORel = new double[2 * edgNumb];
        targets = new bool[vertNumb];

        for (u_int i = 0; i <= vertNumb; i++) {
            KAO[i] = kGraph.KAO[i];
        }
        for (u_int i = 0; i < (2 * edgNumb); i++) {
            FO[i] = kGraph.FO[i];
            FORel[i] = kGraph.FORel[i];
        }

        for (u_int i = 0; i < vertNumb; i++) {
            targets[i] = kGraph.targets[i];
        }

    }


    std::ifstream& operator>>(std::ifstream& input, KGraph& kGraph)
    {

        input >> kGraph.vertNumb;
        input >> kGraph.edgNumb;


        kGraph.KAO = new u_int[kGraph.vertNumb + 1];
        kGraph.FO = new u_int[2 * kGraph.edgNumb];
        kGraph.FORel = new double[2 * kGraph.edgNumb];
        kGraph.targets = new bool[kGraph.vertNumb + 1];

        for (u_int i = 0; i < kGraph.vertNumb + 1; i++) {
            input >> kGraph.KAO[i];
        }

        for (u_int i = 0; i < 2*kGraph.edgNumb; i++) {
            input >> kGraph.FO[i];
            kGraph.FORel[i] = 0.75;
        }

        kGraph.targets[0] = 0;
        for (u_int i = 1; i < kGraph.vertNumb + 1; i++) {
            input >> kGraph.targets[i];
        }

        return input;

    }



    void KGraph::memallocNewGraph(u_int _vertNumb, u_int _edgNumb, u_int*& newKAO, u_int*& newFO, double*& newFORel, bool*& Targets) const {

        newKAO = new u_int[_vertNumb + 1];
        newFO = new u_int[2 * _edgNumb];
        newFORel = new double[2 * _edgNumb];
        Targets = new bool[_vertNumb + 1];
    }



    u_int KGraph::getVertNumb() {
        return vertNumb;
    }
    u_int KGraph::getEdgeNumb() {
        return edgNumb;
    }
    u_int* KGraph::getFO() {
        return FO;
    }

    u_int* KGraph::getKAO() {
        return KAO;
    }

    void KGraph::memupdate(u_int newVertNumb, u_int newEdgeNumb, u_int*& _KAO, u_int*& _FO, double*& _FORel, bool*& Targets) {

        if (newVertNumb > vertNumb) {
            return;             //raise exception here
        }
        if (newEdgeNumb > edgNumb) {
            return;             //raise exception here
        }

        for (u_int i = 0; i < newVertNumb + 1; i++) {
            KAO[i] = _KAO[i];
            targets[i] = Targets[i];
        }

        for (u_int i = 0; i < (2 * newEdgeNumb); i++) {
            FO[i] = _FO[i];
            FORel[i] = _FORel[i];
        }


    }


    void KGraph::memrebase(u_int vertLim, u_int edgLim, u_int*& _KAO, u_int*& _FO, double*& _FORel, bool*& Targets) {
        memclear();
        vertNumb = vertLim;
        edgNumb = edgLim;
        KAO = _KAO;
        FO = _FO;
        FORel = _FORel;
        targets = Targets;
        //memallocNewGraph(vertLim, edgLim, KAO, FO, FORel);
        //memupdate(vertLim, edgLim, _KAO, _FO, _FORel);
    }


    bool* KGraph::getTarget() {
        return targets;
    }

    void KGraph::print()
    {
        if (this == nullptr) return;

        std::cout << vertNumb << std::endl;
        std::cout << edgNumb << std::endl;

        for (u_int i = 0; i < vertNumb + 1; i++) {
            std::cout << KAO[i] << " ";
        }
        std::cout << "\n";
        for (u_int i = 0; i < 2 * edgNumb; i++) {
            std::cout << FO[i] << " ";
        }
        std::cout << "\n";
        for (u_int i = 1; i < vertNumb + 1; i++) {
            std::cout << targets[i] << " ";
        }


        std::cout << "\n";
    }


    u_int KGraph::getPolusNumb () {

        u_int polus_numb = 0;

        for (u_int i = 1; i < vertNumb + 1; i++) {
            if (targets[i] == 1) polus_numb++;
        }
        return polus_numb;
    }

    double KGraph::baseProbabilities() {

        u_int polus_numb = getPolusNumb();

        double Result = 0;

        std::cout << "vertNumb = " << vertNumb << std::endl;
        std::cout << "polusov = " << polus_numb << std::endl;
        std::cout << "edges = " << edgNumb << std::endl;

        if (vertNumb == 1) Result = 1.0;

        else if (vertNumb == 2 && edgNumb == 1) {
            std::cout << "vertNumb = 2\n";
            if (polus_numb == 1) Result = 1.0;
            if (polus_numb == 2) Result = FORel[0];

        }

        else if (vertNumb == 3) {

            if (polus_numb == 2) {
                if (edgNumb == 2) {

                    for (u_int i = 0; i < vertNumb + 1; i++) {
                        if (targets[i] == 1) {

                            for (u_int j = 0; j < KAO[i] - KAO[i - 1]; j++) {

                                if (targets[j] == 1) {
                                    if (isEdge(i, j)) Result = FORel[KAO[i - 1] + j];
                                    else Result = 0;
                                    return Result;
                                }

                            }
                        }
                    }
                }

                if (edgNumb == 3) {
                    double p1, p2, p3;
                    u_int u, v, w;

                    for (u_int i = 0; i < vertNumb + 1; i++) {
                        if (targets[i] == 1) {
                            u = i;
                            for (u_int j = 0; j < KAO[i] - KAO[i - 1]; j++) {

                                if (targets[j] == 1) {
                                    p1 = FORel[KAO[i - 1] + j - 1];
                                    v = j;
                                }

                            }
                        }
                        else {
                            w = i;
                        }
                    }

                    std::cout << "u = " << u << ", v = " << v << ", w = " << w << std::endl;

                    p2 = FORel[FO[KAO[w - 1]]];
                    p3 = FORel[FO[KAO[w - 1] + 1]];

                    Result = p1 + p2 * p3 - p1 * p2 * p3;
                    return Result;
                }

            }

            if (polus_numb == 3) {
                double p1, p2, p3;
                u_int u = 1;
                u_int v = 2;

                p1 = FORel[0];
                p2 = FORel[1];

                for (u_int i = 0; i < 2 * edgNumb; i++) {
                    if (FO[i] != u && FO[i] != v)
                        p3 = FORel[i];
                }

                
                Result = p1 * p2 + p2 * p3 + p1 * p3 - 2 * p1 * p2 * p3;
            }

        }


        return Result;
    }

}





