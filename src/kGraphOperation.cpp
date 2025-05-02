
#include "../headers/KGraph.h"
#include <iostream>

namespace math {

    bool KGraph::Kconnective() {

        if (this == nullptr) {
            return false;
        }


        int sum = 1,sum1 = 0, l = 1, LB;
        std::vector<int> A, B, Spot;

        bool Result;

        sum = 1;
        A.resize(vertNumb + 1);
        Spot.resize(vertNumb + 1);
        B.resize(vertNumb + 1);

        for (int i = 1; i < vertNumb + 1; i++) {
            if (targets[i] == 1) {
                A[0] = i;
                break;
            }

        }


        for (int i = 1; i < vertNumb + 1; i++) {
            Spot[i] = 0;
        }

        Spot[A[0]] = 1;

        while (l > 0) {
            LB = 0;
            for (int i = 0; i <  l; i++) {
                for (int j = KAO[A[i] - 1]; j < KAO[A[i]]; j++) {//!for (u_int j{ G.KAO[A[i] - 1] }; j < G.KAO[A[i]]; j++) {


                    if (Spot[FO[j]] == 0) {
                        LB++;
                        //B.resize(B.size() + 1);
                        B[LB - 1] = FO[j];
                        Spot[FO[j]] = 1;

                        if (targets[FO[j]] == 1) sum++;
                        
                    }
                }

            }//for int i

            l = LB;

            if (l > 0) {
                for (int i = 0; i < l; i++) {
                    A[i] = B[i];
                }
            }

        }//while

        sum1 = 0;
        for (int i = 1; i < vertNumb + 1; i++) {
            if (targets[i] == 1) sum1++;
            if (sum == sum1)  Result = true;
            else Result = false;
            if (vertNumb == 1) Result = true;
        }


        return Result;
    }//KConnective()
    

    KGraph* KGraph::deleteEdgeK(u_int u, u_int v) {

        //std::cout << "delete Edge\n";

        if (u == v) return nullptr; 

        //if ((u == 0) || (u > edgNumb) || (v == 0) || (v > edgNumb) || (!isEdge(u, v))) {
        //    return nullptr;     //raise exception here
        //}

        u_int newVertNumb = vertNumb;
        u_int newEdgNumb = edgNumb - 1;
        u_int* newKAO;
        u_int* newFO;
        double* newFORel;
        bool* Targets;

        memallocNewGraph(newVertNumb, newEdgNumb, newKAO, newFO, newFORel, Targets);

        newKAO[0] = 0;
        u_int iC = 1;
        u_int jC = 0;
        for (u_int i = 1; i <= vertNumb; i++) {
            newKAO[iC] = newKAO[iC - 1] + (KAO[i] - KAO[i - 1]);
            if ((i == u) || (i == v)) {
                for (u_int j = KAO[i - 1]; j < KAO[i]; j++) {
                    if (((i == u) && (FO[j] == v)) || ((i == v) && (FO[j] == u))) {
                        newKAO[iC]--;
                        continue;
                    }
                    else {
                        newFO[jC] = FO[j];
                        newFORel[jC] = FORel[j];
                        jC++;
                    }
                }
            }
            else {
                for (u_int j = KAO[i - 1]; j < KAO[i]; j++) {
                    newFO[jC] = FO[j];
                    newFORel[jC] = FORel[j];
                    jC++;
                }
            }

            iC++;
        }

        KGraph* result = new KGraph(newVertNumb, newEdgNumb, newKAO, newFO, targets, newFORel);
        return result;
    }



    int KGraph::searchEdge(u_int i, u_int j) {
        int result = 2 * edgNumb;

        for (int k = KAO[i - 1]; k < KAO[i]; i++) {
            if (FO[k] == j) return k;
        }

        return result;
    }




    bool KGraph::isEdge(u_int u, u_int v) const {
        if (u == v) {
            return false;
        }

        for (u_int i = KAO[u - 1]; i < KAO[u]; i++) {
            if (FO[i] == v) {
                return true;
            }
        }
        return false;
    }


    u_int KGraph::Boolka(u_int i, u_int cutv1, u_int cutv2) {
        u_int j, sum;

        if (i == cutv1) return 1;
        else if (i == cutv2) return 2;
        else {
            sum = 0;
            for (j = KAO[i - 1]; j < KAO[i] - 1; j++) {
                if (FO[j] == cutv1 || FO[j] == cutv2) sum++;
                if (sum == 2) return 3;
                else return 4;
            }


        }

    }

    //NV - numb vert, NE - numb edge

    KGraph* KGraph::MergeVertex(u_int u, u_int v) {

        //if (this == nullptr) return nullptr;

        //std::cout << "merge\n";

        u_int newVertNumb = vertNumb - 1;
        u_int newEdgNumb = edgNumb - 1;
        u_int* newKAO;
        u_int* newFO;
        double* newFORel;
        bool* Targets;

        
        u_int newVert;      //номер вершины, в которую осуществляется стягивание (выбирается наименьшая по номеру из пары (u, v))
        u_int delVert;      //номер вершины, которая будет удалена из графа

        if (u < v) {
            newVert = u;
            delVert = v;
        }
        else if (v < u) {
            newVert = v;
            delVert = u;
        }
        else {
            return nullptr;         //raise exception here
        }

        for (u_int i = KAO[newVert - 1]; i < KAO[newVert]; i++) {
            if (FO[i] == v) {
                continue;
            }
            if (isEdge(FO[i], delVert)) {
                newEdgNumb--;
            }
        }

        memallocNewGraph(newVertNumb, newEdgNumb, newKAO, newFO, newFORel, Targets);

        newKAO[0] = 0;
        u_int iC = 1;
        u_int jC = 0;
        u_int pos = 0;
        for (u_int i = 1; i <= vertNumb; i++) {
            if (i == delVert) {
                continue;
            }
            else if (i == newVert) {
                newKAO[iC] = newKAO[iC - 1] + (KAO[i] - KAO[i - 1]) + (KAO[delVert] - KAO[delVert - 1]) - 2;
                for (u_int j = KAO[delVert - 1]; j < KAO[delVert]; j++) {
                    if (FO[j] == newVert) {
                        continue;
                    }
                    if (!isEdge(FO[j], newVert)) {
                        newFO[jC] = FO[j];
                        newFORel[jC] = FORel[j];
                        if (FO[j] > delVert) {
                            newFO[jC]--;
                        }
                        jC++;
                    }
                    else {
                        continue;
                    }
                }
                for (u_int j = KAO[i - 1]; j < KAO[i]; j++) {
                    if (FO[j] == delVert) {
                        continue;
                    }
                    if (isEdgeWithPos(FO[j], delVert, pos)) {
                        newKAO[iC]--;
                        newFORel[jC] = FORel[j] + FORel[pos] - FORel[j] * FORel[pos];
                    }
                    else {
                        newFORel[jC] = FORel[j];
                    }
                    newFO[jC] = FO[j];
                    if (FO[j] > delVert) {
                        newFO[jC]--;
                    }
                    jC++;
                }
            }
            else {
                u_int posNew = 0;
                u_int posDel = 0;
                newKAO[iC] = newKAO[iC - 1] + (KAO[i] - KAO[i - 1]);
                bool newVertConnected = isEdgeWithPos(i, newVert, posNew);
                bool delVertConnected = isEdgeWithPos(i, delVert, posDel);

                if (!delVertConnected) {
                    for (u_int j = KAO[i - 1]; j < KAO[i]; j++) {
                        newFO[jC] = FO[j];
                        newFORel[jC] = FORel[j];
                        if (FO[j] > delVert) {
                            newFO[jC]--;
                        }
                        jC++;
                    }
                }
                else if (!newVertConnected && delVertConnected) {
                    for (u_int j = KAO[i - 1]; j < KAO[i]; j++) {
                        if (FO[j] == delVert) {
                            newFO[jC] = newVert;
                        }
                        else {
                            newFO[jC] = FO[j];
                            if (FO[j] > delVert) {
                                newFO[jC]--;
                            }
                        }
                        newFORel[jC] = FORel[j];
                        jC++;
                    }
                }
                else if (newVertConnected && delVertConnected) {
                    newKAO[iC]--;
                    for (u_int j = KAO[i - 1]; j < KAO[i]; j++) {
                        if (FO[j] == delVert) {
                            continue;
                        }
                        newFO[jC] = FO[j];
                        if (FO[j] == newVert) {
                            newFORel[jC] = FORel[j] + FORel[posDel] - FORel[j] * FORel[posDel];
                        }
                        else {
                            newFORel[jC] = FORel[j];
                        }
                        if (FO[j] > delVert) {
                            newFO[jC]--;
                        }
                        jC++;
                    }
                }
            }

            iC++;
        }


        Targets[0] = 0;
        for (u_int i{ 1 }; i < v; i++) {
            Targets[i] = targets[i];
        }
        for (u_int i{ v + 1 }; i < vertNumb + 1; i++) {
            Targets[i - 1] = targets[i];
        }
        if (targets[u] == 1 || targets[v] == 1) Targets[u] = 1;
        KGraph* result = new KGraph(newVertNumb, newEdgNumb, newKAO, newFO,Targets, newFORel);
        return result;
    }



    bool KGraph::isEdgeWithPos(u_int u, u_int v, u_int& pos) const {
        if (u == v) {
            return false;
        }

        for (u_int i = KAO[u - 1]; i < KAO[u]; i++) {
            if (FO[i] == v) {
                pos = i;
                return true;
            }
        }
        return false;
    }





    KGraph* KParallelSeriesTransformation(KGraph*& G, double& p)
    {
        bool b = false;
        u_int  u = 0, v = 0, w, e1, e2;
        double p1, p2, p3;
        p = 1;

        u_int i;
        for (i = 1; i < G->vertNumb + 1; i++) {
            if ((G->KAO[i] - G->KAO[i - 1] == 2) && ((G->targets[i] == 0) ||
                ((G->targets[i] == 1) && (G->targets[G->FO[G->KAO[i - 1]]] == 1) &&
                    (G->targets[G->FO[G->KAO[i] - 1]] == 1)))) {
                b = true;
                v = i;
                u = G->FO[G->KAO[i - 1]];
                w = G->FO[G->KAO[i] - 1];
                p1 = G->FORel[G->KAO[i - 1]];
                p2 = G->FORel[G->KAO[i] - 1];
                break;
            }
        }

        while (b) {
            if (G->targets[i] == 1) {
                p *= (p1 + p2 - p1 * p2);
                p3 = p1 * p2 / (p1 + p2 - p1 * p2); //rc
            }
            else p3 = p1 * p2;
            G = Transformation(G, u, v, w, p3);
            b = false;
            for (i = 1; i < G->vertNumb + 1; i++) {
                if (G->KAO[i] - G->KAO[i - 1] == 2 && ((G->targets[i] == 0) ||
                    ((G->targets[i] == 1) && (G->targets[G->FO[G->KAO[i - 1]]] == 1) &&
                        (G->targets[G->FO[G->KAO[i] - 1]] == 1)))) {
                    b = true;
                    v = i;
                    u = G->FO[G->KAO[i - 1]];
                    w = G->FO[G->KAO[i] - 1];
                    p1 = G->FORel[G->KAO[i - 1]];
                    p2 = G->FORel[G->KAO[i] - 1];
                    break;
                }
            }
        }
        return G;
    }




    KGraph* Transformation(KGraph*& G, u_int u, u_int v, u_int w, double& p)
    {
        int l = 0;

        std::vector<u_int> Numbers;
        Numbers.resize(G->vertNumb + 1);

        KGraph* Result = new KGraph(G->vertNumb - 1, G->edgNumb);

        Result->targets[0] = 0;
        for (u_int i{ 1 }; i < v; i++) {
            Result->targets[i] = G->targets[i];
        }
        for (u_int i{ v }; i < Result->vertNumb + 1; i++) {
            Result->targets[i] = G->targets[i + 1];
        }

        Result->KAO[0] = 0;
        for (u_int i{ 1 }; i < v; i++) {
            Numbers[i] = i;
        }
        Numbers[v] = 0;
        for (u_int i{ v + 1 }; i < G->vertNumb + 1; i++) {
            Numbers[i] = i - 1;
        }

        int SE = FindIndexEdgeForVertex(G, u, w);
        bool bol = true;
        if (SE < G->edgNumb * 2) {
            p = 1 - (1 - G->FORel[SE]) * (1 - p);
            Result->edgNumb = G->edgNumb - 2;
            for (u_int i{ 1 }; i < G->vertNumb + 1; i++) {
                if (i == u || i == w) {
                    Result->KAO[Numbers[i]] = Result->KAO[Numbers[i] - 1];
                    for (u_int j{ G->KAO[i - 1] }; j < G->KAO[i]; j++) {
                        if (G->FO[j] == w || G->FO[j] == u) {
                            Result->FO[l] = Numbers[G->FO[j]];
                            Result->FORel[l] = p;
                            Result->KAO[Numbers[i]]++;
                            l++;
                        }
                        else if (G->FO[j] != v) {
                            Result->FO[l] = Numbers[G->FO[j]];
                            Result->FORel[l] = G->FORel[j];;
                            Result->KAO[Numbers[i]]++;
                            l++;
                        }
                    }
                }
                else if (i != v) {
                    Result->KAO[Numbers[i]] = Result->KAO[Numbers[i] - 1];
                    for (u_int j{ G->KAO[i - 1] }; j < G->KAO[i]; j++) {
                        Result->FO[l] = Numbers[G->FO[j]];
                        Result->FORel[l] = G->FORel[j];;
                        Result->KAO[Numbers[i]]++;
                        l++;
                    }
                }
            }
        }
        else {
            Result->edgNumb = G->edgNumb - 1;
            for (int i{ 1 }; i < G->vertNumb + 1; i++) {
                if (i == u) {
                    Result->KAO[Numbers[i]] = Result->KAO[Numbers[i] - 1];
                    Result->FO[l] = Numbers[w];
                    Result->FORel[l] = p;
                    Result->KAO[Numbers[i]]++;
                    l++;
                    for (u_int j{ G->KAO[i - 1] }; j < G->KAO[i]; j++) { //bpv
                        if (G->FO[j] != v) {
                            Result->FO[l] = Numbers[G->FO[j]];
                            Result->FORel[l] = G->FORel[j];
                            Result->KAO[Numbers[i]]++;
                            l++;
                        }
                    }
                }
                else if (i == w) {
                    Result->KAO[Numbers[i]] = Result->KAO[Numbers[i] - 1];
                    Result->FO[l] = Numbers[u];
                    Result->FORel[l] = p;
                    Result->KAO[Numbers[i]]++;
                    l++;
                    for (u_int j{ G->KAO[i - 1] }; j < G->KAO[i]; j++) {//bpv
                        if (G->FO[j] != v) {
                            Result->FO[l] = Numbers[G->FO[j]];
                            Result->FORel[l] = G->FORel[j];
                            Result->KAO[Numbers[i]]++;
                            l++;
                        }
                    }
                }
                else if (i != v) {
                    Result->KAO[Numbers[i]] = Result->KAO[Numbers[i] - 1];
                    for (u_int j{ G->KAO[i - 1] }; j < G->KAO[i]; j++) {
                        Result->FO[l] = Numbers[G->FO[j]];
                        Result->FORel[l] = G->FORel[j];
                        Result->KAO[Numbers[i]]++;
                        l++;
                    }
                }
            }
        }
        //Result.VertexCount = kaoLength-1;
        //Result.EdgeCount = foLength / 2;
        //std::cout <<" Transform p = " << p << std::endl;
        return Result;
    }



    u_int FindIndexEdgeForVertex(KGraph* G, u_int u, u_int v)
    {
        for (u_int i{ G->KAO[u - 1] }; i < G->KAO[u]; i++) {
            if (G->FO[i] == v) {
                return i;
            }
        }
        return G->edgNumb * 2;
    }


}




