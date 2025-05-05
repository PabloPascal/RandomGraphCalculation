#include "../headers/GraphUtil.h"

#include "../headers/Graph.h"
#include "../headers/KGraph.h"
#include "../headers/GeantGenerator.h"

#include <iostream>
#include <queue>
#include <stack>

#include <cstdlib>
#include <ctime>


namespace math {
    
    Graph *GraphGenerator::generateCompleteGraph(u_int n) {
    u_int gVertNumb = n;
    u_int gEdgNumb = n * (n - 1) / 2;
    u_int *gKAO = new u_int[gVertNumb + 1];
    u_int *gFO = new u_int[2 * gEdgNumb];
    double *gFORel = new double[2 * gEdgNumb];

    u_int iC = 0;
    gKAO[0] = 0;
    for (u_int i = 1; i <= gVertNumb; i++) {
        gKAO[i] = gKAO[i - 1] + (n - 1);
    }
    for (u_int i = 1; i <= gVertNumb; i++) {
        for (u_int j = 1; j <= gVertNumb; j++) {
            if (i == j) {
                continue;
            }
            gFO[iC] = j;
            gFORel[iC] = 0.75;
            iC++;
        }
    }

    Graph *result = new Graph(gVertNumb, gEdgNumb, gKAO, gFO, gFORel);
    return result;
}

Graph *GraphGenerator::generateGrid(u_int n, u_int m) {
    if (m == 0) {
        m = n;
    }


    u_int gVertNumb = (n + 1) * (m + 1);
    u_int gEdgNumb = n * (m + 1) + m * (n + 1);
    u_int *gKAO = new u_int[gVertNumb + 1];
    u_int *gFO = new u_int[2 * gEdgNumb];
    double *gFORel = new double[2 * gEdgNumb];

    for (u_int i = 0; i < (2 * gEdgNumb); i++) {
        gFORel[i] = 0.75;
    }

    u_int iC = 1;
    u_int jC = -1;
    gKAO[0] = 0;
    for (u_int k = 0; k < (m + 1); k++) {       //Ğ³ĞµĞ½ĞµÑ€Ğ°Ñ†Ğ¸Ñ Ñ€ĞµÑˆĞµÑ‚ĞºĞ¸ Ğ¸Ğ´ĞµÑ‚ Ğ¿Ğ¾ Ğ³Ğ¾Ñ€Ğ¸Ğ·Ğ¾Ğ½Ğ°Ñ‚ĞµĞ»ÑŒĞ½Ñ‹Ğ¼ ÑĞ»Ğ¾ÑĞ¼ ÑĞ»ĞµĞ²Ğ° Ğ½Ğ°Ğ¿Ñ€Ğ°Ğ²Ğ¾
        if (k == 0) {
            for (u_int i = 0; i < (n + 1); i++) {
                if (i == 0) {                   //Ğ·Ğ°Ğ¿Ğ¾Ğ»Ğ½ĞµĞ½Ğ¸Ğµ Ğ»ĞµĞ²Ğ¾Ğ³Ğ¾ Ğ²ĞµÑ€Ñ…Ğ½ĞµĞ³Ğ¾ ÑƒĞ·Ğ»Ğ°
                    gKAO[iC] = gKAO[iC - 1] + 2;
                    gFO[++jC] = iC + 1;
                    gFO[++jC] = iC  + (n + 1);                    
                }
                else if (i == n) {              //Ğ·Ğ°Ğ¿Ğ¾Ğ»Ğ½ĞµĞ½Ğ¸Ğµ Ğ¿Ñ€Ğ°Ğ²Ğ¾Ğ³Ğ¾ Ğ²ĞµÑ€Ñ…Ğ½ĞµĞ³Ğ¾ ÑƒĞ·Ğ»Ğ°
                    gKAO[iC] = gKAO[iC - 1] + 2;
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + (n + 1);
                }
                else {
                    gKAO[iC] = gKAO[iC - 1] + 3;
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + 1;
                    gFO[++jC] = iC + (n + 1);
                }
                iC++;
            }
        }
        else if (k == m) {
            for (u_int i = 0; i < (n + 1); i++) {
                if (i == 0) {                   //Ğ·Ğ°Ğ¿Ğ¾Ğ»Ğ½ĞµĞ½Ğ¸Ğµ Ğ»ĞµĞ²Ğ¾Ğ³Ğ¾ Ğ½Ğ¸Ğ¶Ğ½ĞµĞ³Ğ¾ ÑƒĞ·Ğ»Ğ°
                    gKAO[iC] = gKAO[iC - 1] + 2;
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC + 1; 
                }
                else if (i == n) {              //Ğ·Ğ°Ğ¿Ğ¾Ğ»Ğ½ĞµĞ½Ğ¸Ğµ Ğ¿Ñ€Ğ°Ğ²Ğ¾Ğ³Ğ¾ Ğ½Ğ¸Ğ¶Ğ½ĞµĞ³Ğ¾ ÑƒĞ·Ğ»Ğ°
                    gKAO[iC] = gKAO[iC - 1] + 2;
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - 1;
                }
                else {
                    gKAO[iC] = gKAO[iC - 1] + 3;
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + 1;
                }
                iC++;
            }
        }
        else {
            for (u_int i = 0; i < (n + 1); i++) {
                if (i == 0) {
                    gKAO[iC] = gKAO[iC - 1] + 3;
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC + 1;
                    gFO[++jC] = iC + (n + 1);
                }
                else if (i == n) {
                    gKAO[iC] = gKAO[iC - 1] + 3;
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + (n + 1);
                }
                else {
                    gKAO[iC] = gKAO[iC - 1] + 4;
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + 1;
                    gFO[++jC] = iC + (n + 1);
                }
                iC++;
            }
        }
    }

    Graph *result = new Graph(gVertNumb, gEdgNumb, gKAO, gFO, gFORel);
    return result;
}

Graph *GraphGenerator::generateGridWithDiags(u_int n, u_int m) {
    if (m == 0) {
        m = n;
    }
    
    u_int gVertNumb = (n + 1) * (m + 1);
    u_int gEdgNumb = 4 * n * m + n + m;
    u_int *gKAO = new u_int[gVertNumb + 1];
    u_int *gFO = new u_int[2 * gEdgNumb];
    double *gFORel = new double[2 * gEdgNumb];

    for (u_int i = 0; i < (2 * gEdgNumb); i++) {
        gFORel[i] = 0.75;
    }

    u_int iC = 1;
    u_int jC = -1;
    gKAO[0] = 0;
    for (u_int k = 0; k < (m + 1); k++) {
        if (k == 0) {
            for (u_int i = 0; i < (n + 1); i++) {
                if (i == 0) {                   //Ğ·Ğ°Ğ¿Ğ¾Ğ»Ğ½ĞµĞ½Ğ¸Ğµ Ğ»ĞµĞ²Ğ¾Ğ³Ğ¾ Ğ²ĞµÑ€Ñ…Ğ½ĞµĞ³Ğ¾ ÑƒĞ·Ğ»Ğ°
                    gKAO[iC] = gKAO[iC - 1] + 3;
                    gFO[++jC] = iC + 1;
                    gFO[++jC] = iC + (n + 1);
                    gFO[++jC] = iC + (n + 2);
                }
                else if (i == n) {              //Ğ·Ğ°Ğ¿Ğ¾Ğ»Ğ½ĞµĞ½Ğ¸Ğµ Ğ¿Ñ€Ğ°Ğ²Ğ¾Ğ³Ğ¾ Ğ²ĞµÑ€Ñ…Ğ½ĞµĞ³Ğ¾ ÑƒĞ·Ğ»Ğ°
                    gKAO[iC] = gKAO[iC - 1] + 3;
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + n;
                    gFO[++jC] = iC + (n + 1);
                }
                else {
                    gKAO[iC] = gKAO[iC - 1] + 5;
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + 1;
                    gFO[++jC] = iC + n;
                    gFO[++jC] = iC + (n + 1);
                    gFO[++jC] = iC + (n + 2);
                }
                iC++;
            }
        }
        else if (k == m) {
            for (u_int i = 0; i < (n + 1); i++) {
                if (i == 0) {                   //Ğ·Ğ°Ğ¿Ğ¾Ğ»Ğ½ĞµĞ½Ğ¸Ğµ Ğ»ĞµĞ²Ğ¾Ğ³Ğ¾ Ğ½Ğ¸Ğ¶Ğ½ĞµĞ³Ğ¾ ÑƒĞ·Ğ»Ğ°
                    gKAO[iC] = gKAO[iC - 1] + 3;
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - n;
                    gFO[++jC] = iC + 1;
                }
                else if (i == n) {              //Ğ·Ğ°Ğ¿Ğ¾Ğ»Ğ½ĞµĞ½Ğ¸Ğµ Ğ¿Ñ€Ğ°Ğ²Ğ¾Ğ³Ğ¾ Ğ½Ğ¸Ğ¶Ğ½ĞµĞ³Ğ¾ ÑƒĞ·Ğ»Ğ°
                    gKAO[iC] = gKAO[iC - 1] + 3;
                    gFO[++jC] = iC - (n + 2);
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - 1;
                }
                else {
                    gKAO[iC] = gKAO[iC - 1] + 5;
                    gFO[++jC] = iC - (n + 2);
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - n;
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + 1;
                }
                iC++;
            }
        }
        else {
            for (u_int i = 0; i < (n + 1); i++) {
                if (i == 0) {
                    gKAO[iC] = gKAO[iC - 1] + 5;
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - n;
                    gFO[++jC] = iC + 1;
                    gFO[++jC] = iC + (n + 1);
                    gFO[++jC] = iC + (n + 2); 
                }
                else if (i == n) {
                    gKAO[iC] = gKAO[iC - 1] + 5;
                    gFO[++jC] = iC - (n + 2);
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + n;
                    gFO[++jC] = iC + (n + 1); 
                }
                else {
                    gKAO[iC] = gKAO[iC - 1] + 8;
                    gFO[++jC] = iC - (n + 2);
                    gFO[++jC] = iC - (n + 1);
                    gFO[++jC] = iC - n;
                    gFO[++jC] = iC - 1;
                    gFO[++jC] = iC + 1;
                    gFO[++jC] = iC + n;
                    gFO[++jC] = iC + (n + 1);
                    gFO[++jC] = iC + (n + 2);
                }
                iC++;
            }
        }
    }

    Graph *result = new Graph(gVertNumb, gEdgNumb, gKAO, gFO, gFORel);
    return result;
}

Graph *GraphGenerator::generateGraph(u_int vNumb, u_int eNumb) {
    return nullptr;     //not supported yet
}

Graph *GraphGenerator::generateGEANT() {
    u_int gVertNumb = 103;
    u_int gEdgNumb = 127;
    u_int *gKAO = new u_int[gVertNumb + 1];
    u_int *gFO = new u_int[2 * gEdgNumb];
    double *gFORel = new double[2 * gEdgNumb];

    for (u_int i = 0; i < (2 * gEdgNumb); i++) {
        gFORel[i] = 0.75;
    }

    GEANTGenerator::GEANTGenerate(gKAO, gFO);

    Graph *result = new Graph(gVertNumb, gEdgNumb, gKAO, gFO, gFORel);
    return result;
}


//=====================================================ÂÛÁÎĞ ĞÅÁĞÀ============================================================

int GraphFactory::count_of_triangles = 0;
int GraphFactory::count_of_edges_degs_6 = 0;
int GraphFactory::count_first_triangles = 0;


/*
 ÍÀ×ÀËÜÍÛÉ ÌÅÒÎÄ
*/

/*
void GraphFactory::choseVerts(Graph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {
    vC = graph->vertNumb; //vertex count

    u_int deg = 2 * graph->vertNumb + 2;

    for (u_int i = (2 * graph->edgNumb - 1); i >= graph->KAO[1]; i--) {
        uC = graph->FO[i];
        if (i < graph->KAO[vC - 1]) {
            vC--;
        }
        if (vC > uC) {
            if (graph->KAO[uC] - graph->KAO[uC - 1] + graph->KAO[vC] - graph->KAO[vC - 1] < 7) {
                u = uC;
                v = vC;
                p = graph->FORel[i];
                return;
            }
            else if (graph->KAO[uC] - graph->KAO[uC - 1] + graph->KAO[vC] - graph->KAO[vC - 1] < deg) {
                u = uC;
                v = vC;
                p = graph->FORel[i];
                deg = graph->KAO[uC] - graph->KAO[uC - 1] + graph->KAO[vC] - graph->KAO[vC - 1];
            }
        }
    }
}
*/


/* 
ÂÛÁÈĞÀÅÌ ĞÅÁĞÎ ÏÅĞÂÎÅ Ñ ÊÎÍÖÀ 
*/


/*
void GraphFactory::choseVerts(Graph *&graph, double &p, u_int &u, u_int &v, u_int &uC, u_int &vC) {
     
    u = graph->vertNumb;
    int indexU = graph->KAO[u-1];
    int degU = graph->KAO[u] - graph->KAO[u - 1];
    v = graph->FO[indexU + degU - 1];
    uC = u;
    vC = v;

    p = graph->FORel[2 * graph->edgNumb - 1];
}
*/



//ÑËÓ×ÀÉÍÛÉ ÂÛÁÎĞ ĞÅÁĞÀ


/*
#include <cstdlib>
#include <ctime>

void GraphFactory::choseVerts(Graph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {

    std::srand(std::time(0));
    
    u = rand() % (graph->getVertNumb()-1) + 1;

    int deg = graph->KAO[u] - graph->KAO[u-1];

    int k = (rand() % deg);

    v = graph->FO[graph->KAO[u-1] + k];
    p = graph->FORel[u];

    uC = u;
    vC = v;

}
*/


/*
    ÂÛÁÎĞ ÏÅĞÂÎÃÎ ÏÎÏÀÂÍÅÃÎÑß ĞÅÁĞÀ Ñ ÊÎÍÖÀ Ñ ÌÈÍÈÌÀËÜÍÎÉ ÑÒÅÏÅÍÜŞ   
*/
    
/*
void GraphFactory::choseVerts(Graph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {

    int degU = graph->KAO[graph->vertNumb] - graph->KAO[graph->vertNumb-1];
    int min_deg = 100000;
    int degV;

    u = graph->vertNumb;
    uC = u;

    for (size_t i = 2 * graph->edgNumb - 1; i > 2 * graph->edgNumb - 1 - degU; i--) {
        degV = graph->KAO[graph->FO[i]] - graph->KAO[graph->FO[i-1]];
        
        if(degU + degV < min_deg && graph->isEdge(u, i)) 
        {
            v = graph->FO[i];
            p = graph->FORel[u];
            vC = v;
        }
    }

}
*/



/*
    ÂÛÁÎĞ ĞÅÁĞÀ ÑÎ ÑÒÅÏÅÍÜŞ 6 È ÏĞÎÂÅĞÊÀ ÍÀ 'ÍÀ ÒĞÅÓÃÎËÜÍÎÑÒÜ',
    ÈÍÀ×Å ÁÅĞÅÌ ÏÅĞÂÓŞ ÏÎÏÀÂØÓŞÑÜ Ñ ÊÎÍÖÀ ÂÅĞØÈÍÓ ÑÒÅÏÅÍÜŞ 6,
    ÈÍÀ×Å ÁÅĞÅÌ ÂÅĞØÈÍÓ Ñ ÌÈÍÈÌÀËÜÍÎÉ ÑÒÅÏÅÍÜŞ
*/



void GraphFactory::choseVerts(Graph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {

    int min_degSum = 10000000;

    int tmp_u, tmp_v;

    int flag = 0;
    int cft = 0;

    //âíåøíèé öèêë âûáèğàåò ïåğâóş âåğøèíó
    for(size_t i = graph->vertNumb; i > 0; i--)
    {
        int indexU = graph->KAO[i-1]-1;
        int degU = graph->KAO[i] - graph->KAO[i - 1];
        
       // std::cout << "i = " << i << std::endl;

        //ıòîò öèêë âûáèğàåò âòîğóş âåğøèíó
        for (size_t j = degU; j > 0; j--) {
            int V = graph->FO[indexU + j];
            
            //std::cout <<"V = " << V << std::endl;

            int degV = graph->KAO[V] - graph->KAO[V - 1];

            if (degU + degV == 6) {
                
                //std::cout << cft << std::endl;
                //std::cout << "sum deg == 6\n";
                count_of_edges_degs_6++;
                //èùåò ñìåæíóş âåğèøíó ıòèì äâóì
                for (size_t k = degU; k > 0; k--) {

                    int UV = graph->FO[indexU + k];
                    int indexUV = graph->KAO[UV - 1];

                    if (UV != V)
                    {
                        int degUV = graph->KAO[UV] - graph->KAO[UV - 1];

                        //ïğîáåãàåòñÿ ïî ñìåæíûì âåğøèíàì 3é âåğøèíû
                        for (size_t l = 0; l < degUV; l++) {
                            if (graph->FO[indexUV + l] == V) {

                                count_of_triangles++;

                                u = i;
                                v = V;
                                uC = u;
                                vC = v;
                                p = graph->FORel[int(graph->KAO[i] / 2)-1];
                                if (cft == 0)
                                    count_first_triangles++;

                                return;
                            }
                        }

                    }


                }//for k

                cft++;

                if (flag == 0) {
                    tmp_u = i;
                    tmp_v = V;
                    flag++;
                }

            }//if deg1 + deg2 == 6

            else if(min_degSum > degU + degV){
                min_degSum = degU + degV;
                u = i;
                v = V;
                p = graph->FORel[int(graph->KAO[i] / 2)-1];
                uC = u;
                vC = v;
            }

        }//for j


    }//for i 

    if (flag != 0) {
        u = tmp_u;
        v = tmp_v;
        uC = u;
        vC = v;
        p = graph->FORel[int(graph->KAO[u] / 2)-1];
    }

}


//____________________________________________________FACTORING____________________________________________________



int GraphFactory::recursive_deep = 0;


void GraphFactory::branchingWithNoReverse(Graph *&graph, double &rel, double prob) {

    recursive_deep++;

    u_int u, v, uC, vC;
    double pstFactor, p;

    graph->parallelSeriesTransformation(pstFactor);

    if (graph->getVertNumb() > 5) {
        choseVerts(graph, p, u, v, uC, vC);
        Graph *merge = graph->mergeEdge(u, v);
        Graph *cut = graph->cutEdge(u, v);
        if (cut->isConnective()) {
            branchingWithNoReverse(cut, rel, pstFactor * (1 - p) * prob);
        }
        branchingWithNoReverse(merge, rel, pstFactor * p * prob);
        delete merge;
        delete cut;
    }
    else {     
        rel += pstFactor * prob * graph->baseProbabilities();
        return;
    }
}


double GraphFactory::branching(Graph *&graph) {

    //ïîäñ÷åò ãëóáèíû ğåêóğñèè
    recursive_deep++;

    u_int u, v, uC, vC; 
    double pstFactor, p, result;

    graph->parallelSeriesTransformation(pstFactor);

    if (graph->getVertNumb() > 5) {
        choseVerts(graph, p, u, v, uC, vC);
        Graph *merge = graph->mergeEdge(u, v);
        Graph *cut = graph->cutEdge(u, v);
        if (cut->isConnective()) {
            result = pstFactor * (p * branching(merge) + (1 - p) * branching(cut));
        }
        else {
            //delete cut;    Ñ‡Ñ‚Ğ¾Ğ±Ñ‹ Ğ½ĞµĞ½ÑƒĞ¶Ğ½Ñ‹Ğµ Ğ³Ñ€Ğ°Ñ„Ñ‹ Ğ½Ğµ Ğ·Ğ°Ğ½Ğ¸Ğ¼Ğ°Ğ»Ğ¸ Ğ¼ĞµÑÑ‚Ğ¾ Ğ² Ğ¿Ğ°Ğ¼ÑÑ‚Ğ¸, Ğ¼Ğ¾Ğ¶Ğ½Ğ¾ Ğ¿ĞµÑ€ĞµĞ½ĞµÑÑ‚Ğ¸ ÑƒĞ´Ğ°Ğ»ĞµĞ½Ğ¸Ğµ Ğ¾Ñ‚Ñ€Ğ°Ğ±Ğ¾Ñ‚Ğ°Ğ²ÑˆĞ¸Ñ… Ğ³Ñ€Ğ°Ñ„Ğ¾Ğ² Ğ² Ğ±Ğ»Ğ¾Ğº if-else
            result = pstFactor * (p * branching(merge));
        }
        delete merge;
        delete cut;
        return result;
    }
    else {
        return pstFactor * graph->baseProbabilities();
    }
}

double GraphFactory::factoring(Graph *&graph, FactoringType factoringType) {
    double shootsProbability = 1.0;
    double result = 0.0;
    time_t start, end;
    
    time(&start);
    graph->shootsDeletion(shootsProbability);

    if (graph->isConnective()) {
        if (factoringType == FACTORING) {
            result = shootsProbability * branching(graph);
        }
        else if (factoringType == NO_REVERSE_REC_FACTORING) {
            branchingWithNoReverse(graph, result, 1.0);
            result *= shootsProbability;
        }
    }

    time(&end);

    std::cout << "Work time: " << (int)end - (int)start << std::endl;

    return result;
}









/*
            K GRAPH_FACTORING
*/



//=========================================== CHOOSE VERTEX ============================================





u_int LastNotEmptyVertice(KGraph* G);
int TargetVertexQuality(KGraph* G);

void kGraphFactory::choseLastVerts(KGraph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {
    
    v = LastNotEmptyVertice(graph);

    u = graph->FO[graph->edgNumb * 2 - 1];
    p = graph->FORel[graph->edgNumb * 2 - 1];

}



void kGraphFactory::choseVerts2(KGraph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {
    vC = graph->vertNumb; //vertex count

    u_int deg = 2 * graph->vertNumb + 2;

    for (u_int i = (2 * graph->edgNumb - 1); i >= graph->KAO[1]; i--) {
        uC = graph->FO[i];
        if (i < graph->KAO[vC - 1]) {
            vC--;
        }
        if (vC > uC) {
            if (graph->KAO[uC] - graph->KAO[uC - 1] + graph->KAO[vC] - graph->KAO[vC - 1] < 7) {
                u = uC;
                v = vC;
                p = graph->FORel[i];
                return;
            }
            else if (graph->KAO[uC] - graph->KAO[uC - 1] + graph->KAO[vC] - graph->KAO[vC - 1] < deg) {
                u = uC;
                v = vC;
                p = graph->FORel[i];
                deg = graph->KAO[uC] - graph->KAO[uC - 1] + graph->KAO[vC] - graph->KAO[vC - 1];
            }
        }
    }
}


void kGraphFactory::choseVerts3(KGraph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {

    u = graph->vertNumb;
    u_int indexU = graph->KAO[u - 1];
    u_int degU = graph->KAO[u] - graph->KAO[u - 1];
    v = graph->FO[indexU + degU - 1];
    uC = u;
    vC = v;
    p = graph->FORel[2 * graph->edgNumb - 1];


    if (u > v) {
        uC = u;
        vC = v;
        u = vC;
        v = uC;
    }


}



void kGraphFactory::choseVerts4(KGraph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {

    int min_degSum = 10000000;

    int tmp_u, tmp_v;

    int flag = 0;
    int cft = 0;

    //âíåøíèé öèêë âûáèğàåò ïåğâóş âåğøèíó
    for (size_t i = graph->vertNumb; i > 0; i--)
    {
        int indexU = graph->KAO[i - 1] - 1;
        int degU = graph->KAO[i] - graph->KAO[i - 1];

        // std::cout << "i = " << i << std::endl;

         //ıòîò öèêë âûáèğàåò âòîğóş âåğøèíó
        for (size_t j = degU; j > 0; j--) {
            int V = graph->FO[indexU + j];

            //std::cout <<"V = " << V << std::endl;

            int degV = graph->KAO[V] - graph->KAO[V - 1];

            if (degU + degV == 6) {

                //std::cout << cft << std::endl;
                //std::cout << "sum deg == 6\n";
                //count_of_edges_degs_6++;
                //èùåò ñìåæíóş âåğèøíó ıòèì äâóì
                for (size_t k = degU; k > 0; k--) {

                    int UV = graph->FO[indexU + k];
                    int indexUV = graph->KAO[UV - 1];

                    if (UV != V)
                    {
                        int degUV = graph->KAO[UV] - graph->KAO[UV - 1];

                        //ïğîáåãàåòñÿ ïî ñìåæíûì âåğøèíàì 3é âåğøèíû
                        for (size_t l = 0; l < degUV; l++) {
                            if (graph->FO[indexUV + l] == V) {

                                //count_of_triangles++;

                                u = i;
                                v = V;
                                uC = u;
                                vC = v;
                                if (cft == 0);
                                    //count_first_triangles++;

                                 if (u > v) {
                                     uC = u;
                                     vC = v;
                                     u = vC;
                                     v = uC;
                                 }

                                 for (size_t d = 0; d < graph->KAO[u] - graph->KAO[u - 1]; d++) {
                                     if (graph->FO[graph->KAO[v-1] + d] == u) p = graph->FORel[graph->KAO[v-1] + d];
                                 }


                                return;
                            }
                        }

                    }


                }//for k

                cft++;

                if (flag == 0) {
                    tmp_u = i;
                    tmp_v = V;
                    flag++;
                }

            }//if deg1 + deg2 == 6

            else if (min_degSum > degU + degV) {
                min_degSum = degU + degV;
                u = i;
                v = V;
                uC = u;
                vC = v;

                for (u_int k = 0; k < graph->KAO[u] - graph->KAO[u - 1]; k++) {
                    if (graph->FO[graph->KAO[v-1] + k] == u) p = graph->FORel[graph->KAO[v-1] + k];
                }

            }

        }//for j


    }//for i 

    if (flag != 0) {
        u = tmp_u;
        v = tmp_v;
        uC = u;
        vC = v;

        for (u_int k = 0; k < graph->KAO[u] - graph->KAO[u - 1]; k++) {
            if (graph->FO[graph->KAO[v-1] + k] == u) p = graph->FORel[graph->KAO[v-1] + k];
        }
    }


    if (u > v) {
        uC = u;
        vC = v;
        u = vC;
        v = uC;
    }

}

void kGraphFactory::choseVerts5(KGraph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {

    u = graph->vertNumb;
    u_int indexU = graph->KAO[u - 1];
    u_int degU = graph->KAO[u] - graph->KAO[u - 1];
    v = graph->FO[indexU + degU - 1];
    p = graph->FORel[2 * graph->edgNumb - 1];


    degU = graph->KAO[graph->vertNumb] - graph->KAO[graph->vertNumb - 1];

    u_int min_deg = 100000;
    u_int degV;


    for (size_t i = 2 * graph->edgNumb - 1; i > 2 * graph->edgNumb - 1 - degU; i--) {
        degV = graph->KAO[graph->FO[i]] - graph->KAO[graph->FO[i - 1]];

        if ((degU + degV < min_deg) && graph->isEdge(u, i))
        {
            v = i;

            for (u_int k = 0; k < degU; k++) {
                if(graph->FO[graph->KAO[v-1] + k] == u) p = graph->FORel[graph->KAO[v - 1] + k];
            }

        }
    }

    if (u > v) {
        uC = u;
        vC = v;
        u = vC;
        v = uC;
    }


}


void kGraphFactory::choseRandomVerts(KGraph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC) {

    std::srand(std::time(0));

    u = rand() % (graph->getVertNumb() - 1) + 1;

    int deg = graph->KAO[u] - graph->KAO[u - 1];

    int k = (rand() % deg);

    v = graph->FO[graph->KAO[u - 1] + k];
    
    for (u_int k = 0; k < deg; k++) {
        if (graph->FO[graph->KAO[v - 1] + k] == u) p = graph->FORel[graph->KAO[v-1] + k];
    }

    uC = u;
    vC = v;

    if (u > v) {
        uC = u;
        vC = v;
        u = vC;
        v = uC;
    }


}



//  ============================================FACTORIZATION===============================================================







double kGraphFactory::branching(KGraph*& graph, int variant) {

    if (graph == nullptr) return 1;

    int e = 1;
    double p1 = 1;
    double Rel;
    u_int u, v, uC, vC;
    double p;

    if (graph->vertNumb > 3) {
        if (variant == 0) {
            if (graph->Kconnective() == false) return 0;
        }

        if (TargetVertexQuality(graph) == 1) return 1;

        else {
            graph->KParallelSeriesTransformation(p1);

            if (graph->vertNumb == 2) {
                if (graph->targets[1] == 1 && graph->targets[2] == 1 && graph->edgNumb > 0) return p1 * graph->FORel[1];
                else return p1;
            }
            else if (graph->vertNumb < 3) return p1;
        }

        //choseRandomVerts(graph, p, u, v, uC, vC);
        choseVerts2(graph, p, u, v, uC, vC);

        KGraph* merge = graph->MergeVertex(u, v);
        KGraph* cut = graph->deleteEdge(u, v);

        Rel = p1 * p * branching(merge, 1) + p1 * (1 - p) * branching(cut, 0);

        delete merge;
        delete cut;

    }
    else {
        Rel = graph->baseProbabilities();
    }

    return Rel;

}


void kGraphFactory::branchingWithNoReverse(KGraph*& graph, double& rel, double prob) {


    u_int u, v, uC, vC;
    double pstFactor, p;

    graph->KParallelSeriesTransformation(pstFactor);

    if (graph->getVertNumb() > 3) {
        choseLastVerts(graph, p, u, v, uC, vC);
        KGraph* merge = graph->MergeVertex(u, v);
        KGraph* cut = graph->deleteEdge(u, v);
        if (cut->Kconnective()) {
            branchingWithNoReverse(cut, rel, pstFactor * (1 - p) * prob);
        }
        branchingWithNoReverse(merge, rel, pstFactor * p * prob);
        delete merge;
        delete cut;
    }
    else {
        rel += pstFactor * prob * graph->baseProbabilities();
        //std::cout << "rel = " << rel << std::endl;
        return;

    }
}




//count of target vertices (ENG)
//êîë-âî öåëåâûõ âåğøèí (RU)
int TargetVertexQuality(KGraph* G)
{
    int result = 0;
    for (int i = 1; i < G->getVertNumb() + 1; i++) {
        if (G->getTarget()[i] == 1) {
            result++;
        }
    }
    return result;
}


//last not hanging vertex (ENG)
//ïîñëåäíÿÿ íå âèñÿ÷àÿ âåğøèíà (RU)
u_int LastNotEmptyVertice(KGraph* G)
{
    int i = G->getVertNumb();
    while (G->getKAO()[i] - G->getKAO()[i - 1] == 0) {
        i--;
    }
    return i;
}





}//namespace math