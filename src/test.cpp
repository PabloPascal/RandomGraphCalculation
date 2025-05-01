#include "../headers/GraphUtil.h"
#include "../headers/Graph.h"
#include "../headers/KGraph.h"

#include <iomanip>
#include <iostream>
#include <chrono>
#include <string>

void print_graph(math::KGraph* ptrGraph);
void print_graph(math::Graph* ptrGraph);



void test1() {

    std::ifstream input("kGraphFiles/triangle.txt");       //лучше указать абсолютный путь

    if (!input.is_open()) {
        std::cout << "DIDNT OPEN! \n";
        return;
    }

    math::Graph graph;
    input >> graph;


    //math::Graph *graphPtr = math::GraphGenerator::generateCompleteGraph(12);
    math::Graph* graphPtr = &graph;

    print_graph(graphPtr);

    //timer
    auto begin = std::chrono::steady_clock::now();

    double rel = math::GraphFactory::factoring(graphPtr, math::FactoringType::NO_REVERSE_REC_FACTORING);

    auto end = std::chrono::steady_clock::now();

    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "The time: " << elapsed_ms.count() << " ms\n";

    std::cout << std::fixed << std::setprecision(25);
    std::cout << rel << std::endl;
    std::cout << "deep recursive = " << math::GraphFactory::recursive_deep << std::endl;
    std::cout << "count of triangles = " << math::GraphFactory::count_of_triangles << std::endl;
    std::cout << "count_of_edges_degs_6 = " << math::GraphFactory::count_of_edges_degs_6 << std::endl;
    std::cout << "first triangles = " << math::GraphFactory::count_first_triangles << std::endl;

    //delete graphPtr;
}


void test2() {

    std::string path = "../kGraphFiles/";

    //std::ifstream input(path + "kgraph_test2.txt");
    std::ifstream input(path + "/k_grid33.txt");
    //std::ifstream input(path + "/kgraph_test.txt");
    //std::ifstream input(path + "/triangle.txt");
    if(!input.is_open()){
        std::cout << "can't open the file!" << std::endl;
        return;
    }

    math::KGraph graph;
    input >> graph;
    input.close();

    math::KGraph* ptrGraph = &graph;

    std::cout << "==============input Graph==========" << std::endl;
    print_graph(ptrGraph);

    auto begin = std::chrono::steady_clock::now();

    //double rel = math::kGraphFactory::branching(ptrGraph);
    double rel = math::kGraphFactory::branching(ptrGraph,0);

    auto end = std::chrono::steady_clock::now();

    std::cout << "\nRel = " << rel << std::endl;

}


void func_test() {

    //std::ifstream input("src/complete_graph4.txt");
    std::ifstream input("src/k_grid33.txt");
    //std::ifstream input("src/test2.txt");
    if (!input.is_open()) {
        std::cout << "can't open the file!" << std::endl;
        return;
    }

    math::KGraph graph;
    input >> graph;

    input.close();

    math::KGraph* ptrGraph = &graph;

    u_int u, v, uC, vC;
    double pstFactor, p, result;
    
    //math::kGraphFactory::choseVerts2(ptrGraph, p, u, v, uC, vC);
    
    //math::KGraph* merge = ptrGraph->MergeVertex(2, 3);
    //math::KGraph* del = ptrGraph->deleteEdgeK(2, 3);
    //print_graph(merge);
    //print_graph(del);
    math::KGraph* PST = ptrGraph->KParallelSeriesTransformation(ptrGraph, p);
    
    print_graph(ptrGraph);
    //ptrGraph->KParallelSeriesTransformation(p);
    print_graph(PST);


}





void print_graph(math::KGraph* ptrGraph) {

    if (ptrGraph == nullptr) return;

    std::cout << ptrGraph->getVertNumb() << std::endl;
    std::cout << ptrGraph->getEdgeNumb() << std::endl;

    for (u_int i = 0; i < ptrGraph->getVertNumb() + 1; i++) {
        std::cout << ptrGraph->getKAO()[i] << " ";
    }
    std::cout << "\n";
    for (u_int i = 0; i < 2 * ptrGraph->getEdgeNumb(); i++) {
        std::cout << ptrGraph->getFO()[i] << " ";
    }
    std::cout << "\n";
    for (u_int i = 1; i < ptrGraph->getVertNumb()  + 1; i++) {
        std::cout << ptrGraph->getTarget()[i] << " ";
    }


    std::cout << "\n";

}



void print_graph(math::Graph* ptrGraph) {

    if (ptrGraph == nullptr) return;

    std::cout << ptrGraph->getVertNumb() << std::endl;
    std::cout << ptrGraph->getEdgNumb() << std::endl;

    for (u_int i = 0; i < ptrGraph->getVertNumb() + 1; i++) {
        std::cout << ptrGraph->getKAO()[i] << " ";
    }
    std::cout << "\n";
    for (u_int i = 0; i < 2 * ptrGraph->getEdgNumb(); i++) {
        std::cout << ptrGraph->getFO()[i] << " ";
    }
    std::cout << "\n";
    for (u_int i = 0; i < 2*ptrGraph->getEdgNumb(); i++) {
        std::cout << ptrGraph->getForel()[i] << " ";
    }


    std::cout << "\n";

}
