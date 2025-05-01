#ifndef SRC_GRAPHUTIL_H
#define SRC_GRAPHUTIL_H

#include "base_def.h"

namespace math {

class Graph;
class KGraph;


typedef enum {
    FACTORING,
    NO_REVERSE_REC_FACTORING,
   //NO_RECURSION_FACTORING
} FactoringType;

typedef enum {
    SERIAL
    //PARALLEL
} FactoringExecutionType;

/**
 * Класс-генератор различных графов
*/
class GraphGenerator {
    public:
        /**
         * Генерирует полный граф на n вершинах
         * @param n количество вершин
         * @return Указатель на сгенерированный граф
        */
        static Graph *generateCompleteGraph(u_int n);

        /**
         * Генерирует граф-решётку n * m
         * @param n количество рёбер в длину
         * @param m количество рёбер в ширину
         * @return Указатель на сгенерированный граф
        */
        static Graph *generateGrid(u_int n, u_int m = 0);

        /**
         * Генерирует граф-решётку n * m, 
         * где квадратные блоки 1 * 1 соединены диагональными рёбрами
         * @param n количество рёбер в длину
         * @param m количество рёбер в ширину
         * @return Указатель на сгенерированный граф
        */
        static Graph *generateGridWithDiags(u_int n, u_int m = 0);

        /**
         * Генерирует случайный граф (в разработке)
         * @param vNumb количество вершин
         * @param eNumb количество рёбер
         * @return Указатель на сгенерированный граф
        */
        static Graph *generateGraph(u_int vNumb, u_int eNumb);

        /**
         * Генерирует граф модели GEANT
         * @return Указатель на сгенерированный граф
        */
        static Graph *generateGEANT();
};



/**
 * Класс для факторизации графа
 * Расширять можно через добавление новых типов факторизации в enum FactoringType 
 * и добавление соответствующего статического private-метода
*/
class GraphFactory {
    

    /**
     * Выбор ребра в графе для стягивания/удаления
     * @param graph граф рассматриваемого шага факторизации
     * @param prob вероятность присутствия выбранного ребра в графе
     * @param u первая вершина ребра
     * @param v вторая вершина ребра
     * @param uC вспомогательный параметр
     * @param vC вспомогательный параметр
    */
    static void choseVerts(Graph *&graph, double &p, u_int &u, u_int &v, u_int &uC, u_int &vC);
    

    /*
     * Факторизация с передачей вероятности вглубь рекурсии
     * @param graph граф рассматриваемого шага факторизации
     * @param rel общая вероятность
     * @param prob вероятность, передаваемая в следующий шаг рекурсии. Произведение вероятности последовательно-параллельного преобразования, 
     * вероятности присутствия/отсутствия выбранного ребра в графе, вероятности на предыдущем шаге)
    */
    static void branchingWithNoReverse(Graph *&graph, double &rel, double prob);

    /**
     * Обычное ветвление с возвращением при получении графа малой размерности либо связного графа
     * @param graph граф рассматриваемого шага факторизации
     * @return произведение вероятностей на пути к полученному графу, умноженное на его вероятность
    */
    static double branching(Graph *&graph);
    
public:

    /**
     * Стартовая функция для факторизации с возможностью указания типа факторизации
     * @param graph исходный граф
    */
    static double factoring(Graph *&graph, FactoringType factoringType = NO_REVERSE_REC_FACTORING);
    static int count_of_triangles;
    static int count_of_edges_degs_6;
    static int recursive_deep;
    static int count_first_triangles;


};//Graph FACTORY




/*
                                        MY PART OF CODE 
*/

class kGraphFactory {
public:

    static void choseVerts(KGraph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC); 
    static void choseVerts2(KGraph*& graph, double& p, u_int& u, u_int& v, u_int& uC, u_int& vC);
    static double branching(KGraph*& graph);
    static double branching(KGraph*& graph, int variant);

};//kGraphFactory


}//namespace math

#endif