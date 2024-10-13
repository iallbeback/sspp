#include <iostream>
#include <stdlib.h>
#include <stdio.h> 
#include <vector>
#include "papi.h"
#include <fstream>

// Глобальные переменные для хранения графа в формате CSR (Compressed Sparse Row)
unsigned int row_count; // количество вершин в графе
unsigned int col_count; // количество рёбер в графе
std::vector<unsigned int> row_ptr; // указатели на начало рёбер для каждой вершины
std::vector<int> col_ids; // идентификаторы вершин, к которым идут рёбра
std::vector<double> vals; // веса рёбер

// Функция для чтения графа из файла
void read_graph(const char* filename) {
    FILE *graph_file = fopen(filename, "rb");
    fread(reinterpret_cast<char*>(&row_count), sizeof(int), 1, graph_file); // читаем количество вершин
    fread(reinterpret_cast<char*>(&col_count), sizeof(unsigned int), 1, graph_file); // читаем количество рёбер

    std::cout << "Row_count = " << row_count << ", col_count = " << col_count << std::endl;

    row_ptr.resize(row_count + 1); // выделяем память для row_ptr
    col_ids.resize(col_count); // выделяем память для col_ids
    vals.resize(col_count); // выделяем память для весов рёбер

    // Читаем данные графа
    fread(reinterpret_cast<char*>(row_ptr.data()), sizeof(unsigned int), row_count + 1, graph_file);
    fread(reinterpret_cast<char*>(col_ids.data()), sizeof(int), col_count, graph_file);
    fread(reinterpret_cast<char*>(vals.data()), sizeof(double), col_count, graph_file);
    
    fclose(graph_file); // закрываем файл
}

// Функция для вывода рёбер вершины с индексом idx
void print_vertex(int idx) {
    for (unsigned int col = row_ptr[idx]; col < row_ptr[idx + 1]; col++) {
        std::cout << col_ids[col] << " " << vals[col] << std::endl; // выводим id соседней вершины и вес ребра
    }
    std::cout << std::endl;
}

// Функция для нахождения вершины с наибольшей суммой весов рёбер к чётным вершинам
unsigned int max_ves_vertex() {
    double max, res = 0;
    unsigned int n = 0;
    
    // Обрабатываем первую вершину
    for (unsigned int j = row_ptr[0]; j < row_ptr[1]; j++) {
        if (col_ids[j] % 2 == 0) res += vals[j]; // если вершина чётная, добавляем вес ребра
    }
    max = res; // задаём начальное максимальное значение

    // Проходим по оставшимся вершинам
    for (unsigned int i = 1; i < row_count; i++) {
        res = 0;
        for (unsigned int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            if (col_ids[j] % 2 == 0) res += vals[j]; // аналогично, сумма весов рёбер к чётным вершинам
        }
        if (res > max) { // если сумма весов больше текущего максимума, обновляем максимальное значение
            max = res;
            n = i; // сохраняем индекс вершины
        }
    }
    return n; // возвращаем вершину с максимальной суммой весов рёбер
}

// Функция для вычисления дополнительного веса для вершины
double W(int n) {
    double res = 0;
    for (unsigned int j = row_ptr[n]; j < row_ptr[n + 1]; j++) {
        res += vals[j] * (row_ptr[col_ids[j] + 1] - row_ptr[col_ids[j]]); // пример вычисления на основе соседних рёбер
    }
    return res;
}

// Функция для нахождения вершины с наибольшим "рангом"
unsigned int max_rang_vertex() {
    double max, res = 0;
    unsigned int n = 0;

    // Обрабатываем первую вершину
    for (unsigned int j = row_ptr[0]; j < row_ptr[1]; j++) {
        res += vals[j] * W(col_ids[j]); // вычисляем ранг
    }
    max = res; // задаём начальный максимум

    // Проходим по оставшимся вершинам
    for (unsigned int i = 1; i < row_count; i++) {
        res = 0;
        for (unsigned int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            res += vals[j] * W(col_ids[j]); // аналогично, вычисляем ранг для каждой вершины
        }
        if (res > max) { // если ранг больше текущего максимума, обновляем максимальное значение
            n = i;
            max = res;
        }
    }
    return n; // возвращаем вершину с наибольшим рангом
}

// Функция для очистки данных графа
void reset_graph() {
    row_count = 0;
    col_count = 0;
    row_ptr.clear();
    col_ids.clear();
    vals.clear();
}

#define N_TESTS 5 // количество тестов

int main() {
	std::ofstream out;
	out.open("test.txt");
	
    int retval, EventSet = PAPI_NULL; // переменные для работы с PAPI
    long long cm[3]; // массив для хранения результатов замеров
    
    // Инициализация библиотеки PAPI
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        std::cout << "PAPI library init error!" << std::endl;
        exit(1);
    }

    // Определение событий для замеров
    int EventCode0 = PAPI_L1_TCM; // промахи в L1 кэше
    int EventCode1 = PAPI_L2_TCM; // промахи в L2 кэше
    int EventCode2; // общее количество выполненных циклов
	//char event_name[] = "perf::CPU-CYCLES";
	char event_name[] = "CPU-CYCLES";
	PAPI_event_name_to_code(event_name, &EventCode2);
    //PAPI_event_name_to_code("perf::PAPI_LD_INS", &EventCode2);

    // Создание набора событий
    PAPI_create_eventset(&EventSet);
    PAPI_add_event(EventSet, EventCode0); // добавляем событие промахов L1
    PAPI_add_event(EventSet, EventCode1); // добавляем событие промахов L2
    PAPI_add_event(EventSet, EventCode2); // добавляем событие выполнения циклов

    // Массив с именами файлов графов
    const char* filenames[N_TESTS] = {"synt", "road_graph", "stanford", "youtube", "syn_rmat"};

    // Цикл по всем тестам
    for (int n_test = 0; n_test < N_TESTS; n_test++) {
        std::cout << filenames[n_test] << std::endl;
		out << filenames[n_test] << std::endl;
        read_graph(filenames[n_test]); // читаем граф из файла
		out << "Row_count = " << row_count << ", col_count = " << col_count << std::endl;

        std::cout << '\n';
		out << '\n';

        // Начинаем замер производительности для функции max_ves_vertex
        PAPI_start(EventSet);
		
        std::cout << "Вершина с наибольшим весом: " << max_ves_vertex() << std::endl;
		out << "Вершина с наибольшим весом: " << max_ves_vertex() << std::endl;
		
        PAPI_stop(EventSet, cm); // останавливаем замеры и сохраняем результаты
        PAPI_reset(EventSet); // сбрасываем счётчики
		
        std::cout << "Количество промахов L1 кэша: \t" << cm[0] << std::endl;
		out << "Количество промахов L1 кэша: \t" << cm[0] << std::endl;
		
        std::cout << "Количество промахов L2 кэша: \t" << cm[1] << std::endl;
		out << "Количество промахов L2 кэша: \t" << cm[1] << std::endl;
		
        std::cout << "Количество выполненных циклов: \t" << cm[2] << std::endl;
		out << "Количество выполненных циклов: \t" << cm[2] << std::endl;

        std::cout << '\n';
		out << '\n';

        // Начинаем замер производительности для функции max_rang_vertex
        PAPI_start(EventSet);        
        std::cout << "Вершина с наибольшим рангом: " << max_rang_vertex() << std::endl;
		out << "Вершина с наибольшим рангом: " << max_rang_vertex() << std::endl;
		
        PAPI_stop(EventSet, cm); // останавливаем замеры и сохраняем результаты
        PAPI_reset(EventSet); // сбрасываем счётчики
		
        std::cout << "Количество промахов L1 кэша: \t" << cm[0] << std::endl;
		out << "Количество промахов L1 кэша: \t" << cm[0] << std::endl;
        std::cout << "Количество промахов L2 кэша: \t" << cm[1] << std::endl;
		out << "Количество промахов L2 кэша: \t" << cm[1] << std::endl;
        std::cout << "Количество выполненных циклов: \t" << cm[2] << std::endl;
		out << "Количество выполненных циклов: \t" << cm[2] << std::endl;

        std::cout << '\n';
		out << '\n';
    }

    // Завершение работы библиотеки PAPI
    PAPI_shutdown();
	out.close();
}
