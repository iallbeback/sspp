#include <iostream>
#include "pthread.h"
#include <queue>
#include <cassert>

class MyConcurrentQueue {
public:
    void reserve(int n) {
        queue_limit = n;
    }
    
    void put(int value) {
        pthread_mutex_lock(&m);
        while ((int)queue.size() == queue_limit) {
            pthread_cond_wait(&cp, &m);
        }
        queue.push(value);
        pthread_cond_signal(&cg);
        pthread_mutex_unlock(&m);
    }
    
    int get() {
        pthread_mutex_lock(&m);
        while (queue.empty()) {
            pthread_cond_wait(&cg, &m);
        }
        int n = queue.front();
        queue.pop();
        pthread_cond_signal(&cp);
        pthread_mutex_unlock(&m);
        return n;
    }
    
private:
    int queue_limit;
    std::queue<int> queue;
    pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t cp = PTHREAD_COND_INITIALIZER;
    pthread_cond_t cg = PTHREAD_COND_INITIALIZER;
};

MyConcurrentQueue my_queue;

struct argfun {
    int num;
};

int j = 0;
pthread_mutex_t j_mutex = PTHREAD_MUTEX_INITIALIZER;

void* put_fun(void* arg) {
    argfun* par = (struct argfun*)(arg);
    int n = par->num;
    
    for (int i = 0; i < n; i++) {
        pthread_mutex_lock(&j_mutex);
        int local_j = ++j;
        pthread_mutex_unlock(&j_mutex);

        my_queue.put(local_j);
    }
    
    return NULL;
}

void* get_fun(void* arg) {
    argfun* par = (struct argfun*)(arg);
    int n = par->num;
    
    for (int i = 0; i < n; i++) {
        int a = my_queue.get();

        pthread_mutex_lock(&j_mutex);
        assert(j >= a);
        pthread_mutex_unlock(&j_mutex);
    }
    
    return NULL;
}

int main() {
    my_queue.reserve(100);

    // Тест {1 - 1}
    int num_cycle = 1000000;
    argfun t1;
    t1.num = num_cycle;

    pthread_t thr_put;
    pthread_t thr_get;

    pthread_create(&thr_put, NULL, &put_fun, &t1);
    pthread_create(&thr_get, NULL, &get_fun, &t1);

    pthread_join(thr_put, NULL);
    pthread_join(thr_get, NULL);

    std::cout << "Тест 1 - 1 завершен\n";

    // Тест {1 - N}
    num_cycle = 100000;
    int num_get = 100;
    argfun t2;
    t1.num = num_cycle;
    t2.num = num_cycle / num_get;

    pthread_t thr_put1;
    pthread_t thr_get1[num_get];

    pthread_create(&thr_put1, NULL, &put_fun, &t1);
    for (int i = 0; i < num_get; i++) {
        pthread_create(&thr_get1[i], NULL, &get_fun, &t2);
    }

    pthread_join(thr_put1, NULL);
    for (int i = 0; i < num_get; i++) {
        pthread_join(thr_get1[i], NULL);
    }

    std::cout << "Тест 1 - " << num_get << " завершен\n";

    // Тест {N - 1}
    int num_put = 100;
    t2.num = num_cycle / num_put;

    pthread_t thr_put2[num_put];
    pthread_t thr_get2;

    for (int i = 0; i < num_put; i++) {
        pthread_create(&thr_put2[i], NULL, &put_fun, &t2);
    }
    pthread_create(&thr_get2, NULL, &get_fun, &t1);

    for (int i = 0; i < num_put; i++) {
        pthread_join(thr_put2[i], NULL);
    }
    pthread_join(thr_get2, NULL);

    std::cout << "Тест " << num_put << " - 1 завершен\n";

    // Тест {M - N}
    num_cycle = 500000;
    num_put = 200;
    num_get = 50;
    t1.num = num_cycle / num_put;
    t2.num = num_cycle / num_get;

    pthread_t thr_put3[num_put];
    pthread_t thr_get3[num_get];

    for (int i = 0; i < num_put; i++) {
        pthread_create(&thr_put3[i], NULL, &put_fun, &t1);
    }
    for (int i = 0; i < num_get; i++) {
        pthread_create(&thr_get3[i], NULL, &get_fun, &t2);
    }

    for (int i = 0; i < num_put; i++) {
        pthread_join(thr_put3[i], NULL);
    }
    for (int i = 0; i < num_get; i++) {
        pthread_join(thr_get3[i], NULL);
    }

    std::cout << "Тест " << num_put << " - " << num_get << " завершен\n";

    return 0;
}