#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>
#include <time.h>
#include <sys/times.h>

//#define cntThreads 4

struct ArgsThread
{
    long long left,right;
    double step;
    double partialSum;
};

static void *worker(void *ptrArgs)
{
    ArgsThread * args = reinterpret_cast<ArgsThread *>(ptrArgs);
    double x;
    double sum=0.;
    double step=args->step;
    for (long long i=args->left; i<args->right; i++)
    {
        x = (i + .5)*step;
        sum = sum + 4.0/(1.+ x*x);
    }
    args->partialSum=sum*step;
    return NULL;
}

int main(int argc, char** argv)
{
    long long num_steps;
    sscanf(argv[1], "%lld", &num_steps);
    int cntThreads;
    sscanf(argv[2], "%d", &cntThreads);
    
    pthread_t threads[4];
    ArgsThread arrArgsThread[4];
    
    clock_t clockStart, clockStop;
    tms tmsStart, tmsStop;
    
    clockStart = times(&tmsStart);
    double step = 1./(double)num_steps;
    long long cntStepsPerThread= num_steps / cntThreads;
    for (int idThread=0; idThread<cntThreads; idThread++)
    {
        arrArgsThread[idThread].left  = idThread*cntStepsPerThread;
        arrArgsThread[idThread].right = (idThread+1)*cntStepsPerThread;
        arrArgsThread[idThread].step = step;
        if (pthread_create(&threads[idThread], NULL, worker, &arrArgsThread[idThread]) != 0)
        {
            return EXIT_FAILURE;
        }
    }
    double pi=0.;
    for (int idThread=0; idThread<cntThreads; idThread++)
    {
        if (pthread_join(threads[idThread], NULL) != 0)
        {
            return EXIT_FAILURE;
        }
        pi +=arrArgsThread[idThread].partialSum;
    }
    clockStop = times(&tmsStop);
    
    std::cout << "PI\t " << pi << std::endl;
    std::cout << "Time\t " ;
    double sex = (clockStop - clockStart)/static_cast<double>(sysconf(_SC_CLK_TCK));
    std::cout << sex << "s\n" << std::endl;
    
    std::cout << "kers\t" << "times\t" << "boost\t" <<std::endl;
    double first;
    for(int ker = 1; ker<=4; ker++){
    	double avg = 0;
    	for(int n = 0; n<3; n++){
    		    clockStart = times(&tmsStart);
		    double step = 1./(double)num_steps;
		    long long cntStepsPerThread= num_steps / ker;
		    for (int idThread=0; idThread<ker; idThread++)
		    {
			arrArgsThread[idThread].left  = idThread*cntStepsPerThread;
			arrArgsThread[idThread].right = (idThread+1)*cntStepsPerThread;
			arrArgsThread[idThread].step = step;
			if (pthread_create(&threads[idThread], NULL, worker, &arrArgsThread[idThread]) != 0)
			{
			    return EXIT_FAILURE;
			}
		    }
		    double pi=0.;
		    for (int idThread=0; idThread<ker; idThread++)
		    {
			if (pthread_join(threads[idThread], NULL) != 0)
			{
			    return EXIT_FAILURE;
			}
			pi +=arrArgsThread[idThread].partialSum;
		    }
		    clockStop = times(&tmsStop);
		    avg += (clockStop - clockStart)/static_cast<double>(sysconf(_SC_CLK_TCK));
    	}
    	avg /= 3.;
    	if(ker == 1) {first = avg;}
    	std::cout << ker << "\t" << std::fixed << std::setprecision(4) << avg << "\t" << first/avg << std::endl;
    }
    return 0;
}
