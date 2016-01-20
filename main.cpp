/**
 * Usage:
 * ./msg-spr-sim [<no_of_threads_per_block>] [<device_id>] <”<input_file>” >”<output_file>” 2>”<file_with_exec_time>”
 * ./msg-spr-sim [...] | dot -Tpng -o graph.png
 *
 * Input file format:
 * 		| %d			// #V
 * x #V	| %f %f %f %f 	// v_h_i, G_0_i, G_max_i, v_d_i
 * 		| %d			// #E
 * x #V	| %d [%d, ...]	// #edges of v_i [v_j, ...]
 * 		| %d			// source v_i
 * 		| %d %d %d 		// t_c, t_p, t_s
 *
 * Output file format: dot
 * Stderr output: single double number - execution time
 */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <sys/time.h>

//#define _DEBUG
#define ALLOC alloc_if(1)
#define REUSE alloc_if(0)
#define FREE free_if(1)
#define RETAIN free_if(0)

extern double elapsedTime(void) {
    struct timeval t;
    gettimeofday(&t, 0);
    return ((double)t.tv_sec + ((double)t.tv_usec / 1000000.0));
}

template <typename T>
inline void printArray(const T *array, const int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << array[i] << "\t";
    }
    std::cout << std::endl;
}

template <typename T>
inline void arrayMap(T *array, const int size, void (*callback)(T&) ) {
    for (int i = 0; i < size; ++i) {
        callback(array[i]);
    }
}

struct NodeData {
    float v_h;
    float G_0;
    float G_max;
    float v_d;
    float v_r;
    int   last_t;
    bool  send;
};

inline void printG(NodeData &data) {
    std::cout << std::fixed << std::setprecision(2) << data.G_0 << "\t";
}


__attribute__((target(mic))) void recv(const int N, const int *V, NodeData *Vdata, const int Elen, const int *E,
                     int *M, const int t_c, const int t_p, const int threads) {

    const int iter = (N + threads -1) / threads;
    omp_set_dynamic(0);
    #pragma omp parallel num_threads(threads) shared(iter, N, V, Vdata, E, M, t_c, t_p)
    {
        int lastId = (omp_get_thread_num() + 1) * iter,
                tid = omp_get_thread_num() * iter;
        if (lastId > N) {
            lastId = N;
        }

//        #pragma ivdep
//        #pragma vector aligned
        for (; tid < lastId; ++tid) {
            NodeData *data = &Vdata[tid];
            if (data->send) {
                continue;
            }

            int lastTime = data->last_t;
            int msgCount = 0;
            int start = V[tid];
            int end = V[tid + 1];
            int msg;

//            #pragma ivdep
//            #pragma vector aligned
            // reading messages
            for (int i = start; i < end; ++i) {
                msg = M[i];

                if (msg <= 0 || msg <= data->last_t) {
                    continue;
                }

                if (lastTime < msg) {
                    lastTime = msg;
                }
                ++msgCount;
            }

            data->last_t = lastTime;

            // processing messages
            data->G_0 = data->G_max - (data->G_max - data->G_0) * expf(-0.01f * msgCount * t_p);
        }

        //
        #pragma omp barrier
        //

        tid = omp_get_thread_num() * iter;

//        #pragma ivdep
//        #pragma vector aligned
        for (; tid < lastId; ++tid) {
            NodeData *data = &Vdata[tid];
            if (data->send) {
                continue;
            }

            int start = V[tid];
            int end = V[tid + 1];
            int lastTime, start2, end2, v;

            // sending messages
            if (data->G_0 * (data->v_h + data->v_r) >= data->v_d) {
                data->send = true;
                lastTime = data->last_t + t_p + t_c;

//                #pragma ivdep
//                #pragma vector aligned
                for (int i = start; i < end; ++i) {
                    v = E[i];
                    start2 = V[v];
                    end2 = V[v + 1];
//                    #pragma ivdep
//                    #pragma vector aligned
                    for (int j = start2; j < end2; ++j) {
                        if (E[j] == tid) {
                            M[j] = lastTime;
                            break;
                        }
                    }
                }
            }
        } // end of for
    } // end of pragma
}

// those must be globals due to nocopy
__attribute__((target(mic))) int threads;
__attribute__((target(mic))) int N;
__attribute__((target(mic))) int t_c;
__attribute__((target(mic))) int t_p;

int main(int argc, char* argv[]) {

    // obtain available devices
    int deviceCount = 1;

    if (deviceCount <= 0) {
        std::cerr << "No Xeon Phi devices has been found!" << std::endl;
        exit(1);
    }

    // read parameters
    threads = -1;
    int deviceId = -1;

    if (argc >= 2) {
        threads = atoi(argv[1]);
    }
    if (argc >= 3) {
        deviceId = atoi(argv[2]);
    }

    if (deviceId == -1) {
        deviceId = 0;
    }
    if (threads == -1) {
        threads = 240;
    }

    if (deviceId < 0 || deviceId >= deviceCount) {
        std::cerr << "Device id out of range!" << std::endl;
        exit(3);
    }

    // read input data
    //int N;				// no. of vertices (nodes)
    int Elen,			// no. of edges
        Vsrc;			// source vertex
    //int t_c = 3,		// communication time [s]
    //    t_p = 30,		// processing time [s] e.g. short movie
    int t_s = 330;		// max. simulation time [s]
    float v_h = 0,		// reflection potential [-1; 1]
        G_0 = 1,		// initial conductivity [0, G_max)
        G_max = 100,	// max. conductivity (G_0, +inf)
        v_r = 1,		// registration potential (0, +inf)
        v_d = 1;		// decision-making potential (0, +inf)

    int *V,
        *E,
        *M;
    NodeData *Vdata;

#ifdef _DEBUG
    std::cout << "/*\nProgram started." << std::endl;
#endif

    //
    std::cin >> N;

    // N+1 - space for one extra element at the end (for easiest iteration through graph)
    posix_memalign((void**)V, 64, (N+1) * sizeof(int));

    posix_memalign((void**)Vdata, 64, N * sizeof(NodeData));

    for (int i = 0; i < N; ++i) {
        std::cin >> v_h >> G_0 >> G_max >> v_d;
        Vdata[i].v_h = v_h;
        Vdata[i].G_0 = G_0;
        Vdata[i].G_max = G_max;
        Vdata[i].v_d = v_d;
        Vdata[i].v_r = v_r;
        Vdata[i].last_t = 0;
        Vdata[i].send = false;
    }

    //
    std::cin >> Elen;

    posix_memalign((void**)&E, 64, Elen * sizeof(int));

//    M = (int*) calloc(Elen, sizeof(int));	// zero-initialized
    posix_memalign((void**)M, 64, Elen * sizeof(int));
    for (int i = 0; i < Elen; ++i) {
        M[i] = 0;
    }

    V[0] = 0;
    {
        int Elen_i, start, end, e_i;
        for (int i = 0; i < N; ++i) {
            std::cin >> Elen_i;
            start = V[i];
            end = start + Elen_i;
            V[i + 1] = end;

            for (int j = start; j < end; ++j) {
                std::cin >> e_i;
                E[j] = e_i;
            }
        }
    }

    // get source vector and create M(0)
    std::cin >> Vsrc;
    Vdata[Vsrc].send = true;

    {
        int start = V[Vsrc];
        int end = V[Vsrc + 1];
        int start2, end2, v;

        for (int i = start; i < end; ++i) {
            v = E[i];
            start2 = V[v];
            end2 = V[v + 1];
            for (int j = start2; j < end2; ++j) {
                if (E[j] == Vsrc) {
                    M[j] = 1;
                    break;
                }
            }
        }
    }

    //
    std::cin >> t_c >> t_p >> t_s;

#ifdef _DEBUG
    std::cout << "V:\t";
	printArray(V, N);

	std::cout << "E:\t";
	printArray(E, Elen);

	std::cout << "M_0:\t";
	printArray(M, Elen);
#endif

    // capture the start time
    double startTime = elapsedTime();

    #pragma offload target(mic) \
            in(N, t_c, t_p, threads: ALLOC RETAIN) \
            in(V:length(N+1) ALLOC RETAIN) \
            in(E:length(Elen) ALLOC RETAIN) \
            in(Vdata:length(N) ALLOC RETAIN) \
            in(M:length(Elen) ALLOC RETAIN)
    {
        recv(N, V, Vdata, Elen, E, M, t_c, t_p, threads);
    }

    int t_end = t_s - t_c + t_p;
    for (int t = 1 + t_c + t_p; t <= t_end; t += t_c + t_p) {
        #pragma offload target(mic) \
            nocopy(N, t_c, t_p, threads: REUSE RETAIN) \
            nocopy(V:length(N+1) REUSE RETAIN) \
            nocopy(E:length(Elen) REUSE RETAIN) \
            nocopy(Vdata:length(N) REUSE RETAIN) \
            nocopy(M:length(Elen) REUSE RETAIN)
        {
            recv(N, V, Vdata, Elen, E, M, t_c, t_p, threads);
        }
    }

    #pragma offload target(mic) \
            nocopy(N, t_c, t_p, threads: REUSE FREE) \
            nocopy(V:length(N+1) REUSE FREE) \
            nocopy(E:length(Elen) REUSE FREE) \
            nocopy(Vdata:length(N) REUSE FREE) \
            out(M:length(Elen) REUSE FREE)
    {
        recv(N, V, Vdata, Elen, E, M, t_c, t_p, threads);
    }

    // capture the end time
    double endTime = elapsedTime();

#ifdef _DEBUG
    std::cout << "M: ";
	printArray(M, Elen);
	std::cout << "*/\n";
#endif

    // write computing time (to cerr for simplicity)
    std::cerr << std::fixed << (endTime - startTime) << std::endl;

    // generate output
    // 1. how many recipients
    // 2. max. distance from source (range)
    // -> 3. a graph at the end of the simulation in dot format

    std::cout << "digraph G {\n";
    std::cout << "\tnode [fontsize=12]\n";
    std::cout << "\tedge [fontcolor=\"0.5 0.5 0.5\",fontsize=8]\n";
    std::cout << "\t" << Vsrc << " [label=\"src\"]\n\n";

    int start, end;
    for (int i = 0; i < N; ++i) {
        start = V[i];
        end = V[i + 1];
        for (int j = start; j < end; ++j) {
            std::cout << "\t" << E[j] << " -> " << i;
            if (M[j] > 0) {
                std::cout << " [label=" << M[j] << "]\n";
            } else {
                std::cout << " [style=dotted]\n";
            }
        }
    }


    std::cout << "}" << std::endl;

    // free memory on the CPU side
    free(V);
    free(E);
    free(M);
    free(Vdata);
}
