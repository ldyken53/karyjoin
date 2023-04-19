#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <time.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) \
    ( BLOCK_LOW((id)+1,p,n)-1 ) 
#define BLOCK_SIZE(id,p,n) \
    (BLOCK_LOW( (id)+1, p, n) - BLOCK_LOW( (id), p, n  ) )
#define BLOCK_OWNER(index,p,n) \
    ( ( ((p)*(index)+1)-1 ) / (n) )

using namespace std;
int main (int argc, char *argv[]) {
    int N = 100000;
    uint *A = new uint[N];
    uint *B = new uint[N];
    srand (time(NULL));
    for(int i = 0; i < N; ++i) { 
        A[i] = rand() % N; 
        B[i] = rand() % N; 
    }
    int joinCount = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (A[i] == B[j]) {
                ++joinCount;
            }
        }
    }
    cout<<"Total join count: "<<joinCount<<'\n';
    MPI_File Afile, Bfile;
    int id, p;
    MPI_Init (&argc, &argv);
    MPI_File_open(MPI_COMM_WORLD , "A", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &Afile);
    MPI_File_open(MPI_COMM_WORLD , "B", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &Bfile);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    MPI_File_write(Afile, A, N, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    MPI_File_write(Bfile, B, N, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    MPI_File_close(&Afile);
    MPI_File_close(&Bfile);
    uint *C = new uint[N];
    uint *D = new uint[N];
    MPI_File_open(MPI_COMM_WORLD , "A", MPI_MODE_RDONLY, MPI_INFO_NULL, &Afile);
    MPI_File_open(MPI_COMM_WORLD , "B", MPI_MODE_RDONLY, MPI_INFO_NULL, &Bfile);
    MPI_File_read(Afile, C, N, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    MPI_File_read(Bfile, D, N, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    for(int i=0; i<4; i++)  // show values read from file
    cout << C[i] << " ";
    cout<<'\n';
    for(int i=0; i<4; i++)  // show values read from file
    cout << D[i] << " ";
    cout<<'\n';
    MPI_File_close(&Afile);
    MPI_Finalize ();
    return 0;
}