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

struct Row {
    uint x, y;
};

using namespace std;
int main (int argc, char *argv[]) {
    int N = 1000;
    Row *A = new Row[N];
    Row *B = new Row[N];
    srand (time(NULL));
    for(int i = 0; i < N; ++i) { 
        A[i].x = rand() % N; 
        A[i].y = rand() % N; 
        B[i].x = rand() % N; 
        B[i].y = rand() % N; 
    }
    int joinCount = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (A[i].x == B[j].x) {
                ++joinCount;
            }
        }
    }
    cout<<"Total join count: "<<joinCount<<'\n';
    MPI_File Afile, Bfile;
    int id, p;
    MPI_Init (&argc, &argv);
    MPI_Datatype dt_row;
    MPI_Type_contiguous(2, MPI_UNSIGNED, &dt_row);
    MPI_Type_commit(&dt_row);
    MPI_File_open(MPI_COMM_WORLD , "A", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &Afile);
    MPI_File_open(MPI_COMM_WORLD , "B", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &Bfile);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    MPI_File_write(Afile, A, N, dt_row, MPI_STATUS_IGNORE);
    MPI_File_write(Bfile, B, N, dt_row, MPI_STATUS_IGNORE);
    MPI_File_close(&Afile);
    MPI_File_close(&Bfile);
    Row *C = new Row[N];
    Row *D = new Row[N];
    MPI_File_open(MPI_COMM_WORLD , "A", MPI_MODE_RDONLY, MPI_INFO_NULL, &Afile);
    MPI_File_open(MPI_COMM_WORLD , "B", MPI_MODE_RDONLY, MPI_INFO_NULL, &Bfile);
    MPI_File_read(Afile, C, N, dt_row, MPI_STATUS_IGNORE);
    MPI_File_read(Bfile, D, N, dt_row, MPI_STATUS_IGNORE);
    for(int i=0; i<4; i++)  // show values read from file
    cout << C[i].x << " ";
    cout<<'\n';
    for(int i=0; i<4; i++)  // show values read from file
    cout << D[i].y << " ";
    cout<<'\n';
    MPI_File_close(&Afile);
    MPI_Finalize ();
    return 0;
}