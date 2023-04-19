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
    MPI_File Afile, Bfile;
    int id, p;
    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    int local_size = BLOCK_SIZE(id, p, N);
    int offset = BLOCK_LOW(id, p, N);
    // cout<<"Process "<<id<<" "<<local_size<<" "<<offset<<'\n';
    uint *A = new uint[local_size];
    uint *B = new uint[local_size];

    MPI_File_open(MPI_COMM_WORLD , "A", MPI_MODE_RDONLY, MPI_INFO_NULL, &Afile);
    MPI_File_read_at_all(Afile, sizeof(uint) * offset, A, local_size, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    MPI_File_close(&Afile);
    MPI_File_open(MPI_COMM_WORLD , "B", MPI_MODE_RDONLY, MPI_INFO_NULL, &Bfile);
    MPI_File_read_at_all(Bfile, sizeof(uint) * offset, B, local_size, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    MPI_File_close(&Bfile);

    // Need arrays for count of elements to go to each process
    int *Acounts = new int[p];
    int *Bcounts = new int[p];
    // Need arrays for start offsets of elements to go to each process
    int *Aoffsets = new int[p];
    int *Boffsets = new int[p];
    // Need arrays for start offsets of elements to receive from each process
    int *Aoffsetsr = new int[p];
    int *Boffsetsr = new int[p];
    // Need arrays for counters to put elements in right place for send buffer
    int *Acounters = new int[p];
    int *Bcounters = new int[p];

    // Initialize buffers with 0
    for (int i = 0; i < p; ++i) {
        Acounts[i] = 0;
        Bcounts[i] = 0;
        Aoffsets[i] = 0;
        Boffsets[i] = 0;
        Aoffsetsr[i] = 0;
        Boffsetsr[i] = 0;
        Acounters[i] = 0;
        Bcounters[i] = 0;
    }

    // Count the number of elements this process needs to send to others
    for (int i = 0; i < local_size; ++i) {
        ++Acounts[A[i] % p]; 
        ++Bcounts[B[i] % p]; 
    }

    // Prefix sum on send counts to get offsets
    for (int i = 1; i < p; ++i) {
        Aoffsets[i] = Aoffsets[i - 1] + Acounts[i - 1];
        Boffsets[i] = Boffsets[i - 1] + Bcounts[i - 1];
    }

    // Actually divide elements into where they're supposed to go
    uint *sendA = new uint[local_size];
    uint *sendB = new uint[local_size];
    int sendProcA, sendProcB;
    for (int i = 0; i < local_size; ++i) {
        sendProcA = A[i] % p;
        sendA[Aoffsets[sendProcA] + Acounters[sendProcA]] = A[i];
        ++Acounters[sendProcA];
        sendProcB = B[i] % p;
        sendB[Boffsets[sendProcB] + Bcounters[sendProcB]] = B[i];
        ++Bcounters[sendProcB];
    }

    // All to all comm to get how many this process will receive from others
    int *Acountsr = new int[p];
    int *Bcountsr = new int[p];
    MPI_Alltoall(Acounts, 1, MPI_INT, Acountsr, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoall(Bcounts, 1, MPI_INT, Bcountsr, 1, MPI_INT, MPI_COMM_WORLD);

    // Compute how many elements this process is responsible for joining
    int Asize = 0;
    int Bsize = 0;
    for (int i = 0; i < p; ++i) {
        Asize += Acountsr[i];
        Bsize += Bcountsr[i];
    }  

    // Prefix sum on receive counts to get offsets
    for (int i = 1; i < p; ++i) {
        Aoffsetsr[i] = Aoffsetsr[i - 1] + Acountsr[i - 1];
        Boffsetsr[i] = Boffsetsr[i - 1] + Bcountsr[i - 1];
    }

    // Create buffers with correct size and get elements this process will join
    uint *Af = new uint[Asize];
    uint *Bf = new uint[Bsize];
    MPI_Alltoallv(sendA, Acounts, Aoffsets, MPI_UNSIGNED, Af, Acountsr, Aoffsetsr, MPI_UNSIGNED, MPI_COMM_WORLD);
    MPI_Alltoallv(sendB, Bcounts, Boffsets, MPI_UNSIGNED, Bf, Bcountsr, Boffsetsr, MPI_UNSIGNED, MPI_COMM_WORLD);

    // Logging to check correctness
    // cout<<"Process "<<id<<" sizes: ";
    // cout<<"A "<<Asize<<" B "<<Bsize<<"\t";
    // cout<<"\n";
    for (int i = 0; i < Asize; ++i) {
        if (Af[i] % p != id) {
            cout<<"ERROR: process "<<id<<" has element meant for process "<<Af[i] % p<<'\n';
            MPI_Finalize ();
            break;
        }
    }
    for (int i = 0; i < Bsize; ++i) {
        if (Bf[i] % p != id) {
            cout<<"ERROR: process "<<id<<" has element meant for process "<<Bf[i] % p<<'\n';
            break;
        }
    }

    // Count number of joined elements
    int joinCount = 0;
    for (int i = 0; i < Asize; ++i) {
        for (int j = 0; j < Bsize; ++j) {
            if (Af[i] == Bf[j]) {
                ++joinCount;
            }
        }
    }
    cout<<"Process "<<id<<" join count: "<<joinCount<<"\n";
    MPI_Finalize ();
    return 0;
}