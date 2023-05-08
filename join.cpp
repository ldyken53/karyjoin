#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include "btree/btree_set.h"

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
struct Row_Compare {
    bool operator()(const Row a, const Row b) const {
        // Compare based on first column
        if (a.x < b.x) {
            return true;
        } else if (a.x > b.x) {
            return false;
        }
        // If first column equal, use second
        if (a.y < b.y) {
            return true;
        } else if (a.y > b.y) {
            return false;
        }
        // If rows equal, return false
        return false;
    }
};

using namespace std;
int main (int argc, char *argv[]) {
    // Init
    int N = 1000000;
    MPI_File Afile, Bfile;
    int id, p;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    int local_size = BLOCK_SIZE(id, p, N);
    int offset = BLOCK_LOW(id, p, N);
    // cout<<"Process "<<id<<" "<<local_size<<" "<<offset<<'\n';

    // Make MPI datatype for row struct
    MPI_Datatype dt_row;
    MPI_Type_contiguous(2, MPI_UNSIGNED, &dt_row);
    MPI_Type_commit(&dt_row);
    Row *A = new Row[local_size];
    Row *B = new Row[local_size];

    // Open and read input files
    MPI_File_open(MPI_COMM_WORLD , "A", MPI_MODE_RDONLY, MPI_INFO_NULL, &Afile);
    MPI_File_read_at_all(Afile, sizeof(Row) * offset, A, local_size, dt_row, MPI_STATUS_IGNORE);
    MPI_File_close(&Afile);
    MPI_File_open(MPI_COMM_WORLD , "B", MPI_MODE_RDONLY, MPI_INFO_NULL, &Bfile);
    MPI_File_read_at_all(Bfile, sizeof(Row) * offset, B, local_size, dt_row, MPI_STATUS_IGNORE);
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
        ++Acounts[A[i].x % p]; 
        ++Bcounts[B[i].x % p]; 
    }

    // Prefix sum on send counts to get offsets
    for (int i = 1; i < p; ++i) {
        Aoffsets[i] = Aoffsets[i - 1] + Acounts[i - 1];
        Boffsets[i] = Boffsets[i - 1] + Bcounts[i - 1];
    }

    // Actually divide elements into where they're supposed to go
    Row *sendA = new Row[local_size];
    Row *sendB = new Row[local_size];
    int sendProcA, sendProcB;
    for (int i = 0; i < local_size; ++i) {
        sendProcA = A[i].x % p;
        sendA[Aoffsets[sendProcA] + Acounters[sendProcA]] = A[i];
        ++Acounters[sendProcA];
        sendProcB = B[i].x % p;
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
    Row *recvA = new Row[Asize];
    Row *recvB = new Row[Bsize];
    MPI_Alltoallv(sendA, Acounts, Aoffsets, dt_row, recvA, Acountsr, Aoffsetsr, dt_row, MPI_COMM_WORLD);
    MPI_Alltoallv(sendB, Bcounts, Boffsets, dt_row, recvB, Bcountsr, Boffsetsr, dt_row, MPI_COMM_WORLD);

    // Logging to check correctness
    // cout<<"Process "<<id<<" sizes: ";
    // cout<<"A "<<Asize<<" B "<<Bsize<<"\t";
    // cout<<"\n";
    // for (int i = 0; i < Asize; ++i) {
    //     if (recvA[i].x % p != id) {
    //         cout<<"ERROR: process "<<id<<" has element meant for process "<<recvA[i].x % p<<'\n';
    //         MPI_Finalize ();
    //         break;
    //     }
    // }
    // for (int i = 0; i < Bsize; ++i) {
    //     if (recvB[i].x % p != id) {
    //         cout<<"ERROR: process "<<id<<" has element meant for process "<<recvB[i].x % p<<'\n';
    //         break;
    //     }
    // }

    // Insert inner relation into btree_set for hash join 
    // TODO: Should this be multiset to allow duplicate entries in B? or can we assume deduplicated relations
    double start_time_set = MPI_Wtime();
    btree::btree_multiset<Row, Row_Compare> Bset; 
    for (int i = 0; i < Bsize; ++i) {
        Bset.insert(recvB[i]);
    }
    double end_time = MPI_Wtime();
    // cout<<"Time to insert B into set: "<<end_time - start_time_set<<"\n";

    // Count number of joined elements manually
    // start_time = MPI_Wtime();
    // int joinCount = 0;
    // for (int i = 0; i < Asize; ++i) {
    //     for (int j = 0; j < Bsize; ++j) {
    //         if (recvA[i].x == recvB[j].x) {
    //             ++joinCount;
    //         }
    //     }
    // }
    // end_time = MPI_Wtime();
    // cout<<"Process "<<id<<" basic join count: "<<joinCount<<" and time: "<<end_time - start_time<<"\n";

    // Count number of joined elements with btree
    double start_time = MPI_Wtime();
    int joinCount = 0;
    for (int i = 0; i < Asize; ++i) {
        // Get lower and upper bound of elements with the same x column
        Row lower_bound = {recvA[i].x, std::numeric_limits<uint>::min()};
        Row upper_bound = {recvA[i].x, std::numeric_limits<uint>::max()};
        auto lower = Bset.lower_bound(lower_bound);
        auto upper = Bset.upper_bound(upper_bound);
        for (auto it = lower; it != upper; ++it) {
            joinCount++;
        }   
    }
    end_time = MPI_Wtime();
    // cout<<"Process "<<id<<" set join count: "<<joinCount<<" and time: "<<end_time - start_time<<"\n";

    // Get the total join count
    int sumJoinCount;
    MPI_Reduce(&joinCount, &sumJoinCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (id == 0) {
        cout << "Join count: " << sumJoinCount << "\n";
    }

    // Get the total time 
    double localTime = end_time - start_time_set;
    double maxTime;
    MPI_Reduce(&localTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (id == 0) {
        cout << "Join time: " << maxTime << "\n";
    }

    MPI_Finalize ();
    return 0;
}