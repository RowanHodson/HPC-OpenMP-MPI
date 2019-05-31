// Original Author: Daniel Shao
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <sys/time.h>
#include <unistd.h>
#include "mpi.h"
#include <stdbool.h>

int i,j,k;
int N = 1000;


int cmp(const void * a, const void * b) {
  if (*(int*)a < *(int*)b) return -1;
  if (*(int*)a > *(int*)b) return 1;
  else return 0;
}

void phase1(int *array, int N, int startIndex, int subArraySize, int *pivots, int p) {
  // local sorting of subarrays
  qsort(array + startIndex, subArraySize, sizeof(array[0]), cmp);

  // regular sampling
  for (i = 0; i < p; i++) {
    pivots[i] = array[startIndex + (i * (N / (p * p)))];    
  }
  return;
}

void phase2(int *array, int startIndex, int subArraySize, int *pivots, int *partitionSizes, int p, int myId) {
  int *collectedPivots = (int *) malloc(p * p * sizeof(pivots[0]));
  int *phase2Pivots = (int *) malloc((p - 1) * sizeof(pivots[0]));          //主元
  int index = 0;

   // Collecting messages, the root process contains connections to the send buffers of all processes in its accept buffer.
  MPI_Gather(pivots, p, MPI_INT, collectedPivots, p, MPI_INT, 0, MPI_COMM_WORLD);       
  if (myId == 0) {

    qsort(collectedPivots, p * p, sizeof(pivots[0]), cmp);          //Sorting regular sampled samples

    // Principal selection after sampling and sorting
    for (i = 0; i < (p -1); i++) {
      phase2Pivots[i] = collectedPivots[(((i+1) * p) + (p / 2)) - 1];
    }
  }
  //Send broadcast
  MPI_Bcast(phase2Pivots, p - 1, MPI_INT, 0, MPI_COMM_WORLD);
  // Perform the principal division and calculate the size of the division
  for ( i = 0; i < subArraySize; i++) {
    if (array[startIndex + i] > phase2Pivots[index]) {
      //If the size of the current position exceeds the pivot position, proceed to the next division
      index += 1;
    }
    if (index == p) {
      //The last division, the total length of the sub-array minus the current position can get the size of the last sub-array division
      partitionSizes[p - 1] = subArraySize - i + 1;
      break;
    }
    partitionSizes[index]++ ;   // Divided size increase
  }
  free(collectedPivots);
  free(phase2Pivots);
  return;
}

void phase3(int *array, int startIndex, int *partitionSizes, int **newPartitions, int *newPartitionSizes, int p) {
  int totalSize = 0;
  int *sendDisp = (int *) malloc(p * sizeof(int));
  int *recvDisp = (int *) malloc(p * sizeof(int));
    
  //Global to global delivery, each process can send a different amount of data to each recipient.
  MPI_Alltoall(partitionSizes, 1, MPI_INT, newPartitionSizes, 1, MPI_INT, MPI_COMM_WORLD);

  // Calculate the total size of the partition and allocate space for the new partition
  for ( i = 0; i < p; i++) {
    totalSize += newPartitionSizes[i];
  }
  *newPartitions = (int *) malloc(totalSize * sizeof(int));

  // Calculate the displacement relative to sendbuf before sending the partition, this shift stores the data output to the process
  sendDisp[0] = 0;
  recvDisp[0] = 0;      //计算相对于recvbuf的位移，此位移处存放着从进程接受到的数据
  for ( i = 1; i < p; i++) {
    sendDisp[i] = partitionSizes[i - 1] + sendDisp[i - 1];
    recvDisp[i] = newPartitionSizes[i - 1] + recvDisp[i - 1];
  }

  //发送数据，实现n次点对点通信
  MPI_Alltoallv(&(array[startIndex]), partitionSizes, sendDisp, MPI_INT, *newPartitions, newPartitionSizes, recvDisp, MPI_INT, MPI_COMM_WORLD);

  free(sendDisp);
  free(recvDisp);
  return;
}

void phase4(int *partitions, int *partitionSizes, int p, int myId, int *array) {
  int *sortedSubList;
  int *recvDisp, *indexes, *partitionEnds, *subListSizes, totalListSize;

  indexes = (int *) malloc(p * sizeof(int));
  partitionEnds = (int *) malloc(p * sizeof(int));
  indexes[0] = 0;
  totalListSize = partitionSizes[0];
  for ( i = 1; i < p; i++) {
    totalListSize += partitionSizes[i];
    indexes[i] = indexes[i-1] + partitionSizes[i-1];
    partitionEnds[i-1] = indexes[i];
  }
  partitionEnds[p - 1] = totalListSize;

  sortedSubList = (int *) malloc(totalListSize * sizeof(int));
  subListSizes = (int *) malloc(p * sizeof(int));
  recvDisp = (int *) malloc(p * sizeof(int));

  // 归并排序
  for ( i = 0; i < totalListSize; i++) {
    int lowest = INT_MAX;
    int ind = -1;
    for (j = 0; j < p; j++) {
      if ((indexes[j] < partitionEnds[j]) && (partitions[indexes[j]] < lowest)) {
    lowest = partitions[indexes[j]];
    ind = j;
      }
    }
    sortedSubList[i] = lowest;
    indexes[ind] += 1;
  }

  // 发送各子列表的大小回根进程中
  MPI_Gather(&totalListSize, 1, MPI_INT, subListSizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // 计算根进程上的相对于recvbuf的偏移量
  if (myId == 0) {
    recvDisp[0] = 0;
    for ( i = 1; i < p; i++) {
      recvDisp[i] = subListSizes[i - 1] + recvDisp[i - 1];
    }
  }

  //发送各排好序的子列表回根进程中
  MPI_Gatherv(sortedSubList, totalListSize, MPI_INT, array, subListSizes, recvDisp, MPI_INT, 0, MPI_COMM_WORLD);

  free(partitionEnds);
  free(sortedSubList);
  free(indexes);
  free(subListSizes);
  free(recvDisp);
  return;
}

//PSRS排序函数，调用了4个过程函数
bool psrs_mpi(int *array, int N)    
{
    int p, myId, *partitionSizes, *newPartitionSizes, nameLength;
    int subArraySize, startIndex, endIndex, *pivots, *newPartitions;
    char processorName[MPI_MAX_PROCESSOR_NAME];


    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&myId);
    MPI_Get_processor_name(processorName,&nameLength);

    printf("Process %d is on %s\n",myId, processorName);

    pivots = (int *) malloc(p*sizeof(int));
    partitionSizes = (int *) malloc(p*sizeof(int));
    newPartitionSizes = (int *) malloc(p*sizeof(int));
    for ( k = 0; k < p; k++) {
      partitionSizes[k] = 0;
    }

    // 获取起始位置和子数组大小
    startIndex = myId * N / p;
    if (p == (myId + 1)) {
      endIndex = N;
    } 
    else {
      endIndex = (myId + 1) * N / p;
    }
    subArraySize = endIndex - startIndex;

    MPI_Barrier(MPI_COMM_WORLD);
    //调用各阶段函数
    phase1(array, N, startIndex, subArraySize, pivots, p);
    if (p > 1) {
      phase2(array, startIndex, subArraySize, pivots, partitionSizes, p, myId);
      phase3(array, startIndex, partitionSizes, &newPartitions, newPartitionSizes, p);
      phase4(newPartitions, newPartitionSizes, p, myId, array);
    }

    if (myId == 0) 
     for(k = 1; k < N; k++){
        if (array[k] < array[k-1])
        {
          printf("\n");
          if (p > 1) {
            free(newPartitions);
        }
         free(partitionSizes);
        free(newPartitionSizes);
         free(pivots);


          free(array);
           MPI_Finalize();
          return false;
        }
     }
     
     printf("\n");
    if (p > 1) {
      free(newPartitions);
    }
    free(partitionSizes);
    free(newPartitionSizes);
    free(pivots);


  free(array);
  MPI_Finalize();
  return true;
}

int main(int argc, char *argv[]) {
  double etime;
  int *array;
  bool sorting = true;
  if (argv[1] != NULL)
  {
    N = atoi(argv[1]);
  }
  array = (int *) malloc(N*sizeof(int));

    //srand(N);
    for ( k = 0; k < N; k++) {
      array[k] = rand();
    }
    etime = - MPI_Wtime();
    MPI_Init(&argc,&argv);      
    bool result = psrs_mpi(array,N);
    etime += MPI_Wtime();
    if(result == true)
    {
       printf("\n Sorted successfully in %f seconds", etime);
    }
    else
    {
      printf("FAILED");
    }
  return 0;
}
