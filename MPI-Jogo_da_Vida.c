#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>

#define N 2048 
#define IT_MAX 2000

// faz alocacao dinamica da matriz
int** constroiMatriz(int tam){
    int i, j;
    int** matriz = (int** )malloc(tam * sizeof(int*));
    
    for(int i = 0;i<tam;i++){
        matriz[i] = (int*)malloc((N+2) * sizeof(int));
    }

    for(i=0;i<tam;i++){
        for(j=0;j<N+2;j++){
            matriz[i][j] = 0;
        }
    }
    
    return matriz;
}


// inicializa o tabuleiro
void tabuleiroInicial(int** grid){
    int lin = 2, col = 2;
    //inicializa o Glider
    grid[lin  ][col+1] = 1;
    grid[lin+1][col+2] = 1;
    grid[lin+2][col  ] = 1;
    grid[lin+2][col+1] = 1;
    grid[lin+2][col+2] = 1;
    //-------------------
    lin = 11;
    col = 31;
    //inicializa o R-pentonimo
    grid[lin  ][col+1] = 1;
    grid[lin  ][col+2] = 1;
    grid[lin+1][col  ] = 1;
    grid[lin+1][col+1] = 1;
    grid[lin+2][col+1] = 1;
    //-------------------
}


// calcula a quantidade de vizinhos
int getNeighbors(int i, int j, int **grid){
    return grid[i+1][j-1] + grid[i+1][j] + grid[i+1][j+1]
         + grid[i-1][j-1] + grid[i-1][j] + grid[i-1][j+1]
         + grid[i][j-1] + grid[i][j+1] ;
}


// determina o novo status da célula
int defineStatus(int neighbors, int currentStatus){

    if(currentStatus == 1){
        if(neighbors == 2 || neighbors == 3){
            return 1;
        }
    } else {
        if(neighbors == 3){
            return 1;
        }
    }

    return 0;
}


int main(int argc, char * argv[]) {

    /*Buffers de comunicação*/
    int recA[N], recB[N]; 
    int sendA[N], sendB[N];
    
    // limites do tabuleiro
    int end=0, start=0;

    double time1, time2, duration;

    int rankProcess, numberProcess;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numberProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankProcess);

    int beforeProcess = ((rankProcess+numberProcess-1)%numberProcess);
    int afterProcess = ((rankProcess+1)%numberProcess);

    int slice = (N/numberProcess);

    start = rankProcess*slice;

    if((start+slice) >= N){
        slice = N-start+1;
    }

    int** grid = constroiMatriz(slice+2);
    int** newGrid = constroiMatriz(slice+2);

    if(rankProcess == 0){
        tabuleiroInicial(grid);
    }
    
    int neighbors = 0, liveCels = 0;

    time1 = MPI_Wtime();


    for(int k = 0; k < IT_MAX; k++){
        
        for(int j = 0; j < N;j++){
            sendA[j] = grid[1][j+1];
            sendB[j] = grid[slice][j+1];
        }

        MPI_Sendrecv (sendB, N, MPI_INTEGER, afterProcess, 1,
                    recA, N, MPI_INTEGER, beforeProcess, 1,
                    MPI_COMM_WORLD, &status);

        MPI_Sendrecv (sendA, N, MPI_INTEGER, beforeProcess, 1,
                    recB,  N, MPI_INTEGER, afterProcess, 1,
                    MPI_COMM_WORLD, &status);
        
        // copia as linhas recebidas para as linhas de sangria
        for(int j = 0; j < N;j++){
            grid[0][j+1] = recA[j];
            grid[slice+1][j+1] = recB[j];
        }

        // torna o tabuleiro infinito horizontalmente
        for(int i = 0; i < slice+2; i++){
            grid[i][N+1] = grid[i][1];
            grid[i][0] = grid[i][N];
        }

        // gera o novo tabuleiro
        for(int i = 1; i < slice+1; i++){
            for(int j = 1; j < N+1;j++){
                neighbors = getNeighbors(i,j,grid);
                newGrid[i][j] = defineStatus(neighbors, grid[i][j]);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // atualiza a grid e calcula a quanntidade de vivos
        liveCels = 0;
        for(int i = 1; i < slice+1; i++) {
            for(int j = 1; j < N+1;j++) {
                grid[i][j] = newGrid[i][j];
                liveCels += newGrid[i][j];
            }
        } 
    }

    if(rankProcess == 0){

        int results[numberProcess], sum = 0;
        results[0] = liveCels;
        
        for(int i = 1; i < numberProcess; i++){
            MPI_Recv(&results[i], 1, MPI_INTEGER,i,10, MPI_COMM_WORLD, &status);
        }

        for(int i = 0; i < numberProcess; i++){
            sum += results[i];
        }

        printf("\nTotal de celulas vivas: %d\n", sum);

    } else {

        MPI_Send(&liveCels, 1, MPI_INTEGER,0,10, MPI_COMM_WORLD); 
    }

    if(rankProcess == 0){
        time2 = MPI_Wtime();
        duration = time2 - time1;
        printf("Runtime at %d is %f \n", rankProcess,duration);
    }
    
    MPI_Finalize();

    return 0;
}