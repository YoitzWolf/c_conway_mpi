// (c) Razmyslov K. 2025 // Conways Game Of Life with MPI

#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>


// This is unsigned char in a fact :)
typedef uint8_t byte;

// send-recieve neighbours' borderlines
void send_recv_neighbours(int rank, int size, byte** field, size_t my_N, size_t my_M) {
	// printf(">> [%d] send recieve..\n", rank);
	byte* buff = malloc(sizeof(byte) * my_M );
	memcpy(buff, field[my_N-2], sizeof(byte) * my_M); // send bottom border
	byte* burr = malloc(sizeof(byte) * my_M );
	MPI_Status status;
	MPI_Sendrecv(
		  buff, my_M, MPI_BYTE, (rank+1)%size,  1,
                  burr, my_M, MPI_BYTE, (rank-1+size)%size, 1,
                  MPI_COMM_WORLD, &status );
	memcpy(field[0], burr, sizeof(byte)*my_M); // copy to upper border
	// ----------------------------------------------------------
	//printf(">> [%d] send recieve HALF DONE..\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);
	// ----------------------------------------------------------
	memcpy(buff, field[1], sizeof(byte) * my_M); // send upper bottom
	MPI_Sendrecv(
		  buff, my_M, MPI_BYTE, (rank-1+size)%size, 1,
                  burr, my_M, MPI_BYTE, (rank+1)%size, 1,
                  MPI_COMM_WORLD, &status );
	memcpy(field[my_N-1], burr, sizeof(byte)*my_M); // copy to upper border
	//printf(">> [%d] send recieve DONE!..\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);
}

// self field calculation tick
void calculate_self(int rank, int size, byte** field, size_t my_N, size_t my_M) {
	// printf(">> [%d] calculate ..\n", rank);
	byte **newfield = malloc(sizeof(byte*) * my_N);
        for (size_t r = 0; r < my_N; r++) {
        	newfield[r] = malloc(sizeof(byte) * my_M);
		memcpy(newfield[r], field[r], sizeof(byte) * my_M);
        }
	int counter = 0;
	size_t top_idx = 0;
	size_t bot_idx = 0;
	size_t lef_idx = 0;
	size_t rig_idx = 0;
	//printf(">> [%d] temp field done..\n", rank);
	for(size_t i=1; i<=my_N-2;i++) {
		for(size_t j=0; j<my_M;j++) {
			// calculate cords of nearest cells
			top_idx = i - 1;
			bot_idx = i + 1;
			lef_idx = (j - 1 + my_M) % my_M; // make tube
			rig_idx = (j + 1 + my_M) % my_M; // make tube
			counter = (
				field[top_idx][lef_idx] + field[top_idx][j] + field[top_idx][rig_idx] +
				field[i][lef_idx] + field[i][rig_idx] +
				field[bot_idx][lef_idx] + field[bot_idx][j] + field[bot_idx][rig_idx]
			);
			if (field[i][j] == 0) {
				// dead cell
				if (counter==3) {
					newfield[i][j]=1; // born
				}
			} else {
				// alive cell
				if (counter!=3 && counter!=2) {
					newfield[i][j]=0; // die
				}
			}
		}
	}
	// printf(">> [%d] copy newfield..\n", rank);
	for (size_t r = 1; r < my_N-1; r++) {
                // newfield[r] = malloc(sizeof(byte) * my_M);
                memcpy(field[r], newfield[r], sizeof(byte) * my_M);
        }
	//printf(">> [%d] calculate DONE..\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);
}

// root printing function
// awaits for other nodes to call MPI_Send
void root_print_field(int rank, int size, size_t N, size_t M, byte** field, size_t my_N, size_t my_M) {
	byte** world = malloc(N*sizeof(byte*));
	if (rank==0) {
		for(int i=0; i<N; i++){
			world[i] = malloc(sizeof(byte)*M);
		}
		for(int i=0; i<my_N-2; i++) {
			memcpy(world[i], field[i+1], sizeof(byte)*my_M);
		}
		MPI_Status status;
		for(int r=1; r<size; r++) {
			size_t rN = N / size;
		        if ( N % size != 0 ) {
	                	if (size - 1 == r) {
        		                rN -= (size-1);
		                } else {
                		        rN += 1;
        	        	}
			}
			for(int i=0; i<rN; i++) {
				int tag = r*rN + i; // tag is line number in global field
				printf(">> [%d] recieve tag %d ..\n", r, tag);
				MPI_Recv (
					world[tag],
					M,
					MPI_BYTE,
					r,
					tag,
					MPI_COMM_WORLD,
					&status
				);
			}
		}
	}
	if(rank == 0) {
		// print
		for(int i=0; i<N; i++) {
			printf("%d\t\t|", i+1);
			for(int j=0; j<M; j++){
				char c;
				if (world[i][j] == 1) {
					c = '#';
				} else {
					c = '*';
				}
				printf("%c", c);
			}
			printf("\n");
		}
	} else {
		// IGNORE
	}
}

void print_field(int tick, int rank, int size, size_t N, size_t M, byte** field, size_t my_N, size_t my_M)  {
	if (rank==0) {
                printf("\n----------------------------\nDumpingResult (%d)\n----------------------\n", tick);
                root_print_field(rank, size, N, M, field, my_N, my_M);
        } else {
//              MPI_Status status;
                if (size > 1) {
                        for(int i=1; i<=my_N-2; i++) {
                                int tag = rank*(my_N-2) + i - 1; // tag is line numer in global field
                                // printf(">> [%d] send field[%d] with tag %d ..\n", rank, i, tag);
                                MPI_Send(
                                        field[i],
                                        my_M,
                                        MPI_BYTE,
                                        0,
                                        tag,
                                        MPI_COMM_WORLD
                                );
                        }
                }
        }
}

int main( int argn, char **argv ){
  	int size, rank;
  	(void)MPI_Init( &argn, &argv );
  	(void)MPI_Comm_size( MPI_COMM_WORLD, &size );
  	(void)MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  	char nodename[MPI_MAX_PROCESSOR_NAME];
  	int   nodesize;
  	(void)MPI_Get_processor_name((char *) &nodename, &nodesize);
	size_t N, M;
	int T = 7777;
	size_t Freq;
	if (rank==0) {
		// N = 100;
		// M = 100;
		printf("enter [M x N, T, F] field size, ticks to calculate and output frequency: ");
		scanf("%d%d%d%d", &N, &M, &T, &Freq);
	}
	(void)MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD); // send toa ll
	(void)MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD); // send to all
	size_t my_N = N;
	size_t my_M = M;
	size_t nAdd = 2;
	bool rotated = 0;
	if (N < M) {
		// rotate
		rotated = 1;
		size_t b = my_M;
		my_M = my_N;
		M = N;
		my_N = b;
		N = b;
	}
	//cutoff my_N:
	my_N = my_N / size; // <- self size calculation
	if ( N % size != 0 ) { // check if need to chenge node sizes
		if (size - 1 == rank) {
			my_N -= (size-1); //  <- size modification for last node
		} else {
			my_N += 1; // size modification for first node
		}
	}
	// --------------------------------------------------
	// setup field
	// --------------------------------------------------
	MPI_Barrier(MPI_COMM_WORLD);
	printf(">> [%d] field: %d x %d \n", rank, my_N, my_M);
	my_N += 2;
	// my_M;
	// byte** field = malloc(fN*fM*sizeof(byte)); // field is [1:my_N+1] of [0:M]
        // other are borders

	byte **field = malloc(sizeof(byte*) * my_N);
	for (size_t r = 0; r < my_N; r++) {
	  field[r] = malloc(sizeof(byte) * my_M);
	}
	// -------------------------------------------------
	// SETUP FIELD
	// -------------------------------------------------
	if (rank==0) {
		// field[1][2] = 1;
		printf("\nFIELD ROTATED: %d\n", rotated);
		// GLIDER !!!
		field[2][3] = 1;
		field[3][3] = 1;
		field[4][3] = 1;
		field[4][2] = 1;
		field[3][1] = 1;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// -------------------------------------------------
	// MAIN CYCLE
	// -------------------------------------------------
	for (int t=0; t<T; t++) {
		if (size > 1) {
			send_recv_neighbours(rank, size, field, my_N, my_M);
		} else {
			// Do if need self-closure
			// printf("Self-closure");
			// if size is 1 -> only one cell -> need to make torus
			// memcpy 1->my_N+1
			memcpy(field[my_N-1], field[1], sizeof(byte) * my_M);
			// memcpy my_N->0
			memcpy(field[0], field[my_N-2], sizeof(byte) * my_M);
		}
		calculate_self(rank, size, field, my_N, my_M);
		MPI_Barrier(MPI_COMM_WORLD);
		if (t % Freq == 0) {
			print_field(t, rank, size, N, M, field, my_N, my_M);
		}
	}
	// -----------------------------------------------

	if (rank==0) {
		printf("\n---------------------------------\nEnd sumulation\n-----------------------------\n");
	}
	print_field(T-1, rank, size, N, M, field, my_N, my_M);

	MPI_Finalize();
	return 0;
}

