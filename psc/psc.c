#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

////////////////////// F U N C T I O N  P R O T O T Y P E S ////////////////////
void print_usage();
void make_full_path_pdb(char *pdb_file, char *data_dir, char *prot);
int count_ca(char *pdb_file, char *chn);
int use_line(char *line, int *ca_count, char *chn);
void read_ca(char *pdb_file, char *chn, float (*prot_coords)[3]);
void fill_euclid_dist(float (*prot1_coords)[3], float (*prot2_coords)[3],
	float **dist, int *prot1_len, int *prot2_len);
void fill_score(float **dist, float **score, int *prot1_len, int *prot2_len, 
	int *lmin);
////////////////////////////////////////////////////////////////////////////////

#define FNAME_LEN 256
#define NDIM 3
//#define LOG_DEBUG

///////////////////////////// M A C R O S //////////////////////////////////////
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	char *data_dir, *prot1, *prot2, chn1, chn2;
	char pdb_file1[FNAME_LEN], pdb_file2[FNAME_LEN];
	int prot1_len, prot2_len, iter, iter2, min_len, max_aln_len;
	float (*prot1_coords)[NDIM], (*prot2_coords)[NDIM];
	float **dist, **score;
	int **alignment;

	// parse parameters (we only support positional parameters)
	switch(argc) {
		case 6:
			data_dir = argv[1];
			prot1 = argv[2];
			chn1 = *argv[3];
			prot2 = argv[4];
			chn2 = *argv[5];
			break;
		default:
			print_usage();
	}

	// full path to pdb file
	make_full_path_pdb(pdb_file1, data_dir, prot1);
	make_full_path_pdb(pdb_file2, data_dir, prot2);

	// count number of alpha-carbons of given chain type
	prot1_len = count_ca(pdb_file1, &chn1);	
	prot2_len = count_ca(pdb_file2, &chn2);
	min_len = min(prot1_len, prot2_len);
	max_aln_len = min_len;	
	printf("%s:%d, %s:%d\n", prot1, prot1_len, prot2, prot2_len);

	// read ca coordinates from pdb files
	prot1_coords = malloc(prot1_len * NDIM * sizeof(float));
	prot2_coords = malloc(prot2_len * NDIM * sizeof(float));
	read_ca(pdb_file1, &chn1, prot1_coords);
	read_ca(pdb_file2, &chn2, prot2_coords);

	// prepare data structures
	dist = (float**) malloc(sizeof(float*) * prot1_len);
	for(iter = 0; iter < prot1_len; iter++)
		dist[iter] = (float*) malloc(sizeof(float) * prot2_len);
	alignment = (int**) malloc(sizeof(int*) * max_aln_len);
	for(iter = 0; iter < max_aln_len; iter++)
		alignment[iter] = (int*) malloc(sizeof(int) * 2);
	score = (float**) malloc(sizeof(float*) * prot1_len);
	for(iter = 0; iter < prot1_len; iter++)
		score[iter] = (float*) malloc(sizeof(float) * prot2_len);

	//
	fill_euclid_dist(prot1_coords, prot2_coords, dist, &prot1_len, &prot2_len);
	fill_score(dist, score, &prot1_len, &prot2_len, &min_len);

	// free memory
	free(prot1_coords);
	free(prot2_coords);
	for(iter = 0; iter < prot1_len; iter++)
		free(dist[iter]);
	free(dist);
	for(iter = 0; iter < prot1_len; iter++)
		free(score[iter]);
	free(score);
	for(iter = 0; iter < max_aln_len; iter++)
		free(alignment[iter]);
	free(alignment);

	// done
	return 0;
}

//////////////////////////// F U N C T I O N S /////////////////////////////////
void print_usage() {
	printf("Usage: psc [data_dir] [pdb_1] [chain_1] [pdb_2] [chain_2]\n");
	exit(-1);
}

void make_full_path_pdb(char *pdb_file, char *data_dir, char *prot) {
	strcpy(pdb_file, "");
	strcat(pdb_file, data_dir);
	strcat(pdb_file, "/");
	strcat(pdb_file, prot);
	strcat(pdb_file, ".pdb");
}

int use_line(char *line, int *ca_count, char *chn) {
	// tmalign style row selection
	char buffer1[3], buffer2[4], buffer3;
	strncpy(buffer1, line, 3);
	strncpy(buffer2, line + 12, 4);
	strncpy(&buffer3, line + 16, 1);
	if (strncmp(buffer1, "TER", 3) == 0 && ca_count - 1 > 0) return 1;
	if (strncmp(buffer1, "ATO", 3) != 0) return 2;
	if (strncmp(buffer2, "CA  ", 4) != 0 && 
		strncmp(buffer2, " CA ", 4) != 0 && 
		strncmp(buffer2, "  CA", 4) != 0)  return 2;
	if (buffer3 != *chn && buffer3 != ' ') return 2;
	return 0;
}

int count_ca(char *pdf_file, char *chn) {
	FILE *file;
	int ca_count = 0, is_usable;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	//
	file = fopen(pdf_file, "r");
	while((read = getline(&line, &len, file)) != -1) {
		is_usable = use_line(line, &ca_count, chn);
		if(is_usable == 1) break;
		else if(is_usable == 2) continue;
		ca_count++;
	}
	fclose(file);
	return ca_count;
}

void read_ca(char *pdb_file, char *chn, float (*prot_coords)[3]) {
	FILE *file;
	int ca_count = 0, is_usable;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	float x, y, z;
	int t;
	char chain;
	char aaname[3];
	//
	file = fopen(pdb_file, "r");
	while((read = getline(&line, &len, file)) != -1) {
		is_usable = use_line(line, &ca_count, chn);
		if(is_usable == 1) break;
		else if(is_usable == 2) continue;
		sscanf(&line[17], "%s %c%d", aaname, &chain, &t);
		sscanf(&line[30], "%f%f%f", &x, &y, &z);
		prot_coords[ca_count][0] = x;
		prot_coords[ca_count][1] = y;
		prot_coords[ca_count][2] = z;
		ca_count++;
	}
	fclose(file);
}

void fill_euclid_dist(float (*prot1_coords)[3], float (*prot2_coords)[3],
	float **dist, int *prot1_len, int *prot2_len) {
	int iter1, iter2;
	float x1, y1, z1, x2, y2, z2;
	for(iter1 = 0; iter1 < *prot1_len; iter1++) {
		x1 = prot1_coords[iter1][0];
		y1 = prot1_coords[iter1][1];
		z1 = prot1_coords[iter1][2];
		for(iter2 = 0; iter2 < *prot2_len; iter2++) {
			x2 = prot2_coords[iter2][0];
			y2 = prot2_coords[iter2][1];
			z2 = prot2_coords[iter2][2];			
			dist[iter1][iter2] = sqrt((pow(x2 - x1, 2) + 
				pow(y2 - y1, 2) + pow(z2 - z1, 2)));
		}
	}
#ifdef LOG_DEBUG	
	printf("dist\n");
	for(iter1 = 0; iter1 < *prot1_len; iter1++) {
		for(iter2 = 0; iter2 < *prot2_len; iter2++) {
			printf("%f ", dist[iter1][iter2]);
		}
		printf("\n");
	}
#endif
}

void fill_score(float **dist, float **score, int *prot1_len, 
	int *prot2_len, int *lmin) {
	int iter1, iter2;
	float d0Lmin;
	//
	d0Lmin = 1.24 * pow(*lmin - 1.5, (1/3)) - 1.8;
	for(iter1 = 0; iter1 < *prot1_len; iter1++) {
		for(iter2 = 0; iter2 < *prot2_len; iter2++) {
			score[iter1][iter2] = 1 / (1 + pow(dist[iter1][iter2]/d0Lmin, 2));
		}
	}
#ifdef LOG_DEBUG
	printf("score\n");
	for(iter1 = 0; iter1 < *prot1_len; iter1++) {
		for(iter2 = 0; iter2 < *prot2_len; iter2++) {
			printf("%f ", score[iter1][iter2]);
		}
		printf("\n");
	}
#endif
}
////////////////////////////////////////////////////////////////////////////////

