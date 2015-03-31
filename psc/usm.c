#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
//#include <sys/timeb.h>
#include <vector>
#include <string>
#include <sstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/unordered_map.hpp>

#include <boost/iostreams/filter/counter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/filesystem.hpp>


using namespace std;
using namespace boost::iostreams;
namespace io = boost::iostreams;

////////////////////// F U N C T I O N  P R O T O T Y P E S ////////////////////
void print_usage();
void make_full_path_pdb(char *pdb_file, char *data_dir, char *prot);
int count_ca(char *pdb_file, char *chn);
int use_line(char *line, int *ca_count, char *chn);
void read_ca(char *pdb_file, char *chn, float **prot_coords);
void fill_euclid_dist(float **prot1_coords, float **prot2_coords,
	float **dist, int *prot1_len, int *prot2_len);
int fill_contact_map(float **dist, int *prot_len, int **cm);
void print_contact_map(int **cm, int *cm_len);
double calculate_pairwise_distances(std::vector< std::string > &vector1, 
        std::vector< std::string > &vector2);
double compute_distance(int x, int y, int xy, int yx);
////////////////////////////////////////////////////////////////////////////////

#define FNAME_LEN 256
#define NDIM 3
//#define LOG_DEBUG
#define MAX_ATOMS 1000
#define MAX_DIST 6.5

///////////////////////////// M A C R O S //////////////////////////////////////
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
////////////////////////////////////////////////////////////////////////////////

////////////////////////// T E M P L A T E S ///////////////////////////////////
template<typename T>
inline int compressIt(std::vector<T> s){
    std::stringstream uncompressed, compressed;
    for (typename std::vector<T>::iterator it = s.begin();
         it != s.end(); it++)
        uncompressed << *it;

    io::filtering_streambuf<io::input> o;
    o.push(io::gzip_compressor());
    o.push(uncompressed);
    io::copy(o, compressed);

    return compressed.str().length();
}
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	char *data_dir, *prot1, *prot2, chn1, chn2;
	char pdb_file1[FNAME_LEN], pdb_file2[FNAME_LEN];
	int prot1_len, prot2_len, iter, iter2, min_len, max_aln_len;
	float **prot1_coords, **prot2_coords;
	int **cm1, **cm2;
	float **dist1, **dist2;
	int cm1_len, cm2_len;

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

	// prepare data structures
	prot1_coords = (float**) malloc(sizeof(float*) * prot1_len);
	for(iter = 0; iter < prot1_len; iter++)
		prot1_coords[iter] = (float*) malloc(sizeof(float) * NDIM);
	prot2_coords = (float**) malloc(sizeof(float*) * prot2_len);
	for(iter = 0; iter < prot2_len; iter++)
		prot2_coords[iter] = (float*) malloc(sizeof(float) * NDIM);
	dist1 = (float**) malloc(sizeof(float*) * prot1_len);
	for(iter = 0; iter < prot1_len; iter++)
		dist1[iter] = (float*) malloc(sizeof(float) * prot1_len);
	dist2 = (float**) malloc(sizeof(float*) * prot2_len);
	for(iter = 0; iter < prot2_len; iter++)
		dist2[iter] = (float*) malloc(sizeof(float) * prot2_len);
	cm1 = (int**) malloc(sizeof(int*) * prot1_len * prot1_len / 2);
	for(iter = 0; iter < prot1_len * prot1_len / 2; iter++)
		cm1[iter] = (int*) malloc(sizeof(int) * 2);
	cm2 = (int**) malloc(sizeof(int*) * prot2_len * prot2_len / 2);
	for(iter = 0; iter < prot2_len * prot2_len / 2; iter++)
		cm2[iter] = (int*) malloc(sizeof(int) * 2);

	// read ca coordinates from pdb files
	read_ca(pdb_file1, &chn1, prot1_coords);
	read_ca(pdb_file2, &chn2, prot2_coords);

	// get intra-protein inter-residue distances
	fill_euclid_dist(prot1_coords, prot1_coords, dist1, &prot1_len, &prot1_len);
	fill_euclid_dist(prot2_coords, prot2_coords, dist2, &prot2_len, &prot2_len);	
	
	// get contact maps
	cm1_len = fill_contact_map(dist1, &prot1_len, cm1);
	cm2_len = fill_contact_map(dist2, &prot2_len, cm2);
	printf("cm lens: %d %d\n",cm1_len,cm2_len);
#ifdef LOG_DEBUG	
	print_contact_map(cm1, &cm1_len);
#endif

	// get usm distance
	std::vector<std::string> cm1Vec;
	for(iter = 0; iter < cm1_len; iter++) {
		stringstream ss;
		ss << cm1[iter][0] << " " << cm1[iter][1];
		cm1Vec.push_back(ss.str());
	}	
	std::vector<std::string> cm2Vec;
	for(iter = 0; iter < cm2_len; iter++) {
		stringstream ss;
		ss << cm2[iter][0] << " " << cm2[iter][1];
		cm2Vec.push_back(ss.str());
	}	
	printf("usm: %f\n", calculate_pairwise_distances(cm1Vec, cm2Vec));
	
	// free memory
	for(iter = 0; iter < prot1_len; iter++) free(prot1_coords[iter]);
	free(prot1_coords);
	for(iter = 0; iter < prot2_len; iter++) free(prot2_coords[iter]);
	free(prot2_coords);
	for(iter = 0; iter < prot1_len; iter++) free(dist1[iter]);
	free(dist1);
	for(iter = 0; iter < prot2_len; iter++) free(dist2[iter]);
	free(dist2);
	for(iter = 0; iter < prot1_len; iter++) free(cm1[iter]);
	free(cm1);
	for(iter = 0; iter < prot2_len; iter++) free(cm2[iter]);
	free(cm2);

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

void read_ca(char *pdb_file, char *chn, float **prot_coords) {
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

void fill_euclid_dist(float **prot1_coords, float **prot2_coords,
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

int fill_contact_map(float **dist, int *prot_len, int **cm) {
	int i, j, idx = 0;
	int upper = min(*prot_len, MAX_ATOMS);
	for(i = 0; i < upper; i++) {
		for(j = i + 2; j < upper; j++) {
			if (dist[i][j] <= MAX_DIST) {
				cm[idx][0] = i;
				cm[idx++][1] = j;				
			}
		}
	}	
	return idx;
}

void print_contact_map(int **cm, int *cm_len) {
	int i;
	printf("CM\n");
	for(i = 0; i < *cm_len; i++) printf("\t%d %d\n",cm[i][0],cm[i][1]);
	printf("CM End\n");
}

double calculate_pairwise_distances(std::vector< std::string > &vector1, std::vector< std::string > &vector2) {
	int x = compressIt(vector1);
	int y = compressIt(vector2);
	
	std::vector< std::string > v1;
	v1.insert(v1.end(), vector1.begin(), vector1.end());
	v1.insert(v1.end(), vector2.begin(), vector2.end());	
	int xy = compressIt(v1);

	v1.clear();
        v1.insert(v1.end(), vector2.begin(), vector2.end());
	v1.insert(v1.end(), vector1.begin(), vector1.end());
	int yx = compressIt(vector2);

	return compute_distance(x, y, xy, yx);
}

double compute_distance(int x, int y, int xy, int yx) {
	double maxNum = ((yx - y) > (xy - x)) ? yx - y : xy - x;
	double maxDen = (x > y) ? x : y;

        return (maxNum/maxDen);
}
////////////////////////////////////////////////////////////////////////////////

