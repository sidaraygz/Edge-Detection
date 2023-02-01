#include<stdio.h>
#include<stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

typedef struct PGMData {
    char version[3];
    int row;
    int col;
    int max_gray;
    int **matrix;
    int **gx;
    int **gy;
} PGMData;


void init_out_image( PGMData *out, PGMData *data){
    int i, j;
    strcpy(out->version, data->version);
    out->col = data->col;
    out->row = data->row;
    out->max_gray = data->max_gray;

    out->matrix = (int**) calloc(out->row, sizeof(int*));
    for(i = 0; i < out->row; i++) {
        out->matrix[i] = calloc(out->col, sizeof(int));
    }

    out->gx = (int**) calloc(out->row, sizeof(int*));
    for(i = 0; i < out->row; i++) {
        out->gx[i] = calloc(out->col, sizeof(int));
    }

    out->gy = (int**) calloc(out->row, sizeof(int*));
    for(i = 0; i < out->row; i++) {
        out->gy[i] = calloc(out->col, sizeof(int));
    }

    for(i = 0; i < out->row; i++) {
        for(j = 0; j < out->col; j++) {
            out->matrix[i][j] = data->matrix[i][j];
            out->gx[i][j] = data->matrix[i][j];
            out->gy[i][j] = data->matrix[i][j];
        };
    }
}

int **allocate_dynamic_matrix_int(int row, int col)
{
    int **ret_val;
    int i;

    ret_val = (int **)malloc(sizeof(int *) * row);
    if (ret_val == NULL) {
        perror("memory allocation failure");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < row; ++i) {
        ret_val[i] = (int *)malloc(sizeof(int) * col);
        if (ret_val[i] == NULL) {
            perror("memory allocation failure");
            exit(EXIT_FAILURE);
        }
    }

    return ret_val;
}
double **allocate_dynamic_matrix_double(int row, int col)
{
    double **ret_val;
    int i;

    ret_val = (double **)malloc(sizeof(double *) * row);
    if (ret_val == NULL) {
        perror("memory allocation failure");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < row; ++i) {
        ret_val[i] = (double *)malloc(sizeof(double ) * col);
        if (ret_val[i] == NULL) {
            perror("memory allocation failure");
            exit(EXIT_FAILURE);
        }
    }

    return ret_val;
}


double **gaussKernelOlustur(int smooth_kernel_size, double  sigma){
    double **gauss = allocate_dynamic_matrix_double(smooth_kernel_size, smooth_kernel_size);
    double sum = 0;
    int i, j;

    for (i = 0; i < smooth_kernel_size; i++) {
        for (j = 0; j < smooth_kernel_size; j++) {
            double x = i - (smooth_kernel_size - 1) / 2.0;
            double y = j - (smooth_kernel_size - 1) / 2.0;
            gauss[i][j] = exp(((pow(x, 2) + pow(y, 2)) / ((2 * pow(sigma, 2)))) * (-1));
            sum += gauss[i][j];
        }
    }
    for (i = 0; i < smooth_kernel_size; i++) {
        for (j = 0; j < smooth_kernel_size; j++) {
            gauss[i][j] /= sum;
        }
    }
    printf("\n");
    for (i = 0; i < smooth_kernel_size; i++) {
        for (j = 0; j < smooth_kernel_size; j++) {
            printf("%f ", gauss[i][j]);
        }
        printf("\n");
    }
    return gauss;
}

void deallocate_dynamic_matrix(int **matrix, int row)
{
    int i;

    for (i = 0; i < row; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

void skipComments(FILE *fp)
{
    int ch;
    char line[100];
    while ((ch = fgetc(fp)) != EOF && isspace(ch)) {
        ;
    }

    if (ch == '#') {
        fgets(line, sizeof(line), fp);
        skipComments(fp);
    } else {
        fseek(fp, -1, SEEK_CUR);
    }
}
void readPGM(const char *file_name, PGMData *data) {
    FILE *pgmFile;
    int i, j;
    int lo, hi;
    pgmFile = fopen(file_name, "rb");
    if (pgmFile == NULL) {
        perror("cannot open file to read");
        exit(EXIT_FAILURE);
    }
    fgets(data->version, sizeof(data->version), pgmFile);
    skipComments(pgmFile);
    fscanf(pgmFile, "%d %d %d", &data->col, &data->row, &data->max_gray);
    data->matrix = allocate_dynamic_matrix_int(data->row, data->col);

    if (!strcmp(data->version, "P2")) {
        for (i = 0; i < data->row; i++) {
            for (j = 0; j < data->col; j++) {
                fscanf(pgmFile,"%d", &data->matrix[i][j]);
            }
        }
    }
    else if(!strcmp(data->version, "P5")){
        fgetc(pgmFile);
        if (data->max_gray > 255) {
            for (i = 0; i < data->row; ++i) {
                for (j = 0; j < data->col; ++j) {
                    hi = fgetc(pgmFile);
                    lo = fgetc(pgmFile);
                    data->matrix[i][j] = (hi << 8) + lo;
                }
            }
        }
        else {
            for (i = 0; i < data->row; ++i) {
                for (j = 0; j < data->col; ++j) {
                    lo = fgetc(pgmFile);
                    data->matrix[i][j] = lo;
                }
            }
        }
    }
    fclose(pgmFile);
}

void writePGM(const char *filename, const PGMData *data, int **matrix) {
    FILE *out_image;
    int i, j;
    char ver[3] = "P2";
    out_image = fopen(filename, "w");
    fprintf(out_image, "%s\n", ver);
    fprintf(out_image, "%d %d\n", data->col, data->row);
    fprintf(out_image, "%d\n", data->max_gray);


    for (i = 0; i < data->row; i++) {
        for (j = 0; j < data->col; j++) {
            fprintf(out_image, "%d  ", matrix[i][j]);
        }


        fprintf(out_image, "\n");
    }
    deallocate_dynamic_matrix(matrix, data->row);

    fclose(out_image);
}


void padding(PGMData *data, int kernelSize) {
    int i, j;
    for (i = 0; i < data->col; i++) {
        for (j = 0; j < kernelSize; j++) {
            data->matrix[j][i] = 0;
            data->matrix[data->row - (j+1)][i] = 0;
        }
    }

    for (i = 0; i < data->row; i++) {
        for(j = 0; j < kernelSize; j++) {
            data->matrix[i][j] = 0;
            data->matrix[i][data->col - (j+1)] = 0;
        }
    }
}

int convolution(PGMData* image, double **kernel, int row, int col, int m) {
    int i, j, sum = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            sum += image->matrix[i + row][j + col] * kernel[i][j];
        }
    }
    return sum;
}

void gauss(PGMData *outImg, PGMData *data, double **kernel, int m){
    int i, j, k;
    for ( i = m/2 ; i < data->row - (m-1); i++ ) {
        for (j = m/2; j < data->col - (m-1); j++) {
            k = convolution(data, kernel, i, j, m);
            outImg->matrix[i][j] = k;
        }
    }
}
void min_max_normalization(PGMData * image, int** matrix) {
    int min = 1000000, max = 0, i, j;

    for(i = 0; i < image->row; i++) {
        for(j = 0; j < image->col; j++) {
            if (matrix[i][j] < min) {
                min = matrix[i][j];
            }
            else if (matrix[i][j] > max) {
                max = matrix[i][j];
            }
        }
    }

    for(i = 0; i < image->row; i++) {
        for(j = 0; j < image->col; j++) {
            double ratio = (double) (matrix[i][j] - min) / (max - min);
            matrix[i][j] = ratio * 255;
        }
    }
}
void sobel_edge_detector(PGMData *outImg, PGMData *data, int m){
    int i, j, gx, gy;
    double **x = allocate_dynamic_matrix_double(m, m);
    double **y = allocate_dynamic_matrix_double(m, m);


    int mx[3][3] = {
            {-1, 0, 1},
            {-2, 0, 2},
            {-1, 0, 1}
    };
    int my[3][3] = {
            {-1, -2, -1},
            {0, 0, 0},
            {1, 2, 1}
    };
    for(i = 0; i<m; i++){
        for(j = 0; j<m; j++){
            x[i][j] = mx[i][j];
        }
    }
    for(i = 0; i<m; i++){
        for(j = 0; j<m; j++){
            y[i][j] = my[i][j];
        }
    }
    for (i = 1; i < data->row - 2; i++) {
        for (j = 1; j < data->col - 2; j++) {
            gx = convolution(data, x, i, j, m);
            gy = convolution(data,  y, i, j, m);
            outImg->matrix[i][j] = sqrt(gx * gx + gy * gy);;
            outImg->gx[i][j] = gx;
            outImg->gy[i][j] = gy;
        }
    }
    deallocate_dynamic_matrix((int **) x, 3);
    deallocate_dynamic_matrix((int **) y, 3 );

}

void laplace(PGMData *outImg, PGMData *data){
    int m = 3;
    int i, j;
    double **x = allocate_dynamic_matrix_double(m, m);
    double **y = allocate_dynamic_matrix_double(m, m);
    int mx[3][3] = {
            {0, -1, 0},
            {-1, 4, -1},
            {0, -1, 0}
    };
    int my[3][3] = {
            {-1, -1, -1},
            {-1, 8, -1},
            {-1, -1, -1}
    };
    for(i = 0; i<m; i++){
        for(j = 0; j<m; j++){
            x[i][j] = mx[i][j];
        }
    }
    for(i = 0; i<m; i++){
        for(j = 0; j<m; j++){
            y[i][j] = my[i][j];
        }
    }

    for (i = 1; i < data->row - 2; i++) {
        for (j = 1; j < data->col - 2; j++) {
            int k = convolution(data, x, i, j, m);
            int l = convolution(data, y, i, j, m);
            outImg->gx[i][j] = k;
            outImg->gy[i][j] = l;
        }
    }

}

int main() {
    char filename[500], outname[500];
    double **kernel;
    double sigma;
    PGMData *imginfo, *outImg;
    int kernelSize;
    imginfo = (PGMData *) malloc(sizeof(PGMData));
    outImg = (PGMData *) malloc(sizeof(PGMData));


    printf("Dosya adini giriniz: ");
    scanf("%s", filename);
    printf("Cikti dosyasinin adini giriniz: ");
    scanf("%s", outname);
    readPGM(filename, imginfo);
    printf("Kernel boyutunu giriniz:  ");
    scanf("%d", &kernelSize);
    padding(imginfo, kernelSize);
    init_out_image(outImg, imginfo);
    //laplace(outImg, imginfo);
    //sobel_edge_detector(outImg, imginfo, kernelSize);
    //min_max_normalization(outImg, outImg->gx);
    //min_max_normalization(outImg, outImg->gy);
    //writePGM(outname, outImg, outImg->gx);
    //writePGM(outname, outImg, outImg->gx);
    printf("Sigma degerini giriniz: ");
    scanf("%lf", &sigma);
    kernel = gaussKernelOlustur(kernelSize, sigma);
    gauss(outImg, imginfo, kernel, kernelSize);
    writePGM(outname, outImg, outImg->matrix);
}

