//
//  convolution.c
//
//
//  Created by Josep Lluis Lerida on 11/03/15.
//
// This program calculates the convolution for PPM images.
// The program accepts an PPM image file, a text definition of the kernel matrix and the PPM file for storing the convolution results.
// The program allows to define image partitions for processing large images (>500MB)
// The 2D image is represented by 1D vector for chanel R, G and B. The convolution is applied to each chanel separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

// Estructura per emmagatzemar el contingut d'una imatge.
struct imagenppm{
    int altura;
    int ancho;
    char *comentario;
    int maxcolor;
    int P;
    int *R;
    int *G;
    int *B;
    long blockSize;
};
typedef struct imagenppm* ImagenData;

// Estructura per emmagatzemar el contingut d'un kernel.
struct structkernel{
    int kernelX;
    int kernelY;
    float *vkern;
};
typedef struct structkernel* kernelData;

//Functions Definition
ImagenData initimage(char* nombre, FILE **fp, int partitions, int halo);
ImagenData duplicateImageData(ImagenData src, int partitions, int halo);

int readImage(ImagenData Img, FILE **fp, int dim,int halosize, long int *position);
int duplicateImageChunk(ImagenData src, ImagenData dst, int dim);
int initfilestore(ImagenData img, FILE **fp, char* nombre, long *position);
int savingChunk(ImagenData img, FILE **fp, int dim, int offset);
int convolve2D(int* inbuf, int* outbuf, int sizeX, int sizeY, float* kernel, int ksizeX, int ksizeY);
void freeImagestructure(ImagenData *src);

//Open Image file and image struct initialization
ImagenData initimage(char* nombre, FILE **fp,int partitions, int halo){
    char c;
    char comentario[300];
    int i=0,chunk=0;
    ImagenData img=NULL;
    
    /*Se habre el fichero ppm*/

    if ((*fp=fopen(nombre,"r"))==NULL){
        perror("Error: ");
    }
    else{
        //Memory allocation
        img=(ImagenData) malloc(sizeof(struct imagenppm));
        //Reading the first line: Magical Number "P3"
        fscanf(*fp,"%c%d ",&c,&(img->P));
        
        //Reading the image comment
        while((c=fgetc(*fp))!= '\n'){comentario[i]=c;i++;}
        comentario[i]='\0';
        //Allocating information for the image comment
        img->comentario = (char*) malloc(300*sizeof(char));
        strcpy(img->comentario,comentario);
        //Reading image dimensions and color resolution
        fscanf(*fp,"%d %d %d",&img->ancho,&img->altura,&img->maxcolor);
        chunk = img->ancho*img->altura / partitions;
        //We need to read an extra row.

        chunk = chunk + img->ancho * halo;
        img->blockSize = chunk;
        if ((img->R=(int*)malloc(chunk*sizeof(int))) == NULL) {return NULL;}
        if ((img->G=(int*)malloc(chunk*sizeof(int))) == NULL) {return NULL;}
        if ((img->B=(int*)malloc(chunk*sizeof(int))) == NULL) {return NULL;}
    }
    return img;
}

//Duplicate the Image struct for the resulting image
ImagenData duplicateImageData(ImagenData src, int partitions, int halo){
    char c;
    char comentario[300];
    unsigned int imageX, imageY;
    int i=0, chunk=0;
    //Struct memory allocation
    ImagenData dst=(ImagenData) malloc(sizeof(struct imagenppm));

    //Copying the magic number
    dst->P=src->P;
    //Copying the string comment
    dst->comentario = (char*)malloc(300*sizeof(char));
    strcpy(dst->comentario,src->comentario);
    //Copying image dimensions and color resolution
    dst->ancho=src->ancho;
    dst->altura=src->altura;
    dst->maxcolor=src->maxcolor;
    chunk = dst->ancho*dst->altura / partitions;
    //We need to read an extra row.
    chunk = chunk + src->ancho * halo;
    if ((dst->R=(int*)malloc(chunk*sizeof(int))) == NULL) {return NULL;}
    if ((dst->G=(int*)malloc(chunk*sizeof(int))) == NULL) {return NULL;}
    if ((dst->B=(int*)malloc(chunk*sizeof(int))) == NULL) {return NULL;}
    return dst;
}

//Read the corresponding chunk from the source Image
int readImage(ImagenData img, FILE **fp, int dim, int halosize, long *position){
    int i=0,haloposition=0;
    if (fseek(*fp,*position,SEEK_SET))
        perror("Error: ");
    haloposition = dim-(img->ancho*halosize*2);
    for(i=0;i<dim;i++) {
        // When start reading the halo store the position in the image file
        if (halosize != 0 && i == haloposition) *position=ftell(*fp);
        fscanf(*fp,"%d %d %d ",&img->R[i],&img->G[i],&img->B[i]);
        
    }
//    printf ("Readed = %d pixels, posicio=%lu\n",k,*position);
    return 0;
}
/*
//Read the corresponding chunk from the source Image
int readImage(ImagenData img, FILE **fp, int dim, long *position){
    int i=0;
    if (fseek(*fp,*position,SEEK_SET))
        perror("Error: ");
    for(i=0;i<dim;i++) {
        // When start reading the halo store the position in the image file
        fscanf(*fp,"%d %d %d ",&img->R[i],&img->G[i],&img->B[i]);      
    }
//    printf ("Readed = %d pixels, posicio=%lu\n",k,*position);
    return 0;
}*/

//Duplication of the  just readed source chunk to the destiny image struct chunk
int duplicateImageChunk(ImagenData src, ImagenData dst, int dim){
    int i=0;
    
    for(i=0;i<dim;i++){
        dst->R[i] = src->R[i];
        dst->G[i] = src->G[i];
        dst->B[i] = src->B[i];
    }
//    printf ("Duplicated = %d pixels\n",i);
    return 0;
}

// Open kernel file and reading kernel matrix. The kernel matrix 2D is stored in 1D format.
kernelData leerKernel(char* nombre){
    FILE *fp;
    int i=0;
    kernelData kern=NULL;
    
    /*Opening the kernel file*/
    fp=fopen(nombre,"r");
    if(!fp){
        perror("Error: ");
    }
    else{
        //Memory allocation
        kern=(kernelData) malloc(sizeof(struct structkernel));
        
        //Reading kernel matrix dimensions
        fscanf(fp,"%d,%d,", &kern->kernelX, &kern->kernelY);
        kern->vkern = (float *)malloc(kern->kernelX*kern->kernelY*sizeof(float));
        
        // Reading kernel matrix values
        for (i=0;i<(kern->kernelX*kern->kernelY)-1;i++){
            fscanf(fp,"%f,",&kern->vkern[i]);
        }
        fscanf(fp,"%f",&kern->vkern[i]);
        fclose(fp);
    }
    return kern;
}

// Open the image file with the convolution results
int initfilestore(ImagenData img, FILE **fp, char* nombre, long *position){
    /*Se crea el fichero con la imagen resultante*/
    if ( (*fp=fopen(nombre,"w")) == NULL ){
        perror("Error: ");
        return -1;
    }
    /*Writing Image Header*/
    fprintf(*fp,"P%d\n%s\n%d %d\n%d\n",img->P,img->comentario,img->ancho,img->altura,img->maxcolor);
    *position = ftell(*fp);
    return 0;
}

// Writing the image partition to the resulting file. dim is the exact size to write. offset is the displacement for avoid halos.
int savingChunk(ImagenData img, FILE **fp, int dim, int offset){
    int i,k=0;
    //Writing image partition
    for(i=offset;i<dim+offset;i++){
        fprintf(*fp,"%d %d %d ",img->R[i],img->G[i],img->B[i]);
        if ((i+1)%6==0) fprintf(*fp,"\n");
        k++;
    }
    printf ("Writed = %d pixels, dim=%d, offset=%d\n",k,dim, offset);
    return 0;
}

// This function free the space allocated for the image structure.
void freeImagestructure(ImagenData *src){
    
    free((*src)->comentario);
    free((*src)->R);
    free((*src)->G);
    free((*src)->B);
    
    free(*src);
}

///////////////////////////////////////////////////////////////////////////////
// 2D convolution
// 2D data are usually stored in computer memory as contiguous 1D array.
// So, we are using 1D array for 2D data.
// 2D convolution assumes the kernel is center originated, which means, if
// kernel size 3 then, k[-1], k[0], k[1]. The middle of index is always 0.
// The following programming logics are somewhat complicated because of using
// pointer indexing in order to minimize the number of multiplications.
//
//
// signed integer (32bit) version:
///////////////////////////////////////////////////////////////////////////////
int convolve2D(int* in, int* out, int dataSizeX, int dataSizeY,
               float* kernel, int kernelSizeX, int kernelSizeY)
{
    int i, j, m, n;
    int *inPtr, *inPtr2, *outPtr;
    float *kPtr;
    int kCenterX, kCenterY;
    int rowMin, rowMax;                             // to check boundary of input array
    int colMin, colMax;                             //
    float sum;                                      // temp accumulation buffer
    
    // check validity of params
    if(!in || !out || !kernel) return -1;
    if(dataSizeX <= 0 || kernelSizeX <= 0) return -1;
    
    // find center position of kernel (half of kernel size)
    kCenterX = (int)kernelSizeX / 2;
    kCenterY = (int)kernelSizeY / 2;
    
    // init working  pointers
    inPtr = inPtr2 = &in[dataSizeX * kCenterY + kCenterX];  // note that  it is shifted (kCenterX, kCenterY),
    outPtr = out;
    kPtr = kernel;
    
    // start convolution
    for(i= 0; i < dataSizeY; ++i)                   // number of rows
    {
        // compute the range of convolution, the current row of kernel should be between these
        rowMax = i + kCenterY;
        rowMin = i - dataSizeY + kCenterY;
        
        for(j = 0; j < dataSizeX; ++j)              // number of columns
        {
            // compute the range of convolution, the current column of kernel should be between these
            colMax = j + kCenterX;
            colMin = j - dataSizeX + kCenterX;
            
            sum = 0;                                // set to 0 before accumulate
            
            // flip the kernel and traverse all the kernel values
            // multiply each kernel value with underlying input data
            for(m = 0; m < kernelSizeY; ++m)        // kernel rows
            {
                // check if the index is out of bound of input array
                if(m <= rowMax && m > rowMin)
                {
                    for(n = 0; n < kernelSizeX; ++n)
                    {
                        // check the boundary of array
                        if(n <= colMax && n > colMin)
                            sum += *(inPtr - n) * *kPtr;
                        
                        ++kPtr;                     // next kernel
                    }
                }
                else
                    kPtr += kernelSizeX;            // out of bound, move to next row of kernel
                
                inPtr -= dataSizeX;                 // move input data 1 raw up
            }
            
            // convert integer number
            if(sum >= 0) *outPtr = (int)(sum + 0.5f);
//            else *outPtr = (int)(sum - 0.5f)*(-1);
            // For using with image editors like GIMP or others...
            else *outPtr = (int)(sum - 0.5f);
            // For using with a text editor that read ppm images like libreoffice or others...
//            else *outPtr = 0;
            
            kPtr = kernel;                          // reset kernel to (0,0)
            inPtr = ++inPtr2;                       // next input
            ++outPtr;                               // next output
        }
    }
    
    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
//////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int size, rank;

    int i=0,j=0,k=0;

    int *rgb_array;
//    int headstored=0, imagestored=0, stored;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int blockSize_array[size];

    if(argc != 4)
    {
        printf("Usage: %s <image-file> <kernel-file> <result-file> \n", argv[0]);
        
        printf("\n\nError, Missing parameters:\n");
        printf("format: ./serialconvolution image_file kernel_file result_file\n");
        printf("- image_file : source image path (*.ppm)\n");
        printf("- kernel_file: kernel path (text file with 1D kernel matrix)\n");
        printf("- result_file: result image path (*.ppm)\n");
        return -1;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // READING IMAGE HEADERS, KERNEL Matrix, DUPLICATE IMAGE DATA, OPEN RESULTING IMAGE FILE
    //////////////////////////////////////////////////////////////////////////////////////////////////
    int imagesize, partitions, heigth_part, missing_parts, chunksize, halo, halosize;
    long position=0;
    double start, tstart=0, tend=0, tread=0, tcopy=0, tconv=0, tstore=0, treadk=0;
    struct timeval tim;
    FILE *fpsrc=NULL,*fpdst=NULL;
    ImagenData source=NULL, output=NULL;

    char *image_file, *result_file, *kernel_file;

    // Store number of partitions
    partitions = size;
    image_file = argv[1];
    result_file = argv[3];
    kernel_file = argv[2];
    ////////////////////////////////////////

    //Reading kernel matrix
    gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    tstart = start;
    kernelData kern=NULL;
    if ( (kern = leerKernel(kernel_file))==NULL) {
        //        free(source);
        //        free(output);
        return -1;
    }
    //The matrix kernel define the halo size to use with the image. The halo is zero when the image is not partitioned.
    if (partitions==1) halo=0;
    else halo = (kern->kernelY/2)*2;
    gettimeofday(&tim, NULL);
    treadk = treadk + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);

    ////////////////////////////////////////
    //Reading Image Header. Image properties: Magical number, comment, size and color resolution.
    gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    //Memory allocation based on number of partitions and halo size.
    if ( (source = initimage(image_file, &fpsrc, partitions, halo)) == NULL) {
        return -1;
    }

    gettimeofday(&tim, NULL);
    tread = tread + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
    //Duplicate the image struct.
    gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    if ( (output = duplicateImageData(source, partitions, halo)) == NULL) {
        return -1;
    }

    gettimeofday(&tim, NULL);
    tcopy = tcopy + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
    
    ////////////////////////////////////////
    //Initialize Image Storing file. Open the file and store the image header.

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // CHUNK READING
    //////////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Status status;
    int offset=0;
    imagesize = source->altura*source->ancho;
    heigth_part = source->altura/partitions;
    missing_parts= source->altura % partitions;
    rgb_array = (int *) malloc(source->blockSize * 3 * sizeof(int));
//    printf("%s ocupa %dx%d=%d pixels. Partitions=%d, halo=%d, partsize=%d pixels\n", argv[1], source->altura, source->ancho, imagesize, partitions, halo, partsize);
    // Last node make partitions and send to other nodes
    if (rank==size-1)
    {

        if (initfilestore(output, &fpdst, result_file, &position)!=0) {
            perror("Error: ");
            //        free(source);
            //        free(output);
            return -1;
        }
        //For all node read a part of image and send them
        for (int i = 0; i < size-1; ++i)
        {
            
            if (i==0) {
                halosize  = halo/2;
                chunksize = (heigth_part + halosize) * source->ancho ;
            }
            else if(i<size-1) {
                halosize  = halo;
                chunksize = (heigth_part + halosize) * source->ancho ;
            }

            if (readImage(source, &fpsrc, chunksize, halo/2, &position)) {return -1;}


            blockSize_array[i] = source->blockSize;

            //Serialize
            for (int j = 0; j < (source->blockSize*3); ++j)
            {
                if (j<source->blockSize)
                {
                    rgb_array[j] = source->R[j];
                }else if (j>=source->blockSize && j< (source->blockSize*2))
                {
                    rgb_array[j] = source->G[(j-source->blockSize)];
                }else{
                    rgb_array[j] = source->B[(j-source->blockSize*2)];
                }
            }
            
            MPI_Send(rgb_array, source->blockSize * 3 , MPI_INT, i, 1, MPI_COMM_WORLD);
        }

        // Last node do the convolve and of her part and the missing part
        halosize  = halo/2;
        chunksize = (heigth_part + halosize +missing_parts)* source->ancho;

        if (readImage(source, &fpsrc, chunksize, halo/2, &position)) {return -1;}

        blockSize_array[size-1] = source->blockSize;

        convolve2D(source->R, output->R, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(source->G, output->G, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(source->B, output->B, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
       

    }else{

        MPI_Recv(rgb_array, source->blockSize * 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
                MPI_COMM_WORLD, &status);

            if ( (source = initimage(image_file, &fpsrc, partitions, halo)) == NULL) {return -1;}
            
            for (int i = 0; i < (source->blockSize*3); ++i)
            {
                if (i < source->blockSize)
                {
                    source->R[i] = rgb_array[i];
                }else if (i>=source->blockSize && i < (source->blockSize*2))
                {
                    source->G[i - source->blockSize] = rgb_array[i];
                }else{
                    source->B[i - (source->blockSize * 2)] = rgb_array[i];
                
                }
            }
            
            if (rank==0) {
                halosize  = halo/2;
                chunksize = (heigth_part + halo) * source->ancho;
            }
            else if(rank<size-1) {
                halosize  = halo;
                chunksize = (heigth_part + halo) * source->ancho;
            }
            
            if ( duplicateImageChunk(source, output, chunksize) ) {return -1;}

            convolve2D(source->R, output->R, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
            convolve2D(source->G, output->G, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
            convolve2D(source->B, output->B, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
       

       for (int i = 0; i < (source->blockSize*3); ++i)
        {
            if (i<source->blockSize)
            {
                rgb_array[i] = output->R[i];
            }else if (i>=source->blockSize && i< (source->blockSize*2))
            {
                rgb_array[i] = output->G[(i-source->blockSize)];
            }else{
                rgb_array[i] = output->B[(i-source->blockSize*2)];
            }
        }

        MPI_Send(rgb_array, source->blockSize * 3 , MPI_INT, size-1, 1, 
            MPI_COMM_WORLD);
    }

        
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // CHUNK SAVING
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //Storing resulting image partition.

    //Last node recieve the blocks with order and write the final image
    if (rank==size-1)
    {
        ImagenData output_aux=NULL;

        for (int i = 0; i < size-1; ++i)
        {

            MPI_Recv(rgb_array, blockSize_array[i] * 3, MPI_INT, i, MPI_ANY_TAG, 
                MPI_COMM_WORLD, &status);

            //printf("He rebut del fill: %d\n",i );
            int offset_aux;

            if(i==0){
                offset_aux = 0;
            }else{
                offset_aux = (source->ancho*halo/2);
            }

            if ( (output_aux = initimage(image_file, &fpsrc, partitions, halo)) == NULL) {
                return -1;
            }

            for (int j = 0; j < (blockSize_array[i]*3); ++j)
            {
                if (j < blockSize_array[i])
                {
                    output_aux->R[j] = rgb_array[j];
                }else if (j>=blockSize_array[i] && j < (blockSize_array[i]*2))
                {
                    output_aux->G[j - blockSize_array[i]] = rgb_array[j];
                }else{
                    output_aux->B[j - (blockSize_array[i] * 2)] = rgb_array[j];
                }
            }
            printf("Node: %i offset_aux:%i\n", i,offset_aux);
            if (savingChunk(output_aux, &fpdst, heigth_part*source->ancho, offset_aux)) {
                perror("Error: ");
                //        free(source);
                //        free(output);
                return -1;
            }
            freeImagestructure(&output_aux);

            
        }
        offset = (source->ancho*halo/2);

        
        if (savingChunk(output, &fpdst, heigth_part*(source->ancho+missing_parts), offset)) {
                perror("Error: ");
                //        free(source);
                //        free(output);
                return -1;
            }
        
        //printf("He rebut del fill: %d\n",status.MPI_SOURCE );


    }


    
    
    freeImagestructure(&source);
    freeImagestructure(&output);
    
    gettimeofday(&tim, NULL);
    tend = tim.tv_sec+(tim.tv_usec/1000000.0);
    /*
    printf("Imatge: %s\n", argv[1]);
    printf("ISizeX : %d\n", source->ancho);
    printf("ISizeY : %d\n", source->altura);
    printf("kSizeX : %d\n", kern->kernelX);
    printf("kSizeY : %d\n", kern->kernelY);
    printf("%.6lf seconds elapsed for Reading image file.\n", tread);
    printf("%.6lf seconds elapsed for copying image structure.\n", tcopy);
    printf("%.6lf seconds elapsed for Reading kernel matrix.\n", treadk);
    printf("%.6lf seconds elapsed for make the convolution.\n", tconv);
    printf("%.6lf seconds elapsed for writing the resulting image.\n", tstore);
    printf("%.6lf seconds elapsed\n", tend-tstart);
    */

    free(kern->vkern);
    free(kern);
    fclose(fpsrc);
    fclose(fpdst);
    MPI_Finalize();

    //printf("END ==%i===\n",rank);

    return 0;
}
