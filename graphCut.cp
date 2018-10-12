#include <cv.h>
#include <highgui.h>
#include <time.h>
#include <queue>
#include <algorithm>

#define UP      0
#define RIGHT   1
#define DOWN    2
#define LEFT    3

#define ALPHA   0.0
#define BETA    1.0

using namespace std;

int Vel[4][2] = { {-1,0},{0,1},{1,0},{0,-1} };

typedef struct stPixel
{
    int val;
    int prev;
    int visited;
    double f_alpha;
    double f_beta;
    double t_alpha;
    double t_beta;
    double c[4];
    double f[4];
}Pixel;

int nRows, nColumns, nRound;
int colA[3], colB[3];
int **I, **Inoise;
double **Ca, **Fa;
double **Cb, **Fb;
Pixel **Pix;

/************************************************
 *              Function Definitions            *
 ************************************************/

double Ford_Fulkerson();
void Cut_Graph(int, int);
void Clean_Pixels();
double V_Val(int, int, int);
double D_Val(int, int, double);
int ** ConvImgToInt(IplImage *);
void ConvIntToImg(int **, IplImage *);
int ** Create_Std_Image();
int ** Create_Random_Image();
int ** Add_Noise(int **, int);
void Print_Image(int **);

/************************************************
 *                  Main Program                *
 ************************************************/

int main()
{
    int i, j, k;
    double minCut;
    IplImage *imgSrc;
    IplImage *imgDst;
    
    //FILE *fout = fopen("out.txt", "w");
    
    srand(time(NULL));
    
    /************************************************
     * Create a default 2D-array and add noise to
     * it.
     ************************************************/
    
    I = Create_Random_Image();
    Inoise = Create_Random_Image();
    Print_Image(Inoise);
    
    /*
    cvNamedWindow("Corrupted Image");
    cvNamedWindow("Clean Image");
    
    imgSrc = cvLoadImage("landscapeNoise.png");
    imgDst = cvCloneImage(imgSrc);

    nRows = imgSrc->height;
    nColumns = imgSrc->width;
    
    I = ConvImgToInt(imgSrc);
    Inoise = ConvImgToInt(imgSrc);
    
    printf("%d %d %d\n", colA[0], colA[1], colA[2]);
    printf("%d %d %d\n", colB[0], colB[1], colB[2]);
     */
    
    /************************************************
     * Print the corrupted array
     ************************************************/
    
    printf("size: %d x %d\n", nRows, nColumns);
    //cvShowImage("Corrupted Image", imgSrc);
    
    /************************************************
     * Create the Pixel structure that will be used
     * in the Ford-Fulkerson algorithm
     ************************************************/
    
    Pix = new Pixel *[nRows];
    Ca = new double *[nRows];
    Fa = new double *[nRows];
    Cb = new double *[nRows];
    Fb = new double *[nRows];
    
    for(i=0; i<nRows; i++)
    {
        Pix[i] = new Pixel[nColumns];
        
        Ca[i] = new double[nColumns];
        Fa[i] = new double[nColumns];
        Cb[i] = new double[nColumns];
        Fb[i] = new double[nColumns];
        
        for(j=0; j<nColumns; j++)
        {
            Pix[i][j].val = Inoise[i][j];
            Pix[i][j].visited = 0;
        }
    }
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
        {
            for(k=0; k<4; k++)
            {
                Pix[i][j].f[k] = 0.0;
                Pix[i][j].c[k] = V_Val(i,j,k);
            }
            
            Pix[i][j].f_alpha = 0.0;
            Pix[i][j].f_beta = 0.0;
            Pix[i][j].t_alpha = D_Val(i,j, ALPHA);
            Pix[i][j].t_beta = D_Val(i,j, BETA);
            
            Fa[i][j] = 0.0;
            Fb[i][j] = 0.0;
            Ca[i][j] = Pix[i][j].t_alpha;
            Cb[i][j] = Pix[i][j].t_beta;
        }
    }
    
    printf("Running Ford-Fulkerson....\n");
    
    minCut = Ford_Fulkerson();
    
    printf("minCut = %.5lf\n", minCut);
    
    Clean_Pixels();
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
        {
            if(Pix[i][j].visited == 0 && (Ca[i][j] - Fa[i][j] > 0.0))
                Cut_Graph(i,j);
        }
    }
    
    printf("Cutting graph...\n");
    
    for(i=0; i<nRows; i++)
        for(j=0; j<nColumns; j++)
            I[i][j] = Pix[i][j].visited;
    
    Print_Image(I);
    
    /*
    ConvIntToImg(I, imgDst);
    
    cvShowImage("Clean Image", imgDst);
    cvSaveImage("landscapeClean.png", imgDst);
     */
    
    printf("Program Finished!! Press ESC to exit\n");
    
    while(1)
    {
        char c = cvWaitKey(33);
        if(c == 27)
            break;
    }
    
    //fclose(fout);
    
    /*
    cvDestroyWindow("Corrupted Image");
    cvDestroyWindow("Clean Image");
    
    cvReleaseImage(&imgSrc);
    cvReleaseImage(&imgDst);
     */
    
    return 0;
}

void Cut_Graph(int r, int c)
{
    int k, r2, c2;
    
    Pix[r][c].visited = 1;
    
    for(k=0; k<4; k++)
    {
        r2 = r + Vel[k][0];
        c2 = c + Vel[k][1];
        
        if(Pix[r][c].c[k] - Pix[r][c].f[k] > 0.0 && Pix[r2][c2].visited == 0)
            Cut_Graph(r2,c2);
    }
}

double Ford_Fulkerson()
{
    int i, j, k;
    int r, c;
    int r2, c2;
    int pathFounded;
    double minFlow, maxFlow = 0.0;
    queue<int> Q;
    
    do
    {
        pathFounded = 0;
        
        /************************************************
         * Find a pixel where the flow can pass from
         * ALPHA to that pixel.
         ************************************************/
    
        for(i=0; i<nRows; i++)
        {
            if(pathFounded == 1)    
                break;
            
            for(j=0; j<nColumns; j++)
            {
                if(pathFounded == 1)    
                    break;
                
                /************************************************
                 * Once you find it, add it to the queue and
                 * start the BFS
                 ************************************************/
                
                if(Ca[i][j] - Fa[i][j] > 0.0)
                {
                    while(!Q.empty())
                        Q.pop();
                    
                    nRound++;
                    
                    Pix[i][j].visited = nRound;
                    Pix[i][j].prev = -1;
                    Q.push(i*nColumns + j);
                    
                    while(!Q.empty() && pathFounded == 0)
                    {
                        k = Q.front();
                        Q.pop();
                        
                        r = k/nColumns;
                        c = k%nColumns;
                        
                        if(Pix[r][c].t_beta - Pix[r][c].f_beta > 0.0)  //The flow can pass trough BETA?
                        {
                            minFlow = (Pix[r][c].t_beta - Pix[r][c].f_beta);
                            
                            r2 = r;
                            c2 = c;
                            
                            while(Pix[r2][c2].prev != -1)
                            {
                                k = Pix[r2][c2].prev;
                                r2 = r2 + Vel[k][0];
                                c2 = c2 + Vel[k][1];
                                
                                minFlow = min(minFlow, Pix[r2][c2].c[(k + 2)%4] - Pix[r2][c2].f[(k+2)%4]);
                            }
                            
                            minFlow = min(minFlow, Ca[r2][c2] - Fa[r2][c2]);
                            
                            maxFlow = maxFlow + minFlow;
                            
                            Pix[r][c].f_beta = Pix[r][c].f_beta  + minFlow;
                            Fb[r][c] = Fb[r][c] - minFlow;
                            
                            r2 = r;
                            c2 = c;
                            
                            while(Pix[r2][c2].prev != -1)
                            {
                                k = Pix[r2][c2].prev;
                                
                                Pix[r2][c2].f[k] = Pix[r2][c2].f[k] - minFlow;
                                
                                r2 = r2 + Vel[k][0];
                                c2 = c2 + Vel[k][1];
                                
                                Pix[r2][c2].f[(k+2)%4] = Pix[r2][c2].f[(k+2)%4] + minFlow;
                            }
                            
                            Pix[r2][c2].f_alpha = Pix[r2][c2].f_alpha - minFlow;
                            Fa[r2][c2] = Fa[r2][c2] + minFlow;
                            
                            pathFounded = 1;
                            
                        }
                        else
                        {
                            for(k=0; k<4; k++)
                            {
                                if(Pix[r][c].c[k] - Pix[r][c].f[k] > 0.0)
                                {
                                    r2 = r + Vel[k][0];
                                    c2 = c + Vel[k][1];
                                    
                                    if(Pix[r2][c2].visited != nRound)
                                    {
                                        Pix[r2][c2].visited = nRound;
                                        Pix[r2][c2].prev = (k+2)%4;
                                        Q.push(r2*nColumns + c2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    while(pathFounded == 1);
    
    return maxFlow;
}

void Clean_Pixels()
{
    int i, j;
    
    for(i=0; i<nRows; i++)
        for(j=0; j<nColumns; j++)
            Pix[i][j].visited = 0;
}

double V_Val(int r, int c, int dir)
{
    double val;
    
    if(r + Vel[dir][0] < 0 || c + Vel[dir][1] < 0 || r + Vel[dir][0] >= nRows || c + Vel[dir][1] >= nColumns)
        val = 0.0;
    else
        val = fabs(double(Pix[r][c].val) - double(Pix[r+Vel[dir][0]][c+Vel[dir][1]].val));
    
    return val;
}

double D_Val(int r, int c, double D)
{
    int k, r2, c2;
    
    double val = (D - double(Pix[r][c].val))*(D - double(Pix[r][c].val));
    
    for(k=0; k<4; k++)
    {
        r2 = r + Vel[k][0];
        c2 = c + Vel[k][1];
        
        if(r2 >= 0 && r2 < nRows && c2 >= 0 && c2 < nColumns)
            val = val + (D - double(Pix[r2][c2].val))*(D - double(Pix[r2][c2].val));
    }
    
    return val;
}

int ** ConvImgToInt(IplImage *img)
{
    int **T;
    int i, j, aux = 0;
    CvScalar col;

    colA[0] = (int)cvGet2D(img,0,0).val[0];
    colA[1] = (int)cvGet2D(img,0,0).val[1];
    colA[2] = (int)cvGet2D(img,0,0).val[2];
    
    T = new int *[img->height];
    for(i=0; i<img->height; i++)
    {
        T[i] = new int[img->width];
        
        for(j=0; j<img->width; j++)
        {
            col = cvGet2D(img, i, j); 

            if((int)col.val[0] == colA[0] && (int)col.val[1] == colA[1] && (int)col.val[2] == colA[2])
                T[i][j] = 0;
            else
            {
                T[i][j] = 1;
                
                if(aux == 0)
                {
                    colB[0] = (int)col.val[0];
                    colB[1] = (int)col.val[1];
                    colB[2] = (int)col.val[2];
                    
                    aux = 1;
                }
            }
        }
    }
    
    return T;
}

void ConvIntToImg(int **T, IplImage *img)
{
    int i, j;
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
        {
            if(T[i][j] == 0)
                cvSet2D(img, i, j, cvScalar(colA[0], colA[1], colA[2]));
            else
                cvSet2D(img, i, j, cvScalar(colB[0], colB[1], colB[2]));
        }
    }
}

int ** Create_Std_Image()
{
    int i, j;
    int **T;
    
    nRows = 20;
    nColumns = 20;
    
    T = new int *[nRows];
    
    for(i=0; i<nRows; i++)
    {
        T[i] = new int[nColumns];
        
        for(j=0; j<nColumns; j++)
        {
            if(j >= nColumns/2)
                T[i][j] = 1;
            else
                T[i][j] = 0;
        }
    }
    
    return T;
}

int ** Create_Random_Image()
{
    int i, j;
    int **T;
    
    nRows = 20;
    nColumns = 20;
    
    T = new int *[nRows];
    
    for(i=0; i<nRows; i++)
    {
        T[i] = new int[nColumns];
        
        for(j=0; j<nColumns; j++)
            T[i][j] = rand()%2;
    }
    
    return T;
}

int ** Add_Noise(int **img, int p)
{
    int i, j, r;
    int **T;
    
    nRows = 20;
    nColumns = 20;
    
    T = new int *[nRows];
    
    for(i=0; i<nRows; i++)
    {
        T[i] = new int[nColumns];
        
        for(j=0; j<nColumns; j++)
        {
            r = rand()%10 + 1;
            if(r <= p)
                T[i][j] = (img[i][j] + 1)%2;
            else
                T[i][j] = img[i][j];
        }
    }
    
    return T;
}

void Print_Image(int **img)
{
    int i, j;
    
    for(j=0; j<nColumns; j++)
        printf("%s", "---");
    printf("-\n");
    
    for(i=0; i<nRows; i++)
    {
        printf("|");
        
        for(j=0; j<nColumns; j++)
        {
            printf("%2d", img[i][j]);
            if(j < nColumns-1 && img[i][j+1] != img[i][j])
                printf("|");
            else
                printf(" ");
        }
        
        printf("|\n|");

        for(j=0; j<nColumns; j++)
        {
            if(i < nRows - 1 && img[i][j] != img[i+1][j])
                printf("%s", "---");
            else
            {
                if(i == nRows-1)
                    printf("%s", "---");
                else
                    printf("%s", "   ");
            }
        }
        
        printf("|\n");
    }
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
            printf("%3d", img[i][j]);
        printf("\n");
    }
}