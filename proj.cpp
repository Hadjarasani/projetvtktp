#include "exo-vtk-include.h"
#include "config.h"
#include "helpers.h"
#include <mpi.h>


/*  //Enlever le premier slash pour commenter
#define FICHIER  "Frog_CHAR_X_256_Y_256_Z_44.raw"

int gridSize = 256;
int YgridSize = 256;
int ZgridSize = 44;

#define CHAR
#define SMALL

int startexploreval=13;
int endexploreval=20;//*/


/* //Ajouter un slash pour décommenter
 #define FICHIER  "Mystere1_SHORT_X_512_Y_512_Z_134.raw"
 int gridSize = 512;
 int YgridSize = 512;
 int ZgridSize = 134;
 
 #define SHORT
 #define SMALL
 
 int startexploreval=20000;
 int endexploreval=25000;//*/


/*
#define FICHIER  "Mystere2_SHORT_X_512_Y_499_Z_512.raw"

 
 int gridSize = 512;
 int YgridSize = 499;
 int ZgridSize = 512;

 #define SHORT
 #define SMALL

 int startexploreval=25000;
 int endexploreval=65000;//*/


//*
 #define FICHIER  "Mystere5_SHORT_X_2048_Y_2048_Z_756.raw"
 
 int gridSize = 2048;
 int YgridSize = 2048;
 int ZgridSize = 756;
 
 #define SHORT
 #define BIG
 
 int startexploreval=20000;
 int endexploreval=20050;//*/

/*
 #define FICHIER  "Mystere6_CHAR_X_1118_Y_2046_Z_694.raw"
 
 int gridSize = 1118;
 int YgridSize = 2046;
 int ZgridSize = 694;
 
 #define CHAR
 #define BIG
 
 int startexploreval=13;   //var posible 12 9
 int endexploreval=16;//*/  //var posible 255


/*
 #define FICHIER  "Mystere4_SHORT_X_512_Y_512_Z_322.raw"
 int gridSize = 512;
 int YgridSize = 512;
 int ZgridSize = 322;
 
 #define SHORT
 #define SMALL
 
 int startexploreval=50000;
 int endexploreval=65000; //*/

/*
 #define FICHIER "Mystere10_CHAR_X_1204_Y_1296_Z_224.raw"
 
 int gridSize = 1204;
 int YgridSize = 1296;
 int ZgridSize = 224;
 
 #define CHAR
 #define BIG
 
 int startexploreval=50;
 int endexploreval=255;//*/

/*
#define FICHIER "Mystere8_CHAR_X_2048_Y_2048_Z_2048.raw"
int gridSize = 2048;
 int YgridSize = 2048;
 int ZgridSize = 2048;
 
 #define CHAR
 #define BIG
 
 int startexploreval=100;
 int endexploreval=255;//*/

/*
#define FICHIER "Mystere11_SHORT_X_512_Y_512_Z_1024.raw"

#define SHORT
#define SMALL

int gridSize = 512;
 int YgridSize = 512;
 int ZgridSize = 1024;

 int startexploreval=47000;
 int endexploreval=65000;// */


/*
#define FICHIER  "Mystere9_SHORT_X_2048_Y_2048_Z_1444.raw" 
int gridSize = 2048;
int YgridSize = 2048;
int ZgridSize = 1444;

#define SHORT
#define SMALL


int startexploreval=30000;
int endexploreval=350000; //*/


int winSize = 1400;

int numPasses = 10;
int nbimages = 10;


const char *prefix = "";



int passNum = 0;

using std::cerr;
using std::endl;

vtkRectilinearGrid *ReadGrid(int zStart, int zEnd);
void WriteImage(const char *name, const float *rgba, int width, int height);
bool ComposeImageZbuffer(float *rgba_out, float *zbuffer,   int image_width, int image_height);


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int procRank, procSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    
    srand(time(NULL));
    GetMemorySize("initialization"); 
    int theStart = timer->StartTimer();

    int totalPixels = winSize * winSize;
    vtkRectilinearGrid *gridData = NULL;
    vtkLookupTable *lut = vtkLookupTable::New();
    lut->SetHueRange(0.1, 0.0);
    lut->SetSaturationRange(0.0, 1.0);
    lut->SetValueRange(1.0, 255.0);
    lut->SetNumberOfColors(100);
    lut->Build();
    
    vtkRenderer *ren = vtkRenderer::New();
    double bounds[6] = {0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001};
    ren->ResetCamera(bounds);
    
    vtkRenderWindow *renderWindow = vtkRenderWindow::New();
    renderWindow->SetSize(winSize, winSize);
    renderWindow->AddRenderer(ren);
    
    int zBlockSize = ZgridSize / procSize;
    int zStart = procRank * zBlockSize;
    int zEnd = zStart + zBlockSize;
    if (procRank == procSize - 1) {
        zEnd = ZgridSize;
    }
    
    int numPasses = 32;
    int sliceSize = (zEnd - zStart) / numPasses;
    int localStart = zStart;
    int localEnd = localStart + sliceSize;
    
    gridData = ReadGrid(localStart, localEnd);

    
    
    vtkContourFilter *contourFilter = vtkContourFilter::New();
    contourFilter->SetNumberOfContours(1);
    int contourValue = startexploreval;
    contourFilter->SetValue(0, contourValue);
    contourFilter->SetInputData(gridData);
    
    int maxsize = std::max(gridSize, std::max(YgridSize, ZgridSize));
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Scale(gridSize / (float)maxsize, YgridSize / (float)maxsize, ZgridSize / (float)maxsize);
    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputConnection(contourFilter->GetOutputPort());
    transformFilter->SetTransform(transform);
    
    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputConnection(transformFilter->GetOutputPort());
    
    vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);
    mapper->SetScalarRange(startexploreval, endexploreval);
    mapper->SetLookupTable(lut);

    ren->AddActor(actor);
    ren->SetViewport(0, 0, 1, 1);

    renderWindow->SetOffScreenRendering(1);
    vtkCamera *camera = ren->GetActiveCamera();
    renderWindow->Render();
    
    for (int al = 0; al < 1; al++) {
        //camera->Azimuth(150); //variable a modifier selon l'image
        camera->Elevation(195); // variable a modifier selon l'image
       

        renderWindow->Render();
        float *rgba = new float[totalPixels * 4];
        float *zBuffer = new float[totalPixels];
        for (int p = 0; p < totalPixels; p++) {
            rgba[4 * p] = 0;
            rgba[4 * p + 1] = 0;
            rgba[4 * p + 2] = 0;
            rgba[4 * p + 3] = 0;
            zBuffer[p] = 1.0;
        }
        for (int pass = 0; pass < numPasses; pass++) {
            localStart = zStart + pass * sliceSize;
            localEnd = localStart + sliceSize;
            if (pass == numPasses - 1) {
                localEnd = zEnd;
            }
            gridData->Delete();
            gridData = ReadGrid(localStart, localEnd);

            contourFilter->SetInputData(gridData);
            mapper->Update();
            renderWindow->Render();
            float *secondrgba = renderWindow->GetRGBAPixelData(0, 0, winSize - 1, winSize - 1, 0);
            float *secondzBuffer = renderWindow->GetZbufferData(0, 0, winSize - 1, winSize - 1);
            for (int i = 0; i < totalPixels; i++) {
                if (zBuffer[i] > secondzBuffer[i]) {
                    rgba[4 * i] = secondrgba[4 * i];
                    rgba[4 * i + 1] = secondrgba[4 * i + 1];
                    rgba[4 * i + 2] = secondrgba[4 * i + 2];
                    rgba[4 * i + 3] = secondrgba[4 * i + 3];
                    zBuffer[i] = secondzBuffer[i];
                }
            }
            delete[] secondrgba;
            delete[] secondzBuffer;
        }
        /*
        std::string name = std::to_string(procrank)+"_al"+std::to_string(al)+".png";
        WriteImage(name.c_str(), rgba, winSize, winSize);//*/
        float *Rgba3 = nullptr;
        float *ZBuffer3 = nullptr;
        if (procRank == 0) {
            Rgba3 = new float[4 * totalPixels * procSize];
            ZBuffer3 = new float[totalPixels * procSize];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(rgba, 4 * totalPixels, MPI_FLOAT, Rgba3, 4 * totalPixels, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Gather(zBuffer, totalPixels, MPI_FLOAT, ZBuffer3, totalPixels, MPI_FLOAT, 0, MPI_COMM_WORLD);
        delete[] rgba;
        delete[] zBuffer;

        if (procRank == 0) {
            std::vector<float> finalRgba1(totalPixels * 4, 0.0f);
            std::vector<float> finalZBuffer1(totalPixels, 1.0f);

             float* finalRgba = finalRgba1.data();
            float* finalZBuffer = finalZBuffer1.data();

            for (int i = 0; i < procSize; i++) {
                for (int index = 0; index < totalPixels; index++) {
                    int index2 = i * totalPixels + index;
                    if (ZBuffer3[index2] < finalZBuffer[index]) {
                        finalZBuffer[index] = ZBuffer3[index2];
                        finalRgba[index * 4] = Rgba3[index2 * 4];
                        finalRgba[index * 4 + 1] = Rgba3[index2 * 4 + 1];
                        finalRgba[index * 4 + 2] = Rgba3[index2 * 4 + 2];
                        finalRgba[index * 4 + 3] = Rgba3[index2 * 4 + 3];
                    }
                }
            }

             std::string finalfile = "final" + std::to_string(al) + ".png";
            WriteImage(finalfile.c_str(), finalRgba, winSize, winSize);
            GetMemorySize("end");
            timer->StopTimer(theStart, "Time elapsed");
            
            
        }
    }

    gridData->Delete();
    mapper->Delete();
    contourFilter->Delete();
    ren->RemoveActor(actor);
    actor->Delete();
    ren->Delete();
    renderWindow->Delete();
    MPI_Finalize();        
}


// You should not need to modify these routines.

vtkRectilinearGrid *
ReadGrid(int zStart, int zEnd)
{
    int  i;
    std::string file=(MY_DATA_PATH+ (std::string )FICHIER);
    
    ifstream ifile(file.c_str());
    if (ifile.fail())
    {
        cerr << prefix << "Unable to open file: " << MY_DATA_PATH+ (std::string )FICHIER<< "!!" << endl;
        throw std::runtime_error("can't find the file!! Check the name and the path of this file? ");
    }
    
    cerr << prefix << "Reading from " << zStart << " to " << zEnd << endl;
    
    vtkRectilinearGrid *rg = vtkRectilinearGrid::New();
    vtkFloatArray *X = vtkFloatArray::New();
    X->SetNumberOfTuples(gridSize);
    for (i = 0 ; i < gridSize ; i++)
        X->SetTuple1(i, i/(gridSize-1.0));
    rg->SetXCoordinates(X);
    X->Delete();
    vtkFloatArray *Y = vtkFloatArray::New();
    Y->SetNumberOfTuples(YgridSize);
    for (i = 0 ; i < YgridSize ; i++)
        Y->SetTuple1(i, i/(YgridSize-1.0));
    rg->SetYCoordinates(Y);
    Y->Delete();
    vtkFloatArray *Z = vtkFloatArray::New();
    int numSlicesToRead = zEnd-zStart+1;
    Z->SetNumberOfTuples(numSlicesToRead);
    for (i = zStart ; i <= zEnd ; i++)
        Z->SetTuple1(i-zStart, i/(ZgridSize-1.0));
    rg->SetZCoordinates(Z);
    Z->Delete();
    
    rg->SetDimensions(gridSize, YgridSize, numSlicesToRead);
    
    unsigned int valuesPerSlice  = gridSize*YgridSize;
    
#if defined(SHORT)
    unsigned int bytesPerSlice   = sizeof(short)*valuesPerSlice;
    
#elif defined(CHAR)
    unsigned int bytesPerSlice   = sizeof(char)*valuesPerSlice;
    
#elif  defined(FLOAT)
    unsigned int bytesPerSlice   = sizeof(float)*valuesPerSlice;
    
#else
#error Unsupported choice setting
#endif
    
    
#if defined(SMALL)
    unsigned int offset          = (unsigned int)zStart * (unsigned int)bytesPerSlice;
    unsigned int bytesToRead     = bytesPerSlice*numSlicesToRead;
    unsigned int valuesToRead    = valuesPerSlice*numSlicesToRead;
#elif defined(BIG)
    unsigned long long offset          = (unsigned long long)zStart * bytesPerSlice;
    unsigned long long  bytesToRead     = (unsigned long long )bytesPerSlice*numSlicesToRead;
    unsigned int valuesToRead    = (unsigned int )valuesPerSlice*numSlicesToRead;
#else
#error Unsupported choice setting
#endif
    
    
    
#if defined(SHORT)
    vtkUnsignedShortArray *scalars = vtkUnsignedShortArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    unsigned short *arr = scalars->GetPointer(0);
    
#elif defined(CHAR)
    vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    unsigned char *arr = scalars->GetPointer(0);
    
#elif  defined(FLOAT)
    vtkFloatArray *scalars = vtkFloatArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    float *arr = scalars->GetPointer(0);
#else
#error Unsupported choice setting
#endif
    
    
    
    
    
    ifile.seekg(offset, ios::beg);
    ifile.read((char *)arr, bytesToRead);
    ifile.close();
    
    int min=+2147483647;
    int max =0;
    
#if defined(SMALL)
    for (unsigned int i = 0 ; i < valuesToRead ; i++){
#elif defined(BIG)
        for (unsigned long long int i = 0 ; i < valuesToRead ; i++){
#endif
            
            if (min>(scalars->GetPointer(0))[i]) min=(scalars->GetPointer(0))[i];
            if (max<(scalars->GetPointer(0))[i]) max=(scalars->GetPointer(0))[i];
            
            if(rand()%(valuesToRead/20)==0){
#if defined(SHORT)
                std::cout<<(unsigned short)(scalars->GetPointer(0))[i]<<" ";
#elif defined(CHAR)
                std::cout<<(unsigned short)(scalars->GetPointer(0))[i]<<" ";
#elif  defined(FLOAT)
                std::cout<<(float)(scalars->GetPointer(0))[i]<<" ";
#else
#error Unsupported choice setting
#endif
                
            }
        }
        
        
        
        std::cout<<"min value read: "<<min<<endl;
        std::cout<<"max value read: "<<max<<endl;
        std::fflush(stdout);
        
        
        scalars->SetName("entropy");
        rg->GetPointData()->AddArray(scalars);
        scalars->Delete();
        
        vtkFloatArray *pr = vtkFloatArray::New();
        pr->SetNumberOfTuples(valuesToRead);
#if defined(SMALL)
        for (unsigned int i = 0 ; i < valuesToRead ; i++)
#elif defined(BIG)
            for (unsigned long long  i = 0 ; i < valuesToRead ; i++)
#endif
                pr->SetTuple1(i, passNum);
        
        pr->SetName("pass_num");
        rg->GetPointData()->AddArray(pr);
        pr->Delete();
        
        rg->GetPointData()->SetActiveScalars("entropy");
        
        cerr << file << " Done reading" << endl;
        return rg;
    }
    
    
    
    void
    WriteImage(const char *name, const float *rgba, int width, int height)
    {
        vtkImageData *img = vtkImageData::New();
        img->SetDimensions(width, height, 1);
#if VTK_MAJOR_VERSION <= 5
        img->SetNumberOfScalarComponents(3);
        img->SetScalarTypeToUnsignedChar();
#else
        img->AllocateScalars(VTK_UNSIGNED_CHAR,3);
#endif
        
        for (int i = 0 ; i < width ; i++)
            for (int j = 0 ; j < height ; j++)
            {
                unsigned char *ptr = (unsigned char *) img->GetScalarPointer(i, j, 0);
                int idx = j*width + i;
                ptr[0] = (unsigned char) (255*rgba[4*idx + 0]);
                ptr[1] = (unsigned char) (255*rgba[4*idx + 1]);
                ptr[2] = (unsigned char) (255*rgba[4*idx + 2]);
            }
        
        
        vtkPNGWriter *writer = vtkPNGWriter::New();
        writer->SetInputData(img);
        writer->SetFileName(name);
        writer->Write();
        
        img->Delete();
        writer->Delete();
    }
    
    
    bool ComposeImageZbuffer(float *rgba_out, float *zbuffer,   int image_width, int image_height)
    {
        int npixels = image_width*image_height;
        
        float min=1;
        float max=0;
        for (int i = 0 ; i < npixels ; i++){
            if (zbuffer[i]<min) min=zbuffer[i];
            if (zbuffer[i]>max) max=zbuffer[i];
            
        }
        std::cout<<"min:"<<min;
        std::cout<<"max:"<<max<<"  ";
        
        float coef=1.0/((max-min));
        
        std::cout<<"coef:"<<coef<<"  ";
        
        for (int i = 0 ; i < npixels ; i++){
            
            rgba_out[i*4] = (zbuffer[i]==1.0?0:1-coef*(zbuffer[i]-min));
            rgba_out[i*4+1] = (zbuffer[i]==1.0?0:1-coef*(zbuffer[i]-min));
            rgba_out[i*4+2] = (zbuffer[i]==1.0?0:1-coef*(zbuffer[i]-min));
            rgba_out[i*4+3] = 0.0;
        }
        
        
        return false;
    }