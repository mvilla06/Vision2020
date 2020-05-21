/*
 * Title: Window Click Example
 * Class: Vision para Robot
 * Instructor: Dr. Jose Luis Gordillo (http://robvis.mty.itesm.mx/~gordillo/)
 * Code: Manlio Barajas (manlito@gmail.com)
 * Institution: Tec de Monterrey, Campus Monterrey
 * Date: January 28, 2013
 *
 * Description: Shows the most basic interaction with OpenCV HighGui window.
 *
 * This programs uses OpenCV http://www.opencv.org/
 */

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <string.h>
#include <vector>
#include <queue>
//define only for training
//#define TRAIN
#define OBJECTS_TO_TRAIN 4

#ifdef TRAIN
#define OBJECTS_TO_FIND 1
#define SAMPLES_PER_OBJECT 10
int samples = 0;
int next_sample = 0;
int x = 1;
#else
#define OBJECTS_TO_FIND 2
#endif

using namespace cv;
using namespace std;

// Here we will store points
Point p(0, 0), q(0, 0);
unsigned char thVal = 0;

void getRegions(Mat &image, long int * ordinary_moments ){
    Mat color(image.rows, image.cols, CV_8UC1, Scalar(0));
    
    
    memset(ordinary_moments, 0, OBJECTS_TO_FIND*sizeof(long int)*6);
    

    for(int k = 1; k<OBJECTS_TO_FIND+1; k++){
        //printf("\n");
        int x, y, timeout=0;
        do{
        x = rand()%image.cols;
        y = rand()%image.rows;
        timeout++;
        }while((image.at<uchar>(y, x) ==0 || color.at<uchar>(y, x)!=0 )&& timeout<10000);
        
        vector<Point> vecinos;
        Point p(0, 1);
        vecinos.push_back(p);
        p.x = -1;
        p.y = 0;
        vecinos.push_back(p);
        p.x = 0;
        p.y = -1;
        vecinos.push_back(p);
        p.x = 1;
        p.y = 0;
        vecinos.push_back(p);

        Point q(x, y);
        queue<Point> lista;
        lista.push(q);
        color.at<uchar>(q) = k*80;
        Point punto;
        while(lista.empty()==false){
            
            punto.x = lista.front().x;
            punto.y = lista.front().y;

            ordinary_moments[(k-1)*6+0]++;
            ordinary_moments[(k-1)*6+1] += punto.x;
            ordinary_moments[(k-1)*6+2] += punto.y;
            ordinary_moments[(k-1)*6+3] += punto.x*punto.x;
            ordinary_moments[(k-1)*6+4] += punto.x*punto.y;
            ordinary_moments[(k-1)*6+5] += punto.y*punto.y;
            lista.pop();
            for(int i = 0; i<4; i++){
                
                if(image.at<uchar>(punto+vecinos[i]) ==255 && color.at<uchar>(punto+vecinos[i])==0){
                    lista.push(punto+vecinos[i]);
                    color.at<uchar>(punto+vecinos[i]) = k*80;
                }
            }

        }
        
    }
    imshow("Color", color);
    
}


/* This is the callback that will only display mouse coordinates */
void mouseCoordinatesExampleCallback(int event, int x, int y, int flags, void *param);
void mouseCallback(int event, int x, int y, int flags, void* param){
    
    switch (event)
    {
    case CV_EVENT_LBUTTONDOWN:
        
       
        break;
    case CV_EVENT_MOUSEMOVE:
        break;
    case CV_EVENT_LBUTTONUP:
         thVal = x>>1;

        break;
    }
}
void histogram(const Mat &image, char t, unsigned char *threshold, char o, char colorSpace);
void destroyHistogram()
{
    destroyWindow("Blue");
    destroyWindow("Green");
    destroyWindow("Red");
    destroyWindow("GreyScale");
    destroyWindow("Y");
    destroyWindow("In-phase");
    destroyWindow("Quadrature");
    destroyWindow("Hue");
    destroyWindow("Saturation");
    destroyWindow("Value");
}
void bgr2gs(Mat &image);
void bgr2hsv(Mat &image);
void bgr2yiq(Mat &image);

void draw(Mat image, char t, char *o, unsigned char *threshold);

void selection(Mat image, unsigned char *threshold, Mat original, char * r);



int main(int argc, char *argv[])
{
    /* First, open camera device */

    VideoCapture camera;
    camera.open(0);

    /* Create images where captured and transformed frames are going to be stored */
    Mat currentImage;
    Mat originalImage;
    /* Create main OpenCV window to attach callbacks */
    namedWindow("Image");
    setMouseCallback("Image", mouseCoordinatesExampleCallback);

    // Flags: t -> Freeze camera, h -> Build histograms, b -> Convert grayscale, y-> Convert YIQ
    // v -> Convert to HSV
    char t = 0, h = 0, o = 0, b = 0, v = 0, y = 0, f = 0, r = 0;
    char colorSpace = 0; // 0->BGR, 1->BW, 2->YIQ, 3->HSV
    char x;
    camera >> currentImage;

    unsigned char *threshold = new unsigned char[currentImage.channels() * 2];

    while (true)
    {

        /* Obtain a new frame from camera */
        if (!t)
        {
            camera >> currentImage;
            camera >> originalImage;
        }

        if (currentImage.data)
        {

            //Convert Image only with live feed
            if (!t)
            {
                if (y)
                {
                    bgr2yiq(currentImage);
                    b = 0;
                    v = 0;
                }

                if (b)
                {
                    bgr2gs(currentImage);
                    v = 0;
                    y = 0;
                }

                if (v)
                {
                    bgr2hsv(currentImage);
                    y = 0;
                    b = 0;
                }
            }

            if (!o)
            {
                destroyWindow("ROI");
                destroyWindow("Selection");
            }
            /* Show image */
            draw(currentImage, t, &o, threshold);

            if (o && f)
            {
                selection(currentImage, threshold, originalImage, &r);
            }
            else
                destroyWindow("Selection");

            //Show histograms
            if (h)
                histogram(currentImage, t, threshold, o, colorSpace);
            else
            {
                
                        
                    destroyHistogram();   
                
                
            }

            imshow("Image", originalImage);
            imshow("Processed Image", currentImage);
            x = waitKey(3);
            switch (x)
            {
            case ' ': //Freeze image
                if (!t)
                    q = p;
                t = ~t;
                
                break;
            case 'x': //End program
                goto end;
                break;
            case 'h': //Show histograms
                h = ~h;
                break;
            case 'b': //Greyscale
                
                if(!b)
                colorSpace = 1;
                else
                colorSpace = 0;
                b = ~b;
                o = 0;
                y = 0;
                v = 0;
                destroyHistogram();
                destroyWindow("Threshold");
                destroyWindow("Binary");
                break;
            case 'y': //YIQ
                if(!y)
                colorSpace = 2;
                else
                colorSpace = 0;
                y = ~y;
                o = 0;
                b = 0;
                v = 0;
                destroyHistogram();
                break;
            case 'v': //HSV
                if(!v)
                colorSpace = 3;
                else
                colorSpace = 0;
                v = ~v;
                o = 0;
                y = 0;
                b = 0;
                destroyHistogram();
                break;
            case 'f': //Filter image
                f = ~f;
                break;
            case 'r': //Get regions
                if(r){
                    destroyWindow("Color");
                    destroyWindow("Centroid");
                }
                r = ~r;
                break;
            #ifdef TRAIN
            case 'n': //next_sample sample
                next_sample = 1;
                break;
            #endif
            }
            #ifdef TRAIN
            if (samples==OBJECTS_TO_TRAIN*SAMPLES_PER_OBJECT)
            goto end;
            #endif
            
        }
        else
        {
            cout << "No image data.. " << endl;
            return 1;
        }
    }
end:
    delete[] threshold;

    return 0;
}

void draw(Mat image, char t, char *o, unsigned char *threshold)
{
    Mat roi;
    
    if (!t)
    {
        if (p.x > 0 && p.y > 0)
            circle(image, p, 5, Scalar(255, 255, 255), 1, 8, 0);
    }
    else
    {
        if (p != q)
        { //Valid ROI selection
            Rect rec(p, q);
            roi = image(rec);
            imshow("ROI", roi);
            //rectangle(img, p, q, Scalar(0, 0, 255), 1, 8, 0);

            *o = 1;
            memset(threshold, 0xff, image.channels() * sizeof(char));
            memset(threshold + image.channels() * sizeof(char), 0x00, image.channels() * sizeof(char));
            //Calculate the maxima and minima of the ROI
            for (int r = 0; r < roi.rows; r++)
                for (int c = 0; c < roi.cols; c++)
                    for (int ch = 0; ch < roi.channels(); ++ch)
                    {
                        if (roi.at<Vec3b>(r, c)[ch] < threshold[ch])
                            threshold[ch] = roi.at<Vec3b>(r, c)[ch];

                        if (roi.at<Vec3b>(r, c)[ch] > threshold[ch + roi.channels()])
                            threshold[ch + roi.channels()] = roi.at<Vec3b>(r, c)[ch];
                    }
        }
    }
}

void selection(Mat image, unsigned char *threshold, Mat original, char * r)
{
    Mat selectionImg = image.clone();
    Mat filteredImage(image.rows, image.cols, CV_8UC1, Scalar(255));
    for (int r = 0; r < selectionImg.rows; r++)
    {
        for (int c = 0; c < selectionImg.cols; c++)
        {
            bool totalCond = true;
            for (int ch = 0; ch < selectionImg.channels(); ch++)
            {
                bool cond = selectionImg.at<Vec3b>(r, c)[ch] > threshold[ch] && selectionImg.at<Vec3b>(r, c)[ch] < threshold[ch + selectionImg.channels()];
                totalCond = totalCond && cond;
            }
            if (!totalCond)
            {
                
                filteredImage.at<uchar>(r, c) = 0;
                
            }
            else
            {
                
                    filteredImage.at<uchar>(r, c) = 0xFF;
            }
        }
    }

    morphologyEx(filteredImage, filteredImage, MORPH_OPEN, getStructuringElement(MORPH_ELLIPSE, Size (3,3)), Point(-1, -1) , 3); 
    imshow("Selection", filteredImage);
    if(*r){
    static FILE * file = NULL;
    #ifdef TRAIN
    if(!file){
    file = fopen("parameters.txt", "w");
    printf("Opening file\n");
    }
    
    #else
    if(!file)
    file = fopen("parameters.txt", "r");
    #endif
    long int moments[OBJECTS_TO_FIND*6];
    double mu[OBJECTS_TO_FIND*3]; //Only second order centralized moments. Order: mu20, mu11, mu20
    double nu[OBJECTS_TO_FIND*3]; //Only second order normalized moments. Order: nu20, nu11, nu20
    double phi[OBJECTS_TO_FIND*2]; //Only first two moment invariant. Order: phi1, phi2
    double X[OBJECTS_TO_FIND*3]; //Centroid and angle:
    getRegions(filteredImage, moments);

    for(int i=0; i<OBJECTS_TO_FIND; i++){ //calculate m, mu, nu and phi of each object
        
        long int m00 = moments[i*6+0];
        long int m10 = moments[i*6+1];
        long int m01 = moments[i*6+2];
        long int m20 = moments[i*6+3];
        long int m11 = moments[i*6+4];
        long int m02 = moments[i*6+5];
        
        
        
        float x_bar = m10/m00;
        float y_bar = m01/m00;
        double mu20 = m20 - x_bar*m10;
        double mu02 = m02 - y_bar*m01;
        double mu11 = m11 - x_bar*y_bar*m00;

        double nu20 = mu20/(m00*m00);
        double nu02 = mu02/(m00*m00);
        double nu11 = mu11/(m00*m00);

        
        
        double phi1 = nu20 + nu02;
        double phi2 = pow(nu20-nu02, 2)+4*nu11*nu11;
        
        phi[2*i] = phi1;
        phi[2*i+1] = phi2;
        float theta = atan2(2*mu11, mu20 - mu02)/2;
        
        circle(selectionImg, Point((int)x_bar, (int)y_bar), 5, Scalar(0, 0,255), -1, 8, 0);
        arrowedLine(selectionImg, Point((int)x_bar, (int)y_bar), Point((int)(x_bar + 50 * cos(theta)), (int)(y_bar + 50 * sin(theta))),Scalar(0, 0, 255), 1, 8, 0, .1);

        mu[i * 3 + 0] = mu20;
        mu[i * 3 + 1] = mu11;
        mu[i * 3 + 2] = mu02;

        nu[i * 3 + 0] = nu20;
        nu[i * 3 + 1] = nu11;
        nu[i * 3 + 2] = nu02;

        phi[i * 2 + 0] = phi1;
        phi[i * 2 + 1] = phi2;

        X[i * 3 + 0] =  x_bar;
        X[i * 3 + 1] =  y_bar;
        X[i * 3 + 2] =  theta;
        
        
        
    }
    //*r = 0;
    imshow("Centroid", selectionImg);
    #ifdef TRAIN

        static double phi_1[SAMPLES_PER_OBJECT];
        static double phi_2[SAMPLES_PER_OBJECT];
        if(next_sample){
            next_sample = 0;
            phi_1 [samples%SAMPLES_PER_OBJECT]= phi[0];
            phi_2 [samples%SAMPLES_PER_OBJECT]= phi[1];
            samples++;
        }
        
        if(samples%SAMPLES_PER_OBJECT==0 && samples!=0 ){
            if(x){
            x = 0;
            destroyAllWindows();
            double phi1_avg, phi2_avg, phi1_dev, phi2_dev;
            for(int i=0; i<SAMPLES_PER_OBJECT; i++){
                phi1_avg+=phi_1[i]/SAMPLES_PER_OBJECT;
                phi2_avg+=phi_2[i]/SAMPLES_PER_OBJECT;
            }
            for(int i=0; i<SAMPLES_PER_OBJECT; i++){
                phi1_dev += pow(phi_1[i]-phi1_avg, 2);
                phi2_dev += pow(phi_2[i]-phi2_avg, 2);
            }
            phi1_dev /= SAMPLES_PER_OBJECT;
            phi2_dev /=SAMPLES_PER_OBJECT;
            phi1_dev = sqrt(phi1_dev);
            phi2_dev = sqrt(phi2_dev);
            printf("%f", phi1_avg);
            fprintf(file, "%f\n", phi1_avg);
            fprintf(file, "%f\n", phi2_avg);
            fprintf(file, "%f\n", phi1_dev);
            fprintf(file, "%f\n", phi2_dev);

            //Order: phi_1 average, phi_2 average, phi_1 deviation, phi_2 deviation
            }
        } else{
            x=1;
        }
        
    #else
        static float parameters[4*OBJECTS_TO_TRAIN];
        float distances[OBJECTS_TO_TRAIN*OBJECTS_TO_FIND];
        char o1, o2, o3, o4;
        
        double angle;
        int index, long_object;
        int x, y;
        Mat mira(200, 200, CV_8UC1, Scalar(0));
        Mat phis(500, 500, CV_8UC1, Scalar(0));
        
        o1 = o2 = o3 = o4 = 0;
        
        //Read the parameters from the file
        for(int i = 0; i<OBJECTS_TO_TRAIN*4; i++){
            fscanf(file, "%f", &parameters[i]);
        }

        //Calculate the distance from each found object to each trained object
        for(int i = 0; i<OBJECTS_TO_FIND; i++){
            for(int j = 0; j<OBJECTS_TO_TRAIN; j++){
                distances [i*OBJECTS_TO_FIND + j] = pow((phi[i * 2]- parameters[4 * j])/parameters[ 4 * j + 2], 2) + 
                                                    pow((phi[i * 2 + 1]- parameters[ 4 * j + 1])/parameters[4 * j + 3], 2);
            }
        }
        
        //Get the minimum distance and set the corresponding flag
        for(int i = 0; i<OBJECTS_TO_FIND; i++){
            double min = distances[ i *OBJECTS_TO_TRAIN ];
            index = 0;
            for(int j = 1; j<4; j++){
                if(distances[i*OBJECTS_TO_TRAIN + j]< min){
                    min = distances[i * OBJECTS_TO_TRAIN + j];
                    index = j;
                }
            }

            //Check if difference with trained object is less than the deviation of the samples
            if(abs(phi[i*2]-parameters[index*4]) <= parameters[index*4 + 2] &&  abs(phi[i*2+1]-parameters[index*4+1]) <= parameters[index*4 +3]){
                    switch(index){
                        case 0:
                            o1 = 1;
                            long_object = i;
                            break;
                        case 1:
                            o2 = 1;
                            long_object = i;
                            break;
                        case 2:
                            o3 = 1;
                            break;
                        case 3:
                            o4 = 1;
                            break;
                    }
            }
        }
        
        //Get the correct coordinates
        if(o1 && o3){           // Espada y escudo
            x = 100; y = 0;
        }else if(o2 && o3){     // Lanza y escudo
            x = 0; y = 0;
        }else if(o1 && o4){     // Espada y casco
            x = 0; y = 100;
        }else if(o2 && o4){     // Lanza y casco
            x = 100; y = 100;
        }
        if(long_object<=0 && long_object<OBJECTS_TO_FIND)
            angle = X[long_object*3+2];
        else
            angle = (M_PI/3);
        printf("%f\n", angle);

        double x_max, y_max;
        x_max = parameters[0];
        y_max = parameters[1];
        for(int i=1; i<OBJECTS_TO_TRAIN; i++){
            if(parameters[4*i]>x_max)
                x_max = parameters[4*i];
            if(parameters[4*i+1]>y_max)
                y_max = parameters[4*i+1];
        }

        for(int i=0; i<OBJECTS_TO_TRAIN; i++){ //plot trained phis
            double x = 450 * parameters[4*i]/x_max;
            double y = 500 - 450 * parameters[4*i+1]/y_max;
            circle(phis, Point((int)x, (int) y),5, Scalar(255), -1, 8, 0);
        }

        for(int i = 0; i<OBJECTS_TO_FIND; i++){ //plot found phis
            double x = 450 * phi[2*i]/x_max;
            double y = 500 - 450 * phi[2*i+1]/y_max;
            circle(phis, Point((int)x, (int) y),5, Scalar(128), -1, 8, 0);
        }
        imshow("Phi", phis);
        rectangle(mira, Rect(x, y, 100, 100), Scalar(255), -1, 8, 0);
        arrowedLine(mira, Point(100, 100), Point((int)(100 + cos(angle) * 50), (int)(100 + sin(angle) * 50)), Scalar(128), 1, 8, 0, .1);
        imshow("Mira", mira);
    #endif
    }
    
    
}

void mouseCoordinatesExampleCallback(int event, int x, int y, int flags, void *param)
{
    switch (event)
    {
    case CV_EVENT_LBUTTONDOWN:
        cout << "  Mouse X, Y: " << x << ", " << y;
        cout << endl;
        /*  Draw a point */
        p.x = x;
        p.y = y;
        q.x = x;
        q.y = y;
        break;
    case CV_EVENT_MOUSEMOVE:
        break;
    case CV_EVENT_LBUTTONUP:
        q.x = x;
        q.y = y;

        break;
    }
}

void histogram(const Mat &image, char t, unsigned char *threshold, char o, char colorSpace)
{
    int **channels = (int **)malloc(sizeof(int *) * image.channels());
    Mat *mats = new Mat[image.channels()];
    unsigned char *pv = new unsigned char[image.channels()];
    memset(pv, 0, sizeof(unsigned char) * image.channels());
    if (p.x != 0 && p.y != 0)
    {

        //Get the values of the three channels from the point
        for (int i = 0; i < image.channels(); i++)
            pv[i] = image.at<Vec3b>(p.y, p.x)[i];
    }

    for (int i = 0; i < image.channels(); i++)
    {
        channels[i] = (int *)malloc(sizeof(int) * 255);
        memset(channels[i], 0, 255 * sizeof(int));
        mats[i].create(306, 512, CV_8UC1);
        mats[i] = Mat::zeros(306, 512, CV_8UC1);
    }

    for (int i = 0; i < image.rows; i++)
        for (int j = 0; j < image.cols; j++)
            for (int k = 0; k < image.channels(); k++)
                channels[k][image.at<Vec3b>(i, j)[k]]++;

    for (int i = 0; i < image.channels(); i++)
    {
        int max = 0;
        for (int j = 0; j < 255; j++)
            if (channels[i][j] > max)
                max = channels[i][j];

        for (int j = 0; j < 254; j++)
        {
            mats[i].at<char>((max - channels[i][j]) * 255 / max, j * 2) = 0xFF;
            mats[i].at<char>((max - (channels[i][j] + channels[i][j + 1]) / 2) * 255 / max, j * 2 + 1) = 0xff;
            for (int k = 0; k < 40; k++)
            {
                mats[i].at<char>(265 + k, j * 2) = j;
                mats[i].at<char>(265 + k, j * 2 + 1) = j;
            }
        }

        //Plot a vertical bar in the position of the values of the point
        for (int j = 0; j < 255; j++)
        {
            mats[i].at<char>(j, pv[i] * 2) = 0x7f;
            mats[i].at<char>(j, pv[i] * 2 + 1) = 0x7f;
        }

        if (o)

            //Plot a vertical bar at the minimum and maximun of the ROI
            for (int j = 0; j < 255; j++)
            {
                mats[i].at<char>(j, threshold[i] * 2) = 0x3f;
                mats[i].at<char>(j, threshold[i] * 2 + 1) = 0x3f;
                mats[i].at<char>(j, threshold[i + image.channels()] * 2) = 0xff;
                mats[i].at<char>(j, threshold[i + image.channels()] * 2 + 1) = 0xff;
            }


        char title[20];
        switch(colorSpace) {
        case 0: //BGR
            switch (i){
                case 0: 
                sprintf(title, "Blue");
                break;
                case 1:
                sprintf(title, "Green");
                break;
                case 2:
                sprintf(title, "Red");
                break;
            }
            break;
        case 1:
            
                for (int j = 0; j < 255; j++)
            {
                mats[i].at<char>(j, thVal*2) = 0xf0;
                mats[i].at<char>(j, thVal*2 + 1) = 0xf0;
                
            }
                sprintf(title, "GreyScale");
                
            
            break;
        case 2:
            switch (i){
                case 0: 
                sprintf(title, "Y");
                break;
                case 1:
                sprintf(title, "In-phase");
                break;
                case 2:
                sprintf(title, "Quadrature");
                break;
            }
            break;
        case 3:
            switch (i){
                case 0: 
                sprintf(title, "Hue");
                break;
                case 1:
                sprintf(title, "Saturation");
                break;
                case 2:
                sprintf(title, "Value");
                break;
            }
            break;
    }
        
        imshow(title, mats[i]);
    }

    for (int i = 0; i < image.channels(); i++)
    {
        free(channels[i]);
        mats[i].deallocate();
    }
    free(channels);
    delete[] mats;
}

void bgr2yiq(Mat &image)
{
    uchar r, g, b;
    double y, i, q;

    Mat yiq = Mat::zeros(image.rows, image.cols, CV_8UC3);
    /* Convert image from RGB to YIQ */

    for (int row = 0; row < image.rows; ++row)
    {
        for (int col = 0; col < image.cols; ++col)
        {
            r = image.at<Vec3b>(row, col)[2];
            g = image.at<Vec3b>(row, col)[1];
            b = image.at<Vec3b>(row, col)[0];

            y = 0.299 * r + 0.587 * g + 0.114 * b;
            i = 0.596 * r - 0.275 * g - 0.321 * b;
            q = 0.212 * r - 0.523 * g + 0.311 * b;

            i /= .596 * 2;
            q /= .523 * 2;
            i += 128;
            q += 128;
            yiq.at<Vec3b>(row, col)[0] = y;
            yiq.at<Vec3b>(row, col)[1] = i;
            yiq.at<Vec3b>(row, col)[2] = q;
        }
    }
    
    image = yiq;
}

void bgr2gs(Mat &image)
{
    
    Mat gsImage(image.rows, image.cols, CV_8UC1, Scalar(0));

    if (image.data)
    {

        uchar r, g, b;
        double gs;
        uchar val;

        /* Convert image from RGB to Grayscale */

        for (int row = 0; row < image.rows; ++row)
        {
            for (int col = 0; col < image.cols; ++col)
            {
                r = image.at<Vec3b>(row, col)[2];
                g = image.at<Vec3b>(row, col)[1];
                b = image.at<Vec3b>(row, col)[0];

                gs = 0.299 * r + 0.587 * g + 0.114 * b;
                val = (uchar)gs;
                gsImage.at<uchar>(row, col) = val;
                
               /* gsImage.at<Vec3b>(row, col)[2] = gs;
                gsImage.at<Vec3b>(row, col)[1] = gs;
                gsImage.at<Vec3b>(row, col)[0] = gs;*/
            }
        }
        //gsImage.convertTo(gsImage, CV_8UC3);
        namedWindow("Threshold", CV_WINDOW_AUTOSIZE);
        Mat th(100,512, CV_8UC(1), Scalar(255));
        for (int i = 0; i<256; i++)
            for(int j = 0; j<100; j++){
                th.at<char>(j, 2*i) = i;
                th.at<char>(j, 2*i+1) = i;
            }
         for (int j = 0; j < 100; j++)
            {
                th.at<char>(j, thVal*2) = 0xf0;
                th.at<char>(j, thVal*2 + 1) = 0xf0;
                
            }
        imshow("Threshold", th);
        setMouseCallback("Threshold", mouseCallback);
        Mat binary = gsImage.clone();
        for (int row = 0; row<gsImage.rows; row++)
            for (int col = 0; col<gsImage.cols; col++)
                    if(gsImage.at<uchar>(row, col)>thVal)
                        binary.at<uchar>(row,col) = 0xFF;
                    else
                        binary.at<uchar>(row, col) = 0;
         
        imshow("Binary", binary);
                    
        image = gsImage;

    }
}

void bgr2hsv(Mat &image)
{
    if (image.data)
    {
        Mat hsvImage = image.clone();
        cvtColor(image, hsvImage, CV_BGR2HSV);
        //imshow("HSV image", hsvImage);

        image = hsvImage;
    }
}
