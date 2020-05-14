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

#define OBJECTS_TO_FIND 2

using namespace cv;
using namespace std;

// Here we will store points
Point p(0, 0), q(0, 0);
unsigned char thVal = 0;

void getRegions(Mat &image, long int * ordinary_moments ){
    Mat color(image.rows, image.cols, CV_8UC1, Scalar(0));
    
    
    memset(ordinary_moments, 0, OBJECTS_TO_FIND*sizeof(long int)*6);
    

    for(int k = 1; k<OBJECTS_TO_FIND+1; k++){
        printf("\n");
        int x, y;
        do{
        x = rand()%image.cols;
        y = rand()%image.rows;
        
        }while(image.at<uchar>(y, x) ==0 || color.at<uchar>(y, x)!=0);
        
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

void selection(Mat image, unsigned char *threshold, Mat original);



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
    char t = 0, h = 0, o = 0, b = 0, v = 0, y = 0, f = 0;
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
                selection(currentImage, threshold, originalImage);
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
            }
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

void selection(Mat image, unsigned char *threshold, Mat original)
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
    imshow("Selection", filteredImage);
    long int moments[OBJECTS_TO_FIND*6];
    getRegions(filteredImage, moments);

    for(int i=0; i<OBJECTS_TO_FIND; i++){
        
        printf("%d M00: %ld\t",i, moments[i*6+0]);
        printf("%d M10: %ld\t",i, moments[i*6+1]);
        printf("%d M01: %ld\t",i, moments[i*6+2]);
        printf("%d M20: %ld\t",i, moments[i*6+3]);
        printf("%d M11: %ld\t",i, moments[i*6+4]);
        printf("%d M02: %ld\t",i, moments[i*6+5]);
        printf("\n");
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