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



using namespace cv;
using namespace std;

// Here we will store points
Point p(0, 0), q(0, 0);




/* This is the callback that will only display mouse coordinates */
void mouseCoordinatesExampleCallback(int event, int x, int y, int flags, void *param);

void histogram(const Mat &image, char t, unsigned char* threshold, char o);

void bgr2gs(const Mat &image);

void draw(Mat image, char t, char *o, unsigned char*threshold){
    Mat roi;
    Mat img = image.clone();
    if(!t){
        if(p.x>0 && p.y>0)
            circle(img, p, 5, Scalar(255, 255, 255), 1, 8, 0);
    }
    else
    {
        if(p!=q){//Valid ROI selection
                Rect rec(p, q);
                roi = image(rec);
                imshow("ROI", roi);
                rectangle(img, p, q, Scalar(0, 0, 255), 1, 8, 0);
                
                *o = 1;
                //Calculate the maxima and minima of the ROI 
                for (int r = 0; r < roi.rows; r++)
                    for (int c = 0; c < roi.cols; c++)
                        for (int ch = 0; ch < roi.channels(); ch++)
                        {
                            if (roi.at<Vec3b>(r, c)[ch] < threshold[ch])
                                threshold[ch] = roi.at<Vec3b>(r, c)[ch];

                            if (roi.at<Vec3b>(r, c)[ch] > threshold[ch + roi.channels()])
                                threshold[ch + roi.channels()] = roi.at<Vec3b>(r, c)[ch];
                        }
                }
    }
    
    imshow("Image", img);

}
int main(int argc, char *argv[])
{
    /* First, open camera device */
    
    VideoCapture camera;
    camera.open(0);

    /* Create images where captured and transformed frames are going to be stored */
    Mat currentImage;
    
    /* Create main OpenCV window to attach callbacks */
    namedWindow("Image");

    setMouseCallback("Image", mouseCoordinatesExampleCallback);
    char t = 0, h = 0, o = 0, b = 0, i = 1;
    char x;
    camera>>currentImage;
    unsigned char *threshold = new unsigned char[currentImage.channels() * 2];
            memset(threshold, 0xff, currentImage.channels() * sizeof(char));
            memset(threshold + currentImage.channels() * sizeof(char), 0x00, currentImage.channels() * sizeof(char));
    while (true)
    {
        /* Obtain a new frame from camera */
        if (!t)
            camera >> currentImage;

        if (currentImage.data)
        {   
            
            /* Show image */
            
            draw(currentImage, t, &o, threshold);

            //Show histograms
            if (h)
                histogram(currentImage, t, threshold, o);
            else{
                destroyWindow("Channel 0");
                destroyWindow("Channel 1");
                destroyWindow("Channel 2");
            }
            
            if (b)
                bgr2gs(currentImage);
            else {
                destroyWindow("GRAYSCALE_IMAGE");
            }

            /* If 'x' is pressed, exit program */
            
            
            x = waitKey(3);
            if (x == ' ') //Still image
            {
                if(!t)
                    q = p;
                t = ~t;
            } else if(x=='x') {
                break; //End program
            } else if(x=='h') {   
                h = ~h; //Show histogram
            } else if(x == 'b') {
                b = ~b;
            }
             
        }
        else
        {
            cout << "No image data.. " << endl;
        }
    }
    
      delete[] threshold; 
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

void histogram(const Mat &image, char t, unsigned char*threshold, char o)
{
    int **channels = (int **)malloc(sizeof(int *) * image.channels());
    Mat *mats = new Mat[image.channels()];
    unsigned char * pv =  new unsigned char[image.channels()];
    memset(pv, 0, sizeof(unsigned char)*image.channels());
    if (p.x != 0 && p.y != 0)
    {

        
        //Get the values of the three channels from the point
        for(int i=0; i<image.channels(); i++)
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
            mats[i].at<char>(j, threshold[i+image.channels()] * 2) = 0xff;
            mats[i].at<char>(j, threshold[i+image.channels()] * 2 + 1) = 0xff;
            
        }


        char title[20];
        sprintf(title, "Channel %d", i);
        imshow(title, mats[i]);
    }

    for (int i = 0; i < image.channels(); i++)
    {
        free(channels[i]);
        mats[i].deallocate();
    }
    free(channels);
    delete [] mats;
    
}

void bgr2gs(const Mat &image) {
    Mat gsImage = Mat::zeros(image.rows,image.cols, CV_8UC1);
 
    if( image.data ) { 
 
    for (int i=0; i<image.cols ; i++) {
        for (int j=0 ; j<image.rows ; j++) { 
            Vec3b color1 = image.at<Vec3b>(Point(i,j));
            Scalar color2 = gsImage.at<uchar>(Point(i,j));
            color2 = (0.11*color1.val[0]+0.59*color1.val[1]+0.3*color1.val[2]);
    
            gsImage.at<uchar>(Point(i,j)) = color2.val[0];
        }
    }
    
    namedWindow("GRAYSCALE_IMAGE",CV_WINDOW_AUTOSIZE); 
    imshow("GRAYSCALE_IMAGE", gsImage); 
    
    }
}