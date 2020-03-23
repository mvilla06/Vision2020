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

void bgr2gs(Mat &image);
void bgr2hsv(Mat &image);
void bgr2yiq(Mat &image);

void draw(Mat image, char t, char *o, unsigned char*threshold);

void selection(Mat image, unsigned char*threshold, int colorSpace);



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
  

    // Flags: t -> Freeze camera, h -> Build histograms, b -> Convert grayscale, y-> Convert YIQ
    // v -> Convert to HSV
    char t = 0, h = 0, o = 0, b = 0, v = 0, y = 0, i = 1, f = 0;
    int colorSpace = 0; // 0->BGR, 1->BW, 2->YIQ, 3->HSV
    char x;
    camera >> currentImage;

    unsigned char *threshold = new unsigned char[currentImage.channels() * 2];
    
    
    while (true)
    {

        

        /* Obtain a new frame from camera */
        if (!t) {
            camera >> currentImage;
            
        }

        if (currentImage.data)
        {   
            //Convert Image only with live feed
            if(!t){
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


            if(!o){
                destroyWindow("ROI");
                destroyWindow("Selection");
            }
            /* Show image */
            draw(currentImage, t, &o, threshold);

            if (o && f){
                selection(currentImage, threshold, colorSpace);
            }
            else
                destroyWindow("Selection");
            
            //Show histograms
            if (h)
                histogram(currentImage, t, threshold, o);
            else{
                destroyWindow("Channel 0");
                destroyWindow("Channel 1");
                destroyWindow("Channel 2");
            }

            
            imshow("Image", currentImage);
            x = waitKey(3);
            switch (x)
            {
            case ' ': //Freeze image
                if (!t)
                    q = p;
                t = ~t; 
                colorSpace = 0;
                break;
            case 'x': //End program
                goto end;
                break;
            case 'h': //Show histograms
                h = ~h;
                break;
            case 'b': //Greyscale 
                b = ~b;
                colorSpace = 1;
                o = 0;
                break;
            case 'y': //YIQ
                y = ~y; 
                colorSpace = 2;
                o = 0;
                break;
            case 'v': //HSV
                v = ~v;
                colorSpace = 3;
                o = 0;
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

void selection(Mat image, unsigned char*threshold, int colorSpace) {
    Mat selectionImg = image.clone();

    for (int r = 0; r < selectionImg.rows; r++) {
        for (int c = 0; c < selectionImg.cols; c++) {
            bool totalCond = true;
            for (int ch = 0; ch < selectionImg.channels(); ch++) {
                bool cond = selectionImg.at<Vec3b>(r, c)[ch] > threshold[ch] && selectionImg.at<Vec3b>(r,c)[ch] < threshold[ch + selectionImg.channels()];
                totalCond = totalCond && cond;
            }
            if (!totalCond) {
                for (int ch = 0; ch < selectionImg.channels(); ch++) {
                    selectionImg.at<Vec3b>(r, c)[ch] = 0;
                }
            }
        }
    }
    imshow("Selection", selectionImg);
    
    /*switch(colorSpace) {
        case 0:
            imshow("Selection BGR", selectionImg);
            break;
        case 1:
            imshow("Selection BW", selectionImg);
            break;
        case 2:
            imshow("Selection YIQ", selectionImg);
            break;
        case 3:
            imshow("Selection HSV", selectionImg);
            break;
    }*/
    
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

void bgr2yiq(Mat &image) {
    uchar r, g, b;
	double y, i, q;

	Mat yiq = Mat::zeros(image.rows, image.cols, CV_8UC3);
	/* Convert image from RGB to YIQ */

    for (int row = 0; row < image.rows; ++row) {
		for (int col = 0; col < image.cols; ++col){
            r = image.at<Vec3b>(row, col)[2];
            g = image.at<Vec3b>(row, col)[1];
            b = image.at<Vec3b>(row, col)[0];
			
            y = 0.299*r + 0.587*g + 0.114*b;
			i = 0.596*r - 0.275*g - 0.321*b;
			q = 0.212*r - 0.523*g + 0.311*b;

            i/=.596*2;
            q/=.523*2;
            i+=128;
            q+=128;
            yiq.at<Vec3b>(row, col)[0] = y;
            yiq.at<Vec3b>(row, col)[1] = i;
            yiq.at<Vec3b>(row, col)[2] = q;
        }
    }   
    //yiq.convertTo(yiq, CV_8U);
    //namedWindow("YIQ IMAGE", CV_WINDOW_AUTOSIZE);
    //imshow("YIQ IMAGE", yiq);

    image = yiq;

}

void bgr2gs(Mat &image) {
    Mat gsImage = Mat::zeros(image.rows,image.cols, CV_8UC3);
 
    if( image.data ) { 
 
    uchar r, g, b;
	double gs;

	/* Convert image from RGB to Grayscale */

    for (int row = 0; row < image.rows; ++row) {
		for (int col = 0; col < image.cols; ++col){
            r = image.at<Vec3b>(row, col)[2];
            g = image.at<Vec3b>(row, col)[1];
            b = image.at<Vec3b>(row, col)[0];
			
            gs = 0.299*r + 0.587*g + 0.114*b;

            gsImage.at<Vec3b>(row, col)[2] = gs;
            gsImage.at<Vec3b>(row, col)[1] = gs;
            gsImage.at<Vec3b>(row, col)[0] = gs;
        }
    }   
    gsImage.convertTo(gsImage, CV_8UC3);
    //namedWindow("GRAYSCALE IMAGE", CV_WINDOW_AUTOSIZE);
    //imshow("GRAYSCALE IMAGE", gsImage);

    image = gsImage;
    }
}

void bgr2hsv(Mat &image) {
    if (image.data) {
        Mat hsvImage = image.clone();
        cvtColor(image, hsvImage, CV_BGR2HSV);
        //imshow("HSV image", hsvImage);

        image = hsvImage;
    }
}