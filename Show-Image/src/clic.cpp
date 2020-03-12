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
Point p(0, 0);

/* This is the callback that will only display mouse coordinates */
void mouseCoordinatesExampleCallback(int event, int x, int y, int flags, void* param);

void histogram(const Mat &image, int * b, int *g, int *r);



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
    char t = 0;
    while (true)
	{
		/* Obtain a new frame from camera */
        if(!t)
            camera >> currentImage;
            
        
       
		    

		if (currentImage.data) 
		{
            

            int b[255], g[255], r[255];
            
            histogram(currentImage, b, g, r);

            
            /* Show image */
            imshow("Image", currentImage);

			/* If 'x' is pressed, exit program */
			/*if (waitKey(3) == 'x')
				break;*/
            
             if (waitKey(3)=='p'){
                if(t)
                t = 0;
                else
                t=1;
                
            }
                
		}
		else
		{
			cout << "No image data.. " << endl;
		}
	}
}


void mouseCoordinatesExampleCallback(int event, int x, int y, int flags, void* param)
{
    switch (event)
    {
        case CV_EVENT_LBUTTONDOWN:
            cout << "  Mouse X, Y: " << x << ", " << y ;
            cout << endl;
            /*  Draw a point */
            p.x = x; p.y = y;
            break;
        case CV_EVENT_MOUSEMOVE:
            break;
        case CV_EVENT_LBUTTONUP:
            break;
    }
}


void histogram(const Mat &image, int * b, int * g, int * r)
{
    for (int i=0; i<255; i++) {
        b[i] = 0;
        g[i] = 0;
        r[i] = 0;
    }
    
    for (int i = 0; i< image.rows; i++){
        for (int j = 0; j<image.cols ; j++){
            
            b[image.at<Vec3b>(i, j)[0]]++;
            g[image.at<Vec3b>(i, j)[1]]++;
            r[image.at<Vec3b>(i, j)[2]]++;  
        }
    }

    int max = 0;
    for (int i = 0; i<255; i++)
        if (b[i]>max) max= b[i];
    
    Mat blue(306, 512, CV_8UC(1), Scalar(0));
    for (int i=0; i<254; i++){
        blue.at<char>((max-b[i])*255/max, i*2) = 0xFF;
        blue.at<char>((max-(b[i]+b[i+1])/2)*255/max, i*2+1) = 0xff;
        for(int j = 0; j<40; j++){
            blue.at<char>(265+j, i*2) = i;
            blue.at<char>(265+j, i*2+1) = i;
        }
    }
    

    max = 0;
    for (int i = 0; i<255; i++)
        if (g[i]>max) max= g[i];
    
    Mat green(306, 512, CV_8UC(1), Scalar(0));
    for (int i=0; i<254; i++){
        green.at<char>((max-g[i])*255/max, i*2) = 0xFF;
        green.at<char>((max-(g[i]+g[i+1])/2)*255/max, i*2+1) = 0xff;
        for(int j = 0; j<40; j++){
            green.at<char>(265+j, i*2) = i;
            green.at<char>(265+j, i*2+1) = i;
        }
    }
    

    max = 0;
    for (int i = 0; i<255; i++)
        if (r[i]>max) max= r[i];
    
    Mat red(306, 512, CV_8UC(1), Scalar(0));
    for (int i=0; i<254; i++){
        red.at<char>((max-r[i])*255/max, i*2) = 0xFF;
        red.at<char>((max-(r[i]+r[i+1])/2)*255/max, i*2+1) = 0xff;
        for(int j = 0; j<40; j++){
            red.at<char>(265+j, i*2) = i;
            red.at<char>(265+j, i*2+1) = i;
        }
    }

    unsigned char pb = 0, pg = 0, pr = 0;
    
    if(p.x!= 0 && p.y!=0){

        circle(image,p, 5, Scalar(255, 0, 0), 1, 8, 0);
        //Get the values of the three channels from the point
        pb = image.at<Vec3b>(p.y, p.x)[0];
        pg = image.at<Vec3b>(p.y, p.x)[1];
        pr = image.at<Vec3b>(p.y, p.x)[2];

        //Plot a vertical bar in the position of the values of the point
        for (int i=0; i<255; i++){
            blue.at<char>(i, pb*2) = 0xff;
            green.at<char>(i, pg*2) = 0xff;
            red.at<char>(i, pr*2) = 0xff;
            blue.at<char>(i, pb*2+1) = 0xff;
            green.at<char>(i, pg*2+1) = 0xff;
            red.at<char>(i, pr*2+1) = 0xff;
        }
    }
    
    
    imshow("Blue", blue);
    imshow("Green", green);
    imshow("Red", red);
}