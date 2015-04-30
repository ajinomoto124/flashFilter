#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>



/*
* convolve with a box filter
*/
Image* ip_blur_box (Image* src, int size)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    double kernel[size*size];//create kernel with width and height = size
    
    //loop through kernel
    for (int i = 0; i<size; ++i) {
        for (int j = 0; j<size; ++j) {
            kernel[j*size+i]= 1.0/(size*size);//set each value of the kernel to 1/size*size
        }
    }
    
    tempImage = ip_convolve(src, size, kernel); //call convolve with image, size and kernel
    
    return tempImage;
    
}


/*
* convolve with a gaussian filter
*/
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    double kernel[size*size];//create kernel with width and height = size
    
    double accumulator = 0.0;//accumulator value saves the total sum of all values in the kernel
    
    //loop through kernel
    for (int i = 0; i<size; ++i) {
        for (int j = 0; j<size; ++j) {
            kernel[j*size+i]= exp((-((i*i) + (j*j)))/(2.0*(sigma*sigma)));//set each value to guassian formula of that pixel
            accumulator += kernel[j*size+i];//add value to accumulator
        }
    }
    //loop through kernel again
    for (int i = 0; i<size; ++i) {
        for (int j = 0; j<size; ++j) {
            kernel[j*size+i] /= accumulator;//divide each value of the kernel by the total sum
        }
    }
    
    tempImage = ip_convolve(src, size, kernel);//call convolve using guassian kernel
    
    return tempImage;
}


/*
* convolve with a triangle filter
*/
Image* ip_blur_triangle (Image* src, int size)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* interpolate with a black image
*/
Image* ip_brighten (Image* src, double alpha)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
	Image* blackImage = new Image(src->getWidth(), src->getHeight());//image with the same size that is all black
    
    //loop through black image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            blackImage->setPixel(i, j, 0, 0); //set each channel to 0 for black
            blackImage->setPixel(i, j, 1, 0);
            blackImage->setPixel(i, j, 2, 0);
        }
    }
    
    tempImage = ip_interpolate(src, blackImage, alpha);//call interpolate using image and black image
    
    
    return tempImage;
}


/*
* shift colors
*/
Image* ip_color_shift(Image* src)
{
  
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    
    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            tempImage->setPixel(i, j, 0, src->getPixel_(i, j, 1)); //shift red to green
            tempImage->setPixel(i, j, 1, src->getPixel_(i, j, 2)); //shift green to blue
            tempImage->setPixel(i, j, 2, src->getPixel_(i, j, 0)); //shift blue to red
        }
    }
    
    return tempImage;
}


/*
* use a mask image for a per-pixel alpha value to perform
* interpolation with a second image
*/
Image* ip_composite3 (Image* src1, Image* src2,
					 Image* mask)
{
    Image* tempImage = new Image(src1->getWidth(), src2->getHeight());
    
    //loop through image
    for (int i = 0; i < src2->getWidth(); ++i) {
        for (int j = 0; j < src2->getHeight(); ++j) {
            tempImage->setPixel_(i, j, 0, (mask->getPixel(i,j,0)) * src1->getPixel_(i, j, 0) + (1-(mask->getPixel(i,j,0))) * src2->getPixel_(i, j, 0)); //combine pixel with second image's pixel based on interpolation formula
            tempImage->setPixel_(i, j, 1, (mask->getPixel(i,j,1)) * src1->getPixel_(i, j, 1) + (1-(mask->getPixel(i,j,1))) * src2->getPixel_(i, j, 1));
            tempImage->setPixel_(i, j, 2, (mask->getPixel(i,j,2)) * src1->getPixel_(i, j, 2) + (1-(mask->getPixel(i,j,2))) * src2->getPixel_(i, j, 2));
        }
    }
    return tempImage;
}



/*
* interpolate with the average intensity of the src image
*/
Image* ip_contrast (Image* src, double alpha)
{
    
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    Image* greyImage = new Image(src->getWidth(), src->getHeight());//create grey image of same size
    
    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            greyImage->setPixel_(i, j, 0, 0.5); //set each channel value to 0.5
            greyImage->setPixel_(i, j, 1, 0.5);
            greyImage->setPixel_(i, j, 2, 0.5);
        }
    }
    
    tempImage = ip_interpolate(src, greyImage, alpha);//call interpolate using image and grey image
    
    return tempImage;
}



Image* ip_composite (Image* src1, Image* src2,
                      Image* mask)
{
    Image* tempImage = new Image(src1->getWidth(), src2->getHeight());
    
    //get mask M with sg_mask(no flash, flash)
    //get Anr with sg_jBilateral(no flash, flash)
    //get fBase with sg_jBolateral(flash, flash);
    //get fDetail with sg_detail(flash, fBase);
    //get Abase with sg_jBilateral(noflash, noflash)
    //𝐴𝐹𝑖𝑛𝑎𝑙 = (1 − 𝑀)(𝐴𝑁𝑅)𝐹𝐷𝑒𝑡𝑎𝑖𝑙 + (𝑀)𝐴𝐵𝑎𝑠𝑒
    
    return tempImage;
}

Image* sg_mask (Image* src1, Image* src2,
					 Image* mask)
{
    int width = src1->getWidth();
    int height = src1->getHeight();
    Image* tempImage = new Image(width, height);
    Image* mask1 = new Image(width, height);
    
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float pix1 = src1->getPixel_(x, y, 0)+src1->getPixel_(x, y, 1)+src1->getPixel_(x, y, 2);
            float pix2 = src2->getPixel_(x, y, 0)+src2->getPixel_(x, y, 1)+src2->getPixel_(x, y, 2);
            if(pix2-pix1 <= 0.5){
                mask1->setPixel_(x,y, 0, 1);
                mask1->setPixel_(x,y, 1, 1);
                mask1->setPixel_(x,y, 2, 1);
            }
            else{
                mask1->setPixel_(x,y, 0, 0);
                mask1->setPixel_(x,y, 1, 0);
                mask1->setPixel_(x,y, 2, 0);
            }
        }
    }
    
    
    tempImage = ip_composite3(src1, src2, mask1);
    
    return tempImage;
}



Image* sg_jBilateral (Image* src, Image* src2, Image* mask){
    
    float id = 48.0;
    float cd = 0.1;
    
    int width = src->getWidth();
    int height = src->getHeight();
    Image* tempImage = new Image(width, height);
    
    
    int kernelDim = 9;
    int kernelHalf = 4;
    for (int y = 0; y < width; ++y) {
        for (int x = 0; x < height; ++x) {
            
            float sumWeight = 0;
            float sum[3];
            sum[0] = 0;
            sum[1] = 0;
            sum[2] = 0;
            
            int ctrIdx = y*width+x;
            
            Pixel* ctrPix = new Pixel();
            Pixel* ctrPix2 = new Pixel();
            
            ctrPix->setColor(0, src->getPixel_(y, x, 0));
            ctrPix->setColor(1, src->getPixel_(y, x, 1));
            ctrPix->setColor(2, src->getPixel_(y, x, 2));
            
            ctrPix2->setColor(0, src2->getPixel_(y, x, 0));
            ctrPix2->setColor(1, src2->getPixel_(y, x, 1));
            ctrPix2->setColor(2, src2->getPixel_(y, x, 2));
            
            for (int m = 0; m < kernelDim; ++m) {
                for (int n = 0; n < kernelDim; ++n) {
                    
                                   
                    Pixel* curPix = new Pixel();
                    Pixel* curPix2 = new Pixel();
                    
                    curPix->setColor(0, src->getPixel_(y-((kernelDim-1)/2)+m, x-((kernelDim-1)/2)+n, RED));
                    curPix->setColor(1, src->getPixel_(y-((kernelDim-1)/2)+m, x-((kernelDim-1)/2)+n, GREEN));
                    curPix->setColor(2, src->getPixel_(y-((kernelDim-1)/2)+m, x-((kernelDim-1)/2)+n, BLUE));
                    
                    curPix2->setColor(0, src2->getPixel_(y-((kernelDim-1)/2)+m, x-((kernelDim-1)/2)+n, RED));
                    curPix2->setColor(1, src2->getPixel_(y-((kernelDim-1)/2)+m, x-((kernelDim-1)/2)+n, GREEN));
                    curPix2->setColor(2, src2->getPixel_(y-((kernelDim-1)/2)+m, x-((kernelDim-1)/2)+n, BLUE));
                    
                    float currWeight;
                    
                    float imageDist = sqrt( (float)((m-x)*(m-x) + (n-y)*(n-y)) );
                    
                    float colorDist = sqrt( (float)( (curPix2->getColor(0) - ctrPix2->getColor(0))*(curPix2->getColor(0) - ctrPix2->getColor(0)) +
                                                    (curPix2->getColor(1) - ctrPix2->getColor(1))*(curPix2->getColor(1) - ctrPix2->getColor(1)) +
                                                    (curPix2->getColor(2) - ctrPix2->getColor(2))*(curPix2->getColor(2) - ctrPix2->getColor(2)) ) );
                    
                    currWeight = 1.0f/(exp((imageDist/id)*(imageDist/id)*0.5)*exp((colorDist/cd)*(colorDist/cd)*0.5));
                    sumWeight += currWeight;
                    
                    sum[0] += currWeight*curPix->getColor(0);
                    sum[1] += currWeight*curPix->getColor(1);
                    sum[2] += currWeight*curPix->getColor(2);
                }
                
            }
            if(sumWeight> 0){
            sum[0] /= sumWeight;
            sum[1] /= sumWeight;
            sum[2] /= sumWeight;
            }
            
                tempImage->setPixel_(y, x, 0, sum[0]);
                tempImage->setPixel_(y, x, 1, sum[1]);
                tempImage->setPixel_(y, x, 2, sum[2]);
            
        }
    }


    return tempImage;
}

/*
* convolve an image with a kernel
*/
Image* ip_convolve (Image* src, int size, double* kernel )
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    
    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            
            //accumulators store total channel value for current kernel
            double accumulatorRED = 0.0;
            double accumulatorBLUE = 0.0;
            double accumulatorGREEN = 0.0;
            
            //loop through kernel
            for (int m = 0; m < size; ++m) {
                for (int n = 0; n < size; ++n) {
                   
                    //kernel value at this location
                    double tempKernel = kernel[m*size+n];
 
                    //adds the value to the accumulator for each color
                    accumulatorRED += ((src->getPixel_(i-((size-1)/2)+m, j-((size-1)/2)+n, RED)) * tempKernel);
                    accumulatorGREEN += ((src->getPixel_(i-((size-1)/2)+m, j-((size-1)/2)+n, GREEN)) * tempKernel);
                    accumulatorBLUE += ((src->getPixel_(i-((size-1)/2)+m, j-((size-1)/2)+n, BLUE)) * tempKernel);
                }
            }
            
            //set each pixel of the image to the accumulator for that channel from the kernel
            tempImage->setPixel_(i, j, RED, accumulatorRED);
            tempImage->setPixel_(i, j, GREEN, accumulatorGREEN);
            tempImage->setPixel_(i, j, BLUE, accumulatorBLUE);
        }
    }
    return tempImage;
}



/*
*  create cropped version of image
*/
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
    Image* tempImage = new Image(x1-x0, y1-y0);
    for (int i = 0; i < x1-x0; ++i) {
        for (int j = 0; j < y1-y0; ++j) {
            Pixel tempPix = (ip_resample_nearest(src, i, j));
            tempImage->setPixel_(i, j, tempPix);
        }
    }
    return tempImage;
}
/*
* convolve with an edge detection kernel
*/
Image* ip_edge_detect (Image* src)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    double kernel[9];//kernel will always be 3 x 3
    
    //loop through kernel
    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; ++j) {
            kernel[j*3+i]= -1.0;    //set edge values to -1
            if(i == 1 && j == 1){
                kernel[j*3+i] = 8;  //for center, set value to 8
            }
        }
    }
    
    tempImage = ip_convolve(src, 3, kernel);//call convolve with this kernel
    
    return tempImage;
}


/*
* extract channel of input image
*/
Image* ip_extract (Image* src, int channel)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    
    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            if(channel == 0){
            tempImage->setPixel(i, j, 0, src->getPixel_(i, j, 0)); //set channel value to src channel value
            tempImage->setPixel(i, j, 1, 0);    //remove rest of the channel values
            tempImage->setPixel(i, j, 2, 0);
            }
            else if(channel == 1){
                tempImage->setPixel(i, j, 0, 0);
                tempImage->setPixel(i, j, 1, src->getPixel_(i, j, 1));
                tempImage->setPixel(i, j, 2, 0);
            }
            else if(channel == 2){
                tempImage->setPixel(i, j, 0, 0);
                tempImage->setPixel(i, j, 1, 0);
                tempImage->setPixel(i, j, 2, src->getPixel_(i, j, 2));
            }
        }
    }
    
    return tempImage;
}


/*
* create your own fun warp
*/
Image* ip_fun_warp (Image* src)
{

	//  ask user for input parameters here including resampling method and, 
	//  if gaussian resampling is used, its filtersize and sigma
	//  if you implement more than one warp, you should ask the 
	//  user to chose the one to perform here too!
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    
    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            double r = sqrt( ((i-167.5)*(i-167.5)) + ((j - 229.5)*(j - 229.5)));
            double a = atan2(j-229.5, i-167.5);
            double rn = pow(r, 1.1);
            double x = rn*cos(a) + 167.5;
            double y = rn*sin(a) + 229.5;
            Pixel tempPix = ip_resample_nearest(src, i, j);
            tempImage->setPixel_(x, y, tempPix);
//            Pixel tempPix2 = src->getPixel_(i, j);
//            tempImage->setPixel_(x, y, tempPix2);
        }
    }
	return tempImage;
}
/*
* create a new image with values equal to the psychosomatic intensities
* of the source image
*/
Image* ip_grey (Image* src)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());

    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            tempImage->setPixel(i, j, 0, (0.2126 * src->getPixel_(i, j, 0)) + (0.7152 * src->getPixel_(i, j, 1)) +(0.0722 * src->getPixel_(i, j, 2))); //set channel to average value based on formula
            tempImage->setPixel(i, j, 1, (0.2126 * src->getPixel_(i, j, 0)) + (0.7152 * src->getPixel_(i, j, 1)) +(0.0722 * src->getPixel_(i, j, 2)));
            tempImage->setPixel(i, j, 2, (0.2126 * src->getPixel_(i, j, 0)) + (0.7152 * src->getPixel_(i, j, 1)) +(0.0722 * src->getPixel_(i, j, 2)));
        }
    }
    
    return tempImage;
}


/*
*  shift image by dx and dy (modulo width & height)
*/

Image* ip_image_shift (Image* src, double dx, double dy)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    
    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            tempImage->setPixel_(((int)(i + dx) % (int)(src->getWidth())), ((int)(j + dy) % (int)(src->getHeight())), 0, src->getPixel_(i,j, 0)); //shift position of red pixel
            tempImage->setPixel_(((int)(i + dx) % (int)(src->getWidth())), ((int)(j + dy) % (int)(src->getHeight())), 1, src->getPixel_(i,j, 1));
            tempImage->setPixel_(((int)(i + dx) % (int)(src->getWidth())), ((int)(j + dy) % (int)(src->getHeight())), 2, src->getPixel_(i,j, 2));
        }
    }
    return tempImage;
}
/*
* interpolate an image with another image
*/
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
    Image* tempImage = new Image(src1->getWidth(), src2->getHeight());
    
    //loop through image
    for (int i = 0; i < src2->getWidth(); ++i) {
        for (int j = 0; j < src2->getHeight(); ++j) {
            tempImage->setPixel_(i, j, 0, alpha * src1->getPixel_(i, j, 0) + (1-alpha) * src2->getPixel_(i, j, 0)); //combine pixel with second image's pixel based on interpolation formula
            tempImage->setPixel_(i, j, 1, alpha * src1->getPixel_(i, j, 1) + (1-alpha) * src2->getPixel_(i, j, 1));
            tempImage->setPixel_(i, j, 2, alpha * src1->getPixel_(i, j, 2) + (1-alpha) * src2->getPixel_(i, j, 2));
        }
    }
    return tempImage;
}
/*
* invert input image
*/
Image* ip_invert (Image* src)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    Image* greyImage = new Image(src->getWidth(), src->getHeight());//create temp grey image
    
    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            greyImage->setPixel_(i, j, 0, 0.5); //set each channel value to 0.5
            greyImage->setPixel_(i, j, 1, 0.5);
            greyImage->setPixel_(i, j, 2, 0.5);
        }
    }
    
    tempImage = ip_interpolate(src, greyImage, -1); //interpolate with grey image and alpha of -1
    
    return tempImage;
}


/*
* heatmap filter
*/

Image* ip_misc(Image* src)
{
    
//    Image* tempImage = new Image(src->getWidth(), src->getHeight());
//    Image* greyImage = ip_grey(src);//create grey version of source image
//    
//    //loop through image
//    for (int i = 0; i < src->getWidth(); ++i) {
//        for (int j = 0; j < src->getHeight(); ++j) {
//            double distance = greyImage->getPixel_(i, j, 0);//store value for distance from 0
//            if (distance < 85.0/255.0){ //compare to 1/3 of 255
//                tempImage->setPixel_(i, j, 0, 0);
//                tempImage->setPixel_(i, j, 1, 0);
//                tempImage->setPixel_(i, j, 2, 1-(distance/2));//set of blue channel to range of blue based on distance from 0
//            }
//            else if (distance >= 86.0/255.0 && distance <= 170.0/255.0){
//                tempImage->setPixel_(i, j, 0, 0);
//                tempImage->setPixel_(i, j, 1, 1-(distance/2));
//                tempImage->setPixel_(i, j, 2, 0);
//            }
//            else if (distance > 170.0/255.0){
//                tempImage->setPixel_(i, j, 0, 1-(distance/2));
//                tempImage->setPixel_(i, j, 1, 0);
//                tempImage->setPixel_(i, j, 2, 0);
//            }
//        }
//    }
//    
//    return tempImage;
    //gamma
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            
            tempImage->setPixel_(i, j, 0, pow(src->getPixel_(i, j, 0), 2));
            tempImage->setPixel_(i, j, 1, pow(src->getPixel_(i, j, 1), 2));
            tempImage->setPixel_(i, j, 2, pow(src->getPixel_(i, j, 2), 2));
        }
    }
    return tempImage;

}




/*
* round each pixel to the nearest value in the new number of bits
*/
Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight(), bitsPerChannel);

    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            tempImage->setPixel_(i, j, 0, src->getPixel_(i, j, 0));
            tempImage->setPixel_(i, j, 1, src->getPixel_(i, j, 1));
            tempImage->setPixel_(i, j, 2, src->getPixel_(i, j, 2));
        }
    }
    return tempImage;
}


/*
* dither each pixel to the nearest value in the new number of bits
* using a static 2x2 matrix
*/
Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* dither each pixel to the nearest value in the new number of bits
* using error diffusion
*/
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
    Image* tempImage = ip_quantize_simple(src, bitsPerChannel);
    //Image* tempImage = src;
    double errors[4];
    double errorRed = 0.0;
    double errorGreen = 0.0;
    double errorBlue = 0.0;
    
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            errorRed = src->getPixel_(i, j, 0) - tempImage->getPixel_(i,j, 0);
            tempImage->setPixel_(i+1, j, 0, (7.0/16.0*errorRed) + src->getPixel_(i+1,j,0));
            tempImage->setPixel_(i-1, j+1, 0, (3.0/16.0*errorRed) + src->getPixel_(i-1,j+1,0));
            tempImage->setPixel_(i, j+1, 0, (5.0/16.0*errorRed) + src->getPixel_(i,j+1,0));
            tempImage->setPixel_(i+1, j+1, 0, (1.0/16.0*errorRed) + src->getPixel_(i+1,j+1,0));
            
            errorGreen = src->getPixel_(i, j, 1) - tempImage->getPixel_(i,j, 1);
            tempImage->setPixel_(i+1, j, 1, (7.0/16.0*errorGreen) + src->getPixel_(i+1,j,1));
            tempImage->setPixel_(i-1, j+1, 1, (3.0/16.0*errorGreen) + src->getPixel_(i-1,j+1,1));
            tempImage->setPixel_(i, j+1, 1, (5.0/16.0*errorGreen) + src->getPixel_(i,j+1,1));
            tempImage->setPixel_(i+1, j+1, 1, (1.0/16.0*errorGreen) + src->getPixel_(i+1,j+1,1));
            
            errorBlue = src->getPixel_(i, j, 2) - tempImage->getPixel_(i,j, 2);
            tempImage->setPixel_(i+1, j, 2, (7.0/16.0*errorBlue) + src->getPixel_(i+1,j,2));
            tempImage->setPixel_(i-1, j+1, 2, (3.0/16.0*errorBlue) + src->getPixel_(i-1,j+1,2));
            tempImage->setPixel_(i, j+1, 2, (5.0/16.0*errorBlue) + src->getPixel_(i,j+1,2));
            tempImage->setPixel_(i+1, j+1, 2, (1.0/16.0*errorBlue) + src->getPixel_(i+1,j+1,2));
            
            tempImage = ip_quantize_simple(tempImage, bitsPerChannel);
        }
    }
    return tempImage;
}

/*
* nearest neighbor sample
*/
Pixel ip_resample_nearest(Image* src, double x, double y) {
	
	Pixel myPixel(0,0,0);
    double bigX = ((int)x+1.0);
    double smallX = ((int)x);
    double bigY = ((int)y+1.0);
    double smallY = ((int)y);
    
    int finalX = 0;
    int finalY = 0;
    
    if((bigX - x) < (x - smallX)){
        finalX = bigX;
    }
    else{
        finalX = smallX;
    }
    
    if((bigY - y) < (y - smallY)){
        finalY = bigY;
    }
    else{
        finalY = smallY;
    }
    
    myPixel = src->getPixel_(finalX, finalY);
    
	return myPixel;
}

/*
* bilinear sample
*/

Pixel ip_resample_bilinear(Image* src, double x, double y) {
	
	Pixel myPixel(0,0,0);
    double bigX = ((int)x+1.0);
    double smallX = ((int)x);
    double bigY = ((int)y+1.0);
    double smallY = ((int)y);
    
    double redX =((1- x - smallX) * src->getPixel_(smallX, y, 0)) + ((1-bigX-x)* src->getPixel_(bigX, y, 0));
    double redY = ((1- y - smallY) * src->getPixel_(x, smallY, 0)) + ((1-bigY-y)* src->getPixel_(x, bigY, 0));
    double greenX =((1- x - smallX) * src->getPixel_(smallX, y, 1)) + ((1-bigX-x)* src->getPixel_(bigX, y, 1));
    double greenY = ((1- y - smallY) * src->getPixel_(x, smallY, 1)) + ((1-bigY-y)* src->getPixel_(x, bigY, 1));
    double blueX =((1- x - smallX) * src->getPixel_(smallX, y, 2)) + ((1-bigX-x)* src->getPixel_(bigX, y, 2));
    double blueY = ((1- y - smallY) * src->getPixel_(x, smallY, 2)) + ((1-bigY-y)* src->getPixel_(x, bigY, 2));
    
    myPixel.setColor(0, (redX+redY)/2);
    myPixel.setColor(1, (greenX+greenY)/2);
    myPixel.setColor(2, (blueX+blueY)/2);
    
	return myPixel;
}

/*
* gausian sample
*/
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma)
{
	
Pixel myPixel(0,0,0);
//    
//    Image* tempImage = new Image(src->getWidth(), src->getHeight());
//    double filter[size*size];//create kernel with width and height = size
//    double distance[size*size];
//    double dAccumulator = 0.0;//accumulator value saves the total sum of all values in the kernel
//    double accumulator = 0.0;
//    
//    //loop through kernel
//    for (int i = 0; i<size; ++i) {
//        for (int j = 0; j<size; ++j) {
//            int a = (((int)(x)) + (size/2)-1+i);
//            int b = (((int)(y)) + (size/2)-1+j);
//            filter[j*size+i]= exp((-(((x-a)*(x-a)) + ((y-b)*(y-b))))/(2.0*(sigma*sigma))) * src->getPixel_(a, b, 0);
//            distance[j*size+i] = sqrt(((x-a)*(x-a)) + ((y-b)*(y-b)));
//            dAccumulator += distance[j*size+i];
//        }
//    }
//    //loop through kernel again
//    for (int i = 0; i<size; ++i) {
//        for (int j = 0; j<size; ++j) {
//            distance[j*size+i] /= dAccumulator;
//        }
//    }
//    
//    for (int i = 0; i<size; ++i) {
//        for (int j = 0; j<size; ++j) {
//            distance[j*size+i] /= accumulator;
//        }
//    }

return myPixel;
}

/*
* rotate image using one of three sampling techniques
*/
Image* ip_rotate (Image* src, double theta, int x, int y, int samplingMode, 
				  int gaussianFilterSize, double gaussianSigma)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    
    if(samplingMode == 0){
        for (int i = 0; i < src->getWidth(); ++i) {
            for (int j = 0; j < src->getHeight(); ++j) {
                Pixel resampledPixel = ip_resample_nearest(src, i, j);
                tempImage->setPixel_((((i-x)*(cos(theta))-((j-y)*sin(theta)))), (((j-y)*(cos(theta))+ ((i-x)*(sin(theta))))), resampledPixel);
            }
        }
    }
    else if(samplingMode == 1){
        for (int i = 0; i < src->getWidth(); ++i) {
            for (int j = 0; j < src->getHeight(); ++j) {
                Pixel resampledPixel = ip_resample_bilinear(src, i, j);
                tempImage->setPixel_((((i-x)*(cos(theta))-((j-y)*sin(theta)))), (((j-y)*(cos(theta))+ ((i-x)*(sin(theta))))), resampledPixel);
            }
        }
    }
    else if (samplingMode == 2){
        //go die
    }
    
    return tempImage;
}


/*
* change saturation
*/
Image* ip_saturate (Image* src, double alpha)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    Image* greyScaleImage = ip_grey(src);
    
    tempImage = ip_interpolate(src, greyScaleImage, alpha);
    
    return tempImage;
}


/*
* scale image using one of three sampling techniques
*/
Image* ip_scale (Image* src, double xFac, double yFac, int samplingMode, 
				 int gaussianFilterSize, double gaussianSigma)
{
	cerr << "This function is not implemented." << endl;
	return NULL;
}


/*
* threshold image
*/
Image* ip_threshold (Image* src, double cutoff)
{
    Image* tempImage = new Image(src->getWidth(), src->getHeight());
    
    //loop through image
    for (int i = 0; i < src->getWidth(); ++i) {
        for (int j = 0; j < src->getHeight(); ++j) {
            
            //channel compare to cutoff
            if (src->getPixel_(i, j, 0) <= cutoff){
                tempImage->setPixel(i, j, 0, 0.0);//set value to 0 if below cutoff
            }
            else{
                tempImage->setPixel(i, j, 0, 1.0);//set value to 1 if above cutoff
            }
            
            if (src->getPixel_(i, j, 1) <= cutoff){
                tempImage->setPixel(i, j, 1, 0.0);
            }
            else{
                tempImage->setPixel(i, j, 1, 1.0);
            }
            
            if (src->getPixel_(i, j, 2) <= cutoff){
                tempImage->setPixel(i, j, 2, 0.0);
            }
            else{
                tempImage->setPixel(i, j, 2, 1.0);
            }
        }
    }
    
    return tempImage;
}




