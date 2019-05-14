//g++ PeronaMalik.cpp -o PeronaMalik -Icfitsio -Lcfitsio -lcfitsio -lm -lnsl -O2 `pkg-config --cflags`  ###command line to compile
// ./PeronaMalik image_name filter_function_number(INT:1-5) parameter_k no_noise_image_name     ###command line to run

#include <iostream>
#include <fstream>
#include "fitsio.h"
#include <math.h>
#include <cmath>
#include <algorithm>

#define npad (4)

using namespace std;

int main(int argc, char *argv[]){

    fitsfile *fptr;
    fitsfile *fptr2;

    int status = 0;
    int status2 = 0;
    int ni,nj,bp,h=1;
    int ni2,nj2,bp2,h2=1;
    int T = 0;
    int bitpix, naxis;
    int bitpix2, naxis2;
    char* name = argv[1];
    int g = atoi(argv[2]);
    float k = atof(argv[3]);
    char* original = argv[4];
    long naxes[2] = {1,1}; 
    long naxes2[2] = {1,1}; 
    long fpixel[2] = {1,1};
    long fpixel2[2] = {1,1};
    double mse=1000;
    double mse_old=2000;
    double mae_old=2000000000000;
    double eps=1e-6;
    double dt=0.2;
    double min_val=1e10;
    double min_val2 = min_val;
    double max_val=-min_val;
    double max_val2 = max_val;
    double C1,C2;
    
    bool cond,cond2;
    ofstream outfile2 ("PM_mse_T_all.txt");
    fits_open_file(&fptr, name, READONLY, &status);
    fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);

    fits_open_file(&fptr2, original, READONLY, &status);
    fits_get_img_param(fptr2, 2, &bitpix2, &naxis2, naxes2, &status2);

    ni=naxes[1];
    nj=naxes[0];
    ni2=naxes2[1];
    nj2=naxes2[0];
    bp=bitpix;
    bp2=bitpix2;
    double* pixels;
    double* pixels2;
    double** u0;
    double** u;
    double** Orig;
    double** N;
    double** S;
    double** W;
    double** E;
    double** DuN;
    double** DuS;
    double** DuW;
    double** DuE;
    

    pixels = new double[nj];
    pixels2 = new double[nj2];
    u0 = new double*[ni];
    u = new double*[ni];
    Orig = new double*[ni2];
    N = new double*[ni];
    S = new double*[ni];
    W = new double*[ni];
    E = new double*[ni];
    DuN = new double*[ni];
    DuS = new double*[ni];
    DuW = new double*[ni];
    DuE = new double*[ni];
    for (int i=0; i<ni; i++)
    {
	u0[i]=new double[nj+npad];
        u[i]=new double[nj+npad];
        Orig[i] = new double[nj2+npad];
        N[i]=new double[nj+npad];
	S[i]=new double[nj+npad];
        W[i]=new double[nj+npad];
        E[i]=new double[nj+npad];
	DuN[i]=new double[nj+npad];
        DuS[i]=new double[nj+npad];
        DuW[i]=new double[nj+npad];
        DuE[i]=new double[nj+npad];
    }
    
    
    cout << " ni = " << ni 
         << " nj = " << nj 
         << " bp = " << bp 
         << " name = " << name
         << " naxis = " << naxis       
         << " status = " << status 
         << " k = " << k  
         << " g = " << g            
         << endl;
    cout << " ni = " << ni2 
         << " nj = " << nj2 
         << " bp = " << bp2 
         << " name = " << original
         << " naxis = " << naxis2       
          << " status = " << status2 
          << " k = " << k  
          << " g = " << g            
          << endl;
        
    //filling up M with pixels read from the input_image with fits_read_pix
    for(int i=0;i<ni;i++)
    {

        fpixel[0]=1;
        fpixel[1]=ni-i;
        fpixel2[0]=1;
        fpixel2[1]=ni-i;

        fits_read_pix(fptr, TDOUBLE, fpixel, nj, NULL, pixels, NULL, &status);
        fits_read_pix(fptr2, TDOUBLE, fpixel2, nj2, NULL, pixels2, NULL, &status2);

	//finding min and max values in input image
        for(int j=0; j<nj; j++)
        {
            u0[i][j] = pixels[j];
            Orig[i][j] = pixels2[j];
            max_val = max(pixels[j], max_val); 
            min_val = min(pixels[j], min_val);
            max_val2 = max(pixels2[j], max_val2); 
            min_val2 = min(pixels2[j], min_val2); 
        }
    
    }

    //normalization of the input image

    for(int i=0; i<ni; ++i) 
    {
        for(int j=0; j<nj; ++j) 
        {
            u0[i][j] = (u0[i][j]-min_val)/(max_val-min_val);
            Orig[i][j] = (Orig[i][j]-min_val2)/(max_val2-min_val2);
        }
    } 

    //Beginning of T cycle, for each cycle the Perona-Malik algorithm is applied here
        
    while (T<500 and mse_old > mse and (mse_old - mse)/mse_old > 1e-10)
    {   
        cout << "T" << T << endl;
        mse_old = mse;        
        mse=0; 
        mae_old = mae;
        mae=0;       
        T++;

        for(int i=0;i<ni;i++)
        {
            for(int j=0; j<nj; j++)
            {   
                if(i != 0)
                {
                    DuN[i][j] = u0[i-1][j]-u0[i][j];
                }
                else
                {
                    DuN[i][j] = 0;
                }
                if(i != ni-1)
                {
                    DuS[i][j] = u0[i+1][j]-u0[i][j];
                }
                else
                {
                    DuS[i][j] = 0;
                }
                if(j != 0)
                {
                    DuW[i][j] = u0[i][j-1] - u0[i][j];
                }
                else
                {
                    DuW[i][j] = 0;
                }
                if(j != nj-1)
                {
                    DuE[i][j] = u0[i][j+1] - u0[i][j];
                }
                else
                {
                    DuE[i][j] = 0;
                }              
            }
        }

        for(int i=0;i<ni;i++)
        {
            for(int j=0; j<nj; j++)
            {
				//here are defined the filtering functions
                if(g == 1)
                {
                    N[i][j] = 1./(1.+(DuN[i][j]/k)*(DuN[i][j]/k));
                    S[i][j] = 1./(1.+(DuS[i][j]/k)*(DuS[i][j]/k));
                    W[i][j] = 1./(1.+(DuW[i][j]/k)*(DuW[i][j]/k));
                    E[i][j] = 1./(1.+(DuE[i][j]/k)*(DuE[i][j]/k));
                }
                if(g == 2)
                {
                    N[i][j] = exp(-(DuN[i][j]/k)*(DuN[i][j]/k));
                    S[i][j] = exp(-(DuS[i][j]/k)*(DuS[i][j]/k));
                    W[i][j] = exp(-(DuW[i][j]/k)*(DuW[i][j]/k));
                    E[i][j] = exp(-(DuE[i][j]/k)*(DuE[i][j]/k));
                }
                if(g == 3)
                {
                    if (fabs(DuN[i][j]) <= k*sqrt(2))
                    {
                        N[i][j] = 1./2*pow(1.-pow(DuN[i][j]/(k*sqrt(2)),2),2);
                    }                    
                    else
                    {
                        N[i][j] = 0;
                    } 
                    if (fabs(DuS[i][j]) <= k*sqrt(2))
                    {
                        S[i][j] = 1./2*pow(1.-pow(DuS[i][j]/(k*sqrt(2)),2),2);
                    }                    
                    else
                    {
                        S[i][j] = 0;
                    } 
                    if (abs(DuW[i][j]) <= k*sqrt(2))
                    {
                        W[i][j] = 1./2*pow(1.-pow(DuW[i][j]/(k*sqrt(2)),2),2);
                    }                    
                    else
                    {
                        W[i][j] = 0;
                    } 
                    if (abs(DuE[i][j]) <= k*sqrt(2))
                    {
                        E[i][j] = 1./2*pow(1.-pow(DuE[i][j]/(k*sqrt(2)),2),2);
                    }                    
                    else
                    {
                        E[i][j] = 0;
                    } 
                }
                if (g == 4)
                {
                    N[i][j] = 1./(1.+pow((fabs(DuN[i][j])/k),2.-2/(1+pow(fabs(DuN[i][j])/k,2))));
                    S[i][j] = 1./(1.+pow((fabs(DuS[i][j])/k),2.-2/(1+pow(fabs(DuS[i][j])/k,2))));
                    W[i][j] = 1./(1.+pow((fabs(DuW[i][j])/k),2.-2/(1+pow(fabs(DuW[i][j])/k,2))));
                    E[i][j] = 1./(1.+pow((fabs(DuE[i][j])/k),2.-2/(1+pow(fabs(DuE[i][j])/k,2))));                    
                }
                if (g == 5)
                {
                    if (abs(DuN[i][j]) != 0)
                    {
                        N[i][j] = 1- exp(-3.31488*pow(k/DuN[i][j],8));
                    }                    
                    else
                    {
                        N[i][j] = 1;
                    } 
                    if (abs(DuS[i][j]) != 0)
                    {
                        S[i][j] = 1- exp(-3.31488*pow(k/DuS[i][j],8));
                    }                    
                    else
                    {
                        S[i][j] = 1;
                    } 
                    if (abs(DuW[i][j]) != 0)
                    {
                        W[i][j] = 1- exp(-3.31488*pow(k/DuW[i][j],8));
                    }                    
                    else
                    {
                        W[i][j] = 1;
                    } 
                    if (abs(DuE[i][j]) != 0)
                    {
                        E[i][j] = 1- exp(-3.31488*pow(k/DuE[i][j],8));
                    }                    
                    else
                    {
                        E[i][j] = 1;
                    } 
                }
            }
        }
        
        for(int i=0;i<ni;i++)
        {
            for(int j=0; j<nj; j++)
            {           
                u[i][j] = u0[i][j] + dt * (N[i][j]*DuN[i][j] + S[i][j]*DuS[i][j] + W[i][j]*DuW[i][j] + E[i][j]*DuE[i][j]);
            }                                    
        }

        u0 = u;
        
		//here I measure the MSE
        for(int i=0;i<ni;i++)
        {
           for(int j=0; j<nj; j++)
            {
                mse += pow(u[i][j]-Orig[i][j],2);                
            }
        }
        mse /= (ni*nj);
    }
    //end of T cycle and denormalyze back the denoised image using max_val and min_val from input
    for(int i=0; i<ni; ++i) 
    {
        for(int j=0; j<nj; ++j) 
        {
            u0[i][j] = u0[i][j]*(max_val-min_val)+min_val;
        }
    }
      
    //writing the fits file for the output with name DenoisedPM.fits
    double *array = new double[nj];
    fitsfile *outfptr;
    string str = "DenoisedPM.fits";
    name = (char*)str.c_str();
    fits_create_file(&outfptr, name, &status);
    fits_create_img(outfptr, bitpix, naxis, naxes, &status);
   
    for(int i=0;i<ni;i++)
    {
        for(int j=0;j<nj;j++) 
        {
            array[j]=u0[i][j];
        }
        
	    fpixel[0]=1;
	    fpixel[1]=ni-i;
        fits_write_pix(outfptr, TDOUBLE, fpixel, nj, array, &status);
    }
    
    fits_close_file(outfptr,  &status);    
    return 0;
}
