#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>
#include <iostream>  //  C++ stream IO
#include <sstream>
#include <math.h>
#include <string>
#include <iomanip>   //  to format the output
#include <unistd.h>
using namespace std;
#define eulerfile_line 500000
#define DF_MAX -10000
#define DF_MIN -5000
#define DF_STEP 1000	//no longer used


using namespace std;

#include "cfpix.hpp"       /* complex image handler with FFT */
#include "slicelib.hpp"    /* misc. routines for multislice */
#include "floatTIFF.hpp"   /* file I/O routines in TIFF format */
//#include "fftw3.h"
const int MCPIX=   0;   /* coherent image mode */
const int MPCPIX=  1;   /* partial coherent mode */
const int MDIFFP=  2;   /* diffraction mode */

// define more parameter keys not in slicelib.h
// make symbolic offsets to merge with slicelib.h def's 
const int
    pB = 260,
    p51 =236,
    p52 =237,

    p31 =251,
    p32 =252,
    p33 =253,
    p34 =254,
    p35 =255,
    p36 =256,
    p37 =257,
    p38 =258,
    p39 =259;

/*  define subroutines at end of file */
void invert2D( cfpix &pix, int nx, int ny );
void tcross( cfpix &pixin, cfpix &pixo,
        int nx, int ny, float *p );

void display_usage(){
	cout<<"This program produce multislice tif before imaging."<<endl;
	cout<<"-a <apix>"<<endl;
	cout<<"-i <Input list>"<<endl;
	cout<<"-f -l <first to last>"<<endl;
}
static const char *optString = "a:i:f:l:h?";

int main(int argc, char *argv[])
{
	int opt = 0;
	string tmp;
    string filein, fileout,temp0,nname;
    const char version[] = "21-jun-2015 (ejk)";
//	floatTIFF myFile;
    int ix, iy, nx, ny, ixmid, iymid, npix, ns,
        itens, mode, nsum, NPARAM,number_of_diff_df_steps,i;
    int lcenter=0, lapert=0;

    float kx,ky, kx2,ky2, k2, k2max, v0, wavlen, scale, pixc,
        ax, by, cz, rx, ry, pi, rmin, rmax, aimin, aimax,
        Cs3, Cs5, df,  alpha0, ddf, objlx, objly,
        tr, ti, wr, wi, dfa2, dfa2phi, dfa3, dfa3phi,df_max,df_min,df_step;
    float *param;
    float  **pixr, **pixi;
	float amplitude=1.0,apix=1.0;
    double sum, chi, chi1, chi2, chi3, phi, clog;
	std::string filename[2];
//	std::stringstream stream0;
    cfpix cpix, cpix2;     //  complex floating point image
	string list_filein;
	std::string str[eulerfile_line+2];
	int jmin,jmax;
//	cout<<"Enter jmin jmax"<<endl;
//	cin>>jmin>>jmax;

/*  echo version date */

//    cout << "image version dated " << version << endl;
 //   cout << "Copyright (C) 1998-2014 Earl J. Kirkland" << endl;
 //   cout <<  "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
 //       << " under the GNU general public license\n" << endl;

 //   cout <<  "calculate TEM images with defocus, using FFTW\n"<< endl;
	opt = getopt(argc,argv,optString);
	while(opt!=-1) {
		switch(opt) {
			case 'a':
				tmp=optarg;
				apix=atof(tmp.c_str());
				break;
			case 'i':
				list_filein=optarg;
				break;
			case 'f':
				tmp=optarg;
				jmin=atoi(tmp.c_str());
				break;
			case 'l':
				tmp=optarg;
				jmax=atoi(tmp.c_str());
				break;
			case '?':
                display_usage();
				return 0;
                break;
			case 'h':
				display_usage();
				return 0;
                break;
		}
		opt = getopt( argc, argv, optString );
	}
    pi = (float) (4.0 * atan( 1.0 ));
	if(jmax<jmin)
	{
		cout<<"max smaller than min"<<endl;
		return 0;
	}
/*  get input file name */
	srand((unsigned)time(0));

//	cout << "Name of list of autosliced files" << endl;
 //   cin >> list_filein;
	std::ifstream infile(list_filein.c_str());
	int countN=0;
	while(!infile.eof()){
		getline(infile, str[countN],'\n');
		countN++;
    }
	if(countN-1<jmax){
		cout<<"jmax too large"<<endl;
		return 0;
	}
//	cout<<"Enter amplitude"<<endl;
//	cin>>amplitude;
//	amplitude=1.0;
//	cout<<"Enter apix"<<endl;
//	cin>>apix;
	int j=0,just_count=0;
	df_min=DF_MIN;
	df_max=DF_MAX;
	df_step=DF_STEP;
	if(fabs(df_min)>fabs(df_max)){
		printf("min over max,exit\n");
		return 0;
	}

//	number_of_diff_df_steps=floor((fabs(df_max)-fabs(df_min))/fabs(df_step))+1;
	number_of_diff_df_steps=1;
	int just_total;
	just_total=(jmax-jmin+1)*number_of_diff_df_steps;
	mode = MCPIX;
	ofstream fpout;
//	fpout.open("ctfparm.txt",ios::app);

for(j=jmin-1;j<jmax;j++){
    floatTIFF myFile;
	filename[0]=str[j];
    filein = filename[0];
	filename[1]=filein;
/*  open input file and check sizes etc. */

    if( myFile.tFloatTest( filein.c_str() ) != 1 ) {
        cout << "Cannot open input file " << filein << endl;
        exit( 0 );
    }

	for(i=0;i<number_of_diff_df_steps;i++){
//		cfpix cpix;  
		Cs3 = 0.0F;
		Cs5 = 0.0F;
		df = 0.0F;
		k2max = 0.0F;
		objlx = 0.0F;
		objly = 0.0F;
   
        Cs3=-2.7;
		Cs5=0;
        Cs3 *= 1.0e7F;   /*  convert to Angstroms */
        Cs5 *= 1.0e7F;
		float rnd=1*floor(((df_max-df_min)/1)*rand()/double(RAND_MAX));
        df=df_min+rnd;
		ostringstream oss; 
		oss << -df; 
		string temp0 = oss.str(); 
		oss.str("");
		fileout = "Image_"+filename[1]+"_df_"+temp0+".tif";
		nname="Image_"+filename[1]+"_df_"+temp0;
        k2max=1000;

		dfa2=dfa2phi=0;
        dfa2phi = dfa2phi * pi /180.0F;
        dfa3=dfa3phi=0;
        dfa3phi = dfa3phi * pi /180.0F;
        objlx =objly=0;
       

/*  get diffraction pattern parameters if required  */

    
    /* ---- read in specimen parameters and multislice result ---- */

//    time = cputim();    /* get CPU time for comparison */

    NPARAM = myFile.maxParam();
    param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
    for( ix=0; ix<NPARAM; ix++) param[ix] = 0.0F;

    if( myFile.read( filein.c_str() ) != 1 ) {
        cout << "Cannot open input file " << filein  << endl;
        exit( 0 );
    }
    nx = myFile.nx();
    nx = nx/2;      /* should be complex */
    ny = myFile.ny();
    npix = myFile.getnpix();
    ixmid = nx/2;
    iymid = ny/2;

    if( npix != 2 ) {
        cout << "Input file " << filein << " must be complex, can't continue."
              << endl;
        exit( 0 );
    }

    pixr = (float**) malloc2D( 2*nx, ny, sizeof(float), "pixr" );
    pixi = pixr + nx;

    //  crude port of old code using param[]
    for( ix=0; ix<NPARAM; ix++) param[ix] = myFile.getParam(ix);

    ax = param[pDX] * ((float) nx);
    by = param[pDY] * ((float) ny);
    if( param[pDY] <= 0.0F ) by = param[pDX] * ((float) ny);
    cz = param[pC];
    rx  = 1.0F / ax;
    ry  = 1.0F / by;

    v0 = param[pENERGY];
    wavlen = (float) wavelength( v0 );
//    cout << "Starting pix energy = " << v0 << " keV"  << endl;

    rmax  = param[pRMAX];
    aimax = param[pIMAX];
    rmin  = param[pRMIN];
    aimin = param[pIMIN];
//    cout << "Starting pix range " << rmin << "  " << rmax << " real\n"
//            "                   " << aimin << "  " << aimax << " imag." << endl;

    param[pDEFOCUS] = df;
    param[pASTIG] = 0.0F;
    param[pTHETA] = 0.0F;
    param[pOAPERT] = k2max * 0.001F;
    param[pCS] = Cs3;
    param[pCS5] = Cs5;
    param[pWAVEL] = wavlen;

    k2max = k2max*0.001F/wavlen;
    k2max = k2max * k2max;

    objlx = objlx * 0.001F / wavlen;   /* convert to spatial freq. */
    objly = objly * 0.001F / wavlen;

    cpix.resize( nx, ny );
    cpix.init( 1 );   //  fast init but slow execution is better here

    /*------  copy to FFTW style array -------
        not very elegant but... */
    for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            cpix.re(ix,iy) = myFile(ix,iy);    // real
            cpix.im(ix,iy) = myFile(ix+nx,iy); // imag
    }

    cpix.fft();   //  every option uses FT so do it here

/*------  coherent real space image mode ---------------------

    Convolute multislice result with objective lens
    aberration function (coherent case)
    with shifted objective lens for dark field

   NOTE zero freq. is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for FFT
*/

    if( mode == MCPIX ) {
//        cout << "calculate coherent image...." << endl;
        chi1  = pi * wavlen;
        chi2  = 0.5 * Cs3 * wavlen * wavlen;
        chi3  = Cs5 * wavlen*wavlen*wavlen*wavlen /3.0;

        ns = 0;
        for( ix=0; ix<nx; ix++) {
            kx = (float) ix;
            if( ix > ixmid ) kx = (float) (ix-nx);
            kx = kx*rx - objlx;
            kx2 = kx*kx;
            for( iy=0; iy<ny; iy++) {
                ky = (float) iy;
                if( iy > iymid ) ky = (float)(iy-ny);
                ky = ky*ry - objly;
                k2 = ky*ky + kx2;
                if ( k2 < k2max ) {
                    ns++;
                    phi = atan2( ky, kx );
                    chi = chi1*k2* ( (chi2 + chi3*k2)*k2 - df 
                        + dfa2*sin( 2.0*(phi-dfa2phi) ) 
                        + 2.0F*dfa3*wavlen*sqrt(k2)*
                        sin( 3.0*(phi-dfa3phi) )/3.0 );
                    tr = (float) amplitude*cos( chi );			//	cos( chi );
                    ti = (float) -1*amplitude*sin( chi );		//	-sin( chi );
                    wr = cpix.re(ix,iy);
                    wi = cpix.im(ix,iy);
                    cpix.re(ix,iy) = tr*wr - ti*wi;  /* real */
                    cpix.im(ix,iy) = tr*wi + ti*wr;  /* imag */
                 } else {
                    cpix.re(ix,iy) = 0.0F;  /* real */
                    cpix.im(ix,iy) = 0.0F;  /* imag */
                }
            }  /* end for( ix.... */
        }  /* end for(iy... */

        cpix.ifft();
//        cout << "There were " << ns << " pixels inside the obj. apert." << endl;

        /*   copy back to old style array for output */
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            pixr[ix][iy] = cpix.re(ix,iy);
            pixi[ix][iy] = cpix.im(ix,iy);
        }

    } 

/*  Output results and find min and max to echo */

//    cout << "output image" << endl;
    scale = 1.0F / ( ((float)nx) * ((float)ny) );
    
    myFile.resize( nx, ny );
    myFile.setnpix( 1 );
    sum = 0.0;
    nsum = 0;
    for( iy=0; iy<ny; iy++) {
        for( ix=0; ix<nx; ix++) {
            tr = pixr[ix][iy];
            ti = pixi[ix][iy];
            if( mode == MPCPIX ) pixc = tr;
            else  pixc = tr*tr + ti*ti;
            if ( (mode == MDIFFP) &&  (itens == 1)) {
                if( pixc > 1.e-30F)  pixc = (float) log( (double) pixc );
                else pixc = -30.0F;
            } else if ( (mode == MDIFFP) && (itens == 2)) {
                pixc = (float) log( 1.0 + clog*sqrt((double)pixc) );
            } 

            if( (ix == 0) && (iy == 0) ) {
                rmin = pixc;
                rmax = rmin;
            } else if( (ix != ixmid) && (iy != iymid) ) {
                if( pixc < rmin ) rmin = pixc;
                if( pixc > rmax ) rmax = pixc;
            }
            if( (ix>(3*nx)/8) && (ix<(5*nx)/8) &&
                (iy>(3*ny)/8) && (iy<(5*ny)/8) ) {
                sum = sum + pixc;
                nsum += 1;
            }
            myFile(ix,iy) = pixc;   //pixr[ix][iy] = pixc;
        }  /* end for ix... */
    } /* end for iy... */

    param[pRMAX] = rmax;
    param[pIMAX] = 0.0F;
    param[pRMIN] = rmin;
    param[pIMIN] = 0.0F;
    if ( mode == MDIFFP ) {
        param[pRMIN] = (float) (0.05*rmin + 0.95*sum/nsum);
        param[pDX] = 1.0F / (nx*param[pDX]);
        param[pDY] = 1.0F / (ny*param[pDY]);
    }

    //  crude port of old code using param[]
    for(ix=0; ix<NPARAM; ix++) myFile.setParam(ix, param[ix] );

    myFile.write( fileout.c_str(), rmin, rmax, 0.0F, 0.0F, param[pDX], param[pDY] );
//	fpout<<nname<<"\t"<<(df/10000)<<",10,"<<amplitude<<",0,0,0,0,0,300,2.7,"<<apix<<",(null)"<<endl;
//    cout << "Pix range " << param[pRMIN] << " to " << rmax << endl;

//    time = cputim() - time;
//	cout << "Elapsed time = " << time << " sec." << endl;
	just_count++;
//	cpix.~cfpix();
//	myFile.~floatTIFF();
//	fftwf_destroy_plan(cpix.planTi);
//	fftwf_destroy_plan(cpix.planTf);
//	fftwf_cleanup();
	if(just_count%100==0)
	printf("finish iteration %d of %d\n\n",just_count,just_total);
}


}
    return EXIT_SUCCESS;

} /* end main() */

/*------------------- tcross() ------------------------------*/
/*  this version uses cfpix style arrays

 Subroutines to do exact nonlinear partial coherence transfer as in
      [1]  M.A. O'Keefe, 37th EMSA (1979) p.556
      [2]  K. Ishizuka, Ultramicroscopy 5 (1980) p.55-65

  both input and output are in Fourier space with zero frequency
  in the center (NX/2 , NY/2)

 perform 2D weighted convolution with transmission cross coefficient
  TXCOEF to calculate partially coherent CTEM images
    NOTE: this version uses Friedels law

  pixin  : complex input pix array
             range used = P(17) = aperture
  pixo   : complex output pix array
             range output= 2.0*P(17) = 2 x aperture
  nx, ny : actual image size to work with
  p[43]  : parameter array ( dx,dy, defocus etc.) explained in TXCOEF

  started 14-NOV-85 on CIPRES2 Earl J. Kirkland
      started by modifying XFERFN.FTN (from RSX CIPRES)
  changed TXCOEF to be more general (i.e. include alignment and leading
      factors so that it is complete - little or no speed lose)
     28-NOV-85 EJK
  fixed typo (DT3 to D3T) and sign (-D4) in TXCOEF 17-jun-1986 ejk
  changed sign convention to be consistent with forward propagation
    of electrons (see Self et.al. UM 11 (1983) p.35-52 
      - only TXCOEF effected   EJK 2-JUL-86
  converted to C 8-nov-1995 ejk
  fixed sign error 23-jan-1998 ejk
  convert to FFTW style arrays 8-may-2010 ejk
  fix order of psi/psi* indexes to fix flipped image error
      17-jan-2011 ejk
  convert to cfpix style arrays 3-nov-2012 ejk
  remove extra factor pi in p[p52] in tcross() 21-jun-2015 ejk

*/

void tcross( cfpix &pixin, cfpix &pixo,
        int nx, int ny, float *p )
{
    int j, ix, iy, ixmid, iymid, ixmin, ixmax, iymin, iymax,
        ix2min, ix2max, iy0, iy2min, iy2max, ix1, ix2, iy1, iy2,
        *ixi, *ixf, ixi0, ixf0;
    float scale, k2p, k2pp, rx, ry, txcoefr, txcoefi,
        xr, xi, yr, yi, zr, zi, PI, *kxp, *kyp, *kxp2, *kyp2;

    void txcoef( float kxp, float kyp, float kxpp, float kypp,
         float p[], float *txcoefr, float *txcoefi );


/*  get scratch arrays and init params */

    ixi = (int*) malloc1D( nx, sizeof(int), "ixi" );
    ixf = (int*) malloc1D( nx, sizeof(int), "ixf" );
    kxp = (float*) malloc1D( nx, sizeof(float), "kxp" );
    kyp = (float*) malloc1D( ny, sizeof(float), "kyp" );
    kxp2 = (float*) malloc1D( nx, sizeof(float), "kxp" );
    kyp2 = (float*) malloc1D( ny, sizeof(float), "kyp2" );

    PI = (float) (4.0 * atan( 1.0 ));
    xr = p[pOAPERT]/p[pWAVEL];
    p[p31] = xr*xr;
    xi = p[pWAVEL];
    p[p32] = 0.5F*PI* p[pCS]* xi*xi*xi;
    p[p33] = PI* p[pWAVEL] * p[pDEFOCUS];

    p[p51] = PI * p[pCS5]* xi*xi*xi*xi*xi/3.0F;
    p[p52] = p[pCS5]* xi*xi*xi*xi;

    xr = PI*p[pCAPERT]* p[pDDF];
    p[p34] = xr*xr;
    p[p35] = p[pCS]* p[pWAVEL]* p[pWAVEL];
    xr = PI*p[pCAPERT];
    p[p36] = xr*xr;
    xr =  PI* p[pWAVEL]* p[pDDF];
    p[p37] = xr*xr/4.0F;
    xr = PI*PI*p[pCAPERT]*p[pCAPERT]*p[pDDF] ;
    p[p38] = xr*xr;
    xr =  PI*p[pCAPERT]*p[pDDF] ;
    p[p39] = PI*( xr*xr )*p[pWAVEL];

    /* initialize  */

    rx = p[pDX]* nx;
    ry = p[pDY]* ny;
    if( ry <= 0.0F) ry = rx;
    rx = 1.0F/rx;
    ry = 1.0F/ry;

    scale=1.0F/( ((float)nx) * ((float)ny) );

    ixmid = nx/2;
    iymid = ny/2;

    /* find range of convolution */

    j = (int) ( p[pOAPERT] / (rx*p[pWAVEL]) +1.0F );
    ixmin  = ixmid - j;
    ixmax  = ixmid + j;
    ix2min = ixmid - 2*j;
    ix2max = ixmid + 2*j;
    j = (int) ( p[pOAPERT] / (ry*p[pWAVEL]) +1.0F );
    iymin  = iymid - j;
    iymax  = iymid + j;
    iy2min = iymid - 2*j;
    iy2max = iymid + 2*j;
    
    if( ixmin < 0 ) ixmin = 0;
    if( ixmax > (nx-1) ) ixmax = (nx-1);
    if( ix2min < 0 ) ix2min = 0;
    if( ix2max > (nx-1) ) ix2max = (nx-1);
    if( iymin < 0 ) iymin = 0;
    if( iymax > (ny-1) ) iymax = (ny-1);
    if( iy2min < 0 ) iy2min = 0;
    if( iy2max > (ny-1) ) iy2max = (ny-1);

    for ( ix=ix2min; ix<=ix2max; ix++) {
       ix1= ixmin - ix + ixmid;
       if( ixmin > ix1 ) ixi[ix] = ixmin;  else  ixi[ix] = ix1;
       ix1= ixmax - ix + ixmid;
       if( ixmax < ix1 ) ixf[ix] = ixmax;   else  ixf[ix] = ix1;
    }

    for( ix=ix2min; ix<=ix2max; ix++) {
        kxp[ix] = (ix-ixmid) * rx;
        kxp2[ix]= kxp[ix] * kxp[ix];
    }
    for( iy=iy2min; iy<=iy2max; iy++) {
        kyp[iy] = (iy-iymid) * ry;
        kyp2[iy]= kyp[iy] * kyp[iy];
    }

    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
        pixo.re(ix,iy) = pixo.im(ix,iy) = 0.0F;

/*  do convolution in bottom half plane 
    the origin is at (ixmid,iymid)
    (ix,iy)   = output point = k
    (ix1,iy1) = input point from psi* = kp
    (ix2,iy2) = input point from psi  = kp + k
*/

    for( iy=iy2min; iy<=iymid; iy++ )
        for( iy1=iymin; iy1<=iymax; iy1++ ) {
            iy2 = iy1 + (iy - iymid);
            if( (iy2>=iymin) && (iy2<=iymax)) {
                for( ix=ix2min; ix<=ix2max; ix++ ) {
                    ixi0 = ixi[ix];
                    ixf0 = ixf[ix];
                    for( ix1=ixi0; ix1<=ixf0; ix1++ ) {
                        ix2= ix1 + (ix - ixmid);
                        k2p = kxp2[ix1] + kyp2[iy1];
                        k2pp= kxp2[ix2] + kyp2[iy2];
                        if ( (k2p<=p[p31]) && (k2pp<=p[p31]) ) {
                            txcoef( kxp[ix1],kyp[iy1],kxp[ix2],kyp[iy2],
                                p, &txcoefr, &txcoefi);
                            xr = pixin.re(ix1,iy1);  /* rpixin[ix1][iy1]; */
                            xi = -pixin.im(ix1,iy1); /* -ipixin[ix1][iy1]; */
                            yr = pixin.re(ix2,iy2);  /* rpixin[ix2][iy2]; */
                            yi = pixin.im(ix2,iy2);  /* ipixin[ix2][iy2]; */
                            zr = xr*yr - xi*yi;
                            zi = xr*yi + xi*yr;
                            pixo.re(ix,iy) += zr*txcoefr - zi*txcoefi;
                            pixo.im(ix,iy) += zr*txcoefi + zi*txcoefr;
                        }  /* end if( ( k2p... */
                    }  /* end for( ix1=... */
                }  /* end for( ix=... */
            }  /* end if( (iy2... */
        }  /* end for iy1... */

/*  scale result */

    for( ix= ix2min; ix<=ix2max; ix++) 
    for( iy= iy2min; iy<=iymid; iy++) {
        pixo.re(ix,iy) *= scale;
        pixo.im(ix,iy) *= scale;
    }

/*  Invoke Friedel's law to get top half plane  */

    iy0 = iymid+1;
    for( iy=iy0; iy<=iy2max; iy++) {
        iy1= iymid - (iy-iymid);
        for( ix=ix2min; ix<=ix2max; ix++) {
            ix1= ixmid - (ix-ixmid);
            pixo.re(ix,iy) =  pixo.re(ix1,iy1); /* rpixo[ix1][iy1]; */
            pixo.im(ix,iy) = -pixo.im(ix1,iy1); /* -ipixo[ix1][iy1]; */
        }
       pixo.re(0,iy) = pixo.im(0,iy) = 0.0F;  /* ??? */
    }

    free( ixi );
    free( ixf );
    free( kxp );
    free( kyp );
    free( kxp2 );
    free( kyp2 );

    return;
}     /*   end tcross() */

/*------------------------ txcoef -------------------------*/
/*
   The cross spectral transfer function as in:
    [1] M.A. O'Keefe, 37th EMSA (1979) p.556 EJK
    [2] K. Ishizuka, Ultramic. 5 (1980) p.55-65.

   switched to my derivation 5-DEC-83 EJK
   changed sign convention to be consistent with forward propagation
     of electrons (see Self et.al. UM 11 (1983) p.35-52 
       EJK 2-JUL-86
   converted to C  6-nov-1995 ejk
   fixed sign error 23-jan-1998 ejk
   fix sign to fix flipped image error
      17-jan-2011 ejk

  the following are the array index definitions
  (the order is historical - don't ask why! )

   P(11) = defocus
   P(12) = defocus astigmatism
   P(13) = angle of astigmatism
   P(14) = dx in A pixel dimension
   P(15) = dy in A pixel dimension
   P(17) = aperture in radians
   P(18) = Cs3 in A
   P(19) = wavelength in A
   P(21) = illumination semiangle in rad
   P(22) = defocus spread in A
   P(24) = Debye Waller temp factor/4 in A
   P(31) = P(17)/P(19) **2 maximum k^2 value
 
   P(32) = PI*P(18)* (P(19)**3) /2    ; coherent part
   P(33) = PI* P(19) * P(11)

   p(pCS5) = Cs5 in A         ; add 5th order 21-dec-2010 ejk
   p(p51) = PI*P(pCS5)* (P(19)**5) /3.0 
 
   P(34) = (PI* P(21)* P(22))**2
   P(35) = P(18)* P(19)* P(19)
   P(36) = (PI* P(21)) **2
   P(37) = (PI* P(19)* P(22))**2 /4.
   P(38) = PI^4 * P(21)^4 *P(22)^2
   P(39) = PI^3 * P(22)^2 * P(21)^2 *P(19)
 
   P(p52) = p(pCS5) * p(pWAVEL)^4 
     
   kp  = kprime      = argument of psi*
   kpp = kprime + k  = argument of psi
 
   txcoefr, txcoefi = returned real and imag. parts of function
        (i.e. C can't return a complex value )
*/

void txcoef( float kxp, float kyp, float kxpp, float kypp,
         float p[], float *txcoefr, float *txcoefi )
{
    double kx, ky, k2 ,kp2, kpp2, chip, chipp,
        d1, d2, d3, d3t, d4, d4t, vx, vy, w, u, xr;

    kp2  = kxp*kxp + kyp*kyp;
    kpp2 = kxpp*kxpp + kypp*kypp;

    if( ( kp2 > p[p31] ) || ( kpp2 > p[p31] ) ) {
      *txcoefr = *txcoefi = 0.0F;
      return;
    }

    kx = kxpp - kxp;   /* output freq. k */
    ky = kypp - kyp;
    k2 = kx*kx + ky*ky;

    /*  v = Wc1/(2*pi*lambda)
        w = Wc2/(-pi*lambda)
        u = 1 + pi^2 * beta^2 * delta0^2 k^2
    */
/* -- add Cs5 below 21-dec-2010 ejk */
    vx = p[p52]*( kp2*kp2*kxp - kpp2*kpp2*kxpp )
        + p[p35]*( kp2*kxp - kpp2*kxpp ) + p[pDEFOCUS]*kx;
    vy = p[p52]*( kp2*kp2*kyp - kpp2*kpp2*kypp )
        + p[p35]*( kp2*kyp - kpp2*kypp ) + p[pDEFOCUS]*ky;
    w  = kp2 - kpp2;
    u  = 1.0 + p[p34]*k2;

    d1  = p[p36] * (vx*vx+vy*vy);
    d2  = p[p37] *w*w /u;
    d3t = vx*kx + vy*ky;
    d3  = p[p38]*d3t*d3t/u;
    d4t = p[p39]*w/u;
    d4  = d4t*d3t;

    /* ---- coherent chi  --------------------
         add 5th order Cs5 21-dec-2010 ejk */
    chip  = ( (p[p51]*kp2  + p[p32] )*kp2  - p[p33] ) *kp2;
    chipp = ( (p[p51]*kpp2 + p[p32] )*kpp2 - p[p33] ) *kpp2;

    //chip = chip - chipp - d4;  fix sign 21-jun-2015 ejk
    chip = chip - chipp + d4;
    xr = exp( -d1 -d2 +d3 -p[pB]*(kp2+kpp2) ) / sqrt(u);
    *txcoefr = (float) ( cos(chip) * xr );
    *txcoefi = (float) ( sin(chip) * xr );

    return;

}  /* end txcoef() */


/*------------------------- invert2D() ----------------------*/
/*
    rearrange pix with corners moved to center (a la FFT's)

    pix = complex image
    nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void invert2D( cfpix &pix, int nx, int ny )
{
#define SWAP(a,b)       {t=a; a=b; b=t;}

    int ix, iy, ixmid, iymid;
    float t;

    ixmid = nx/2;
    iymid = ny/2;

    for( ix=0; ix<nx; ix++) 
    for( iy=0; iy<iymid; iy++) {
                SWAP( pix.re(ix,iy), pix.re(ix,iy+iymid) );
                SWAP( pix.im(ix,iy), pix.im(ix,iy+iymid) );
    }

    for( ix=0; ix<ixmid; ix++) 
    for( iy=0; iy<ny; iy++) {
                SWAP( pix.re(ix,iy), pix.re(ix+ixmid,iy) );
                SWAP( pix.im(ix,iy), pix.im(ix+ixmid,iy) );
    }

#undef SWAP
}  /* end invert2D()  */

