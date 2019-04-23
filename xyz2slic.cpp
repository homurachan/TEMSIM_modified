#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output

using namespace std;

#include "cfpix.hpp"       // complex image handler with FFT
#include "slicelib.hpp"    // misc. routines for multislice
#include "floatTIFF.hpp"   // file I/O routines in TIFF format
#include "autoslic.hpp"    //  the calculation engine

#define MANY_ABERR      //  define to include many aberrations
#define eulerfile_line 500000

#ifdef USE_OPENMP
#include <omp.h>
/*  get wall time for benchmarking openMP */
#define walltim() ( omp_get_wtime() )
double walltimer;
#endif

const int NSMAX= 1000;   // max number of slices
const int NCMAX= 1024;   // max characters in file names
const int NZMAX= 103;    // max atomic number Z

int main()
{
    string filein, fileout, filestart, filebeam, filecross, cline, description;
	string filein0;
	std::stringstream stream0;
    const char version[] = "21-feb-2015 (ejk)";

	//  change labPhase to calculate the aberration phase image
    //int lstart=0, lpartl=0, lbeams=0, lwobble=0, lcross=0, nwobble=1, labPhase=1;
    int lstart=0, lpartl=0, lbeams=0, lwobble=0, lcross=0, nwobble=1, labPhase=0;
    int ix, iy, iz, nx, ny, nzout, i, nslic0, islice,
        ndf, nbout, ib, ncellx, ncelly, ncellz, NPARAM;
    int nillum, nzbeams;
    int *hbeam, *kbeam;
    int natom, done, status, multiMode;
    long  ltime;

    unsigned long iseed;

    float v0, mm0, wavlen, rx, ry, ax, by, cz, pi,
        rmin, rmax, aimin, aimax, ctiltx, ctilty,
        acmin, acmax, Cs3, Cs5, df0, sigmaf, dfdelt, aobj,
        temperature, ycross, dx, dy;
   
    float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax,wmin_copy, wmax_copy, xmin_copy,xmax_copy, ymin_copy, ymax_copy, zmin_copy, zmax_copy;
    float *param, *sparam;
    double timer, deltaz, vz;
    ofstream fp1;
    cfpix pix;		// to get results of calculation
    cfpix wave0;	// initial wavefunction (if used)
    cfpix depthpix;	// to get xz cross section results
    cfpix beams;	// to get valuse of requested beams during propagation
    floatTIFF myFile;	//  file input/output
    autoslic aslice;	// has the calculation engine
	std::string str[eulerfile_line+1]; 

	int jmin,jmax;
	cout<<"Enter jmin jmax"<<endl;
	cin >>jmin>>jmax;
    cout << "Enter name of tuned xyz file" << endl;
    cin >> filein0;
	cout << "Enter name of out files (without .tif)"<<endl;
	cin>>filein;
	std::ifstream infile("eulerfile.txt");
	int countN=0;
	float theta[eulerfile_line+1],phy[eulerfile_line+1],res[eulerfile_line+1];
	while(!infile.eof()){
		getline(infile, str[countN],'\n');
		sscanf(((str[countN]).c_str()),"%f\t%f\t%f",&theta[countN],&phy[countN],&res[countN]);
		countN++;
    }
	pi = (float) (4.0 * atan( 1.0 ));
	int j=0;
	v0=300;
    cout << "Wavefunction size in pixels, Nx,Ny:" << endl;
    cin >> nx >> ny;
	cout << "Slice thickness (in Angstroms):" << endl;
    cin >> deltaz;
    if( deltaz < 1.0 ) {
        cout << "WARNING: this slice thickness is probably too thin"
            << " for autoslice to work properly." << endl;
    }
	
	float *x, *y, *z, *occ, *wobble;
	
	int *Znum;
    NPARAM = myFile.maxParam();
    param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
    sparam = (float*) malloc1D( NPARAM, sizeof(float), "sparam" );
    for( ix=0; ix<NPARAM; ix++ ) 
		param[ix] = 0.0F;

	ncellx = 1;
	ncelly = 1;
	ncellz = 1;
	lpartl=0;
    acmin = acmax = 0;
	lstart=0;
    ctiltx = 0;
    ctilty = 0;
    if( lpartl == 0 ) {
		lbeams=0;
    }
	lwobble=0;
	temperature = 0.0F;
    if( lpartl == 0 ) {
		lcross=0;
    }
	nslic0 = 0;   

    mm0 = 1.0F + v0/511.0F;
    wavlen = (float) wavelength( v0 );
    natom = ReadXYZcoord( filein0.c_str(), ncellx, ncelly, ncellz,
        &ax, &by, &cz, &Znum, &x, &y, &z, &occ, &wobble,
        description );
//	cout<<"DEBUG1"<<endl;
	float *x_copy,*y_copy,*z_copy,*occ_copy,*wobble_copy;	
	int *Znum_copy; 
/*	memcpy(&x_copy,&x[0],natom*sizeof(float));
	memcpy(&y_copy,&y[0],natom*sizeof(float));
	memcpy(&z_copy,&z[0],natom*sizeof(float));
	memcpy(&occ_copy,&occ[0],natom*sizeof(float));
	memcpy(&wobble_copy,&wobble[0],natom*sizeof(float));
	memcpy(&Znum_copy,&Znum[0],natom*sizeof(int));*/
	x_copy=(float*) malloc1D( natom, sizeof(float), "x_copy" );
	y_copy=(float*) malloc1D( natom, sizeof(float), "y_copy" );
	z_copy=(float*) malloc1D( natom, sizeof(float), "z_copy" );
	occ_copy=(float*) malloc1D( natom, sizeof(float), "occ_copy" );
	wobble_copy=(float*) malloc1D( natom, sizeof(float), "wobble_copy" );
	Znum_copy=(int*) malloc1D( natom, sizeof(int), "Znum_copy" );
	for(i=0;i<natom;i++){
		x_copy[i]=x[i];
		y_copy[i]=y[i];
		z_copy[i]=z[i];
		occ_copy[i]=occ[i];
		wobble_copy[i]=wobble[i];
		Znum_copy[i]=Znum[i];
	}
	
    cout << natom << " atomic coordinates read in"  << endl;
    cout << description << endl;

    cout <<"Size in pixels Nx, Ny= " << nx << " x " << ny << " = " << nx*ny 
        << " beams" << endl;
    cout <<"Lattice constant a,b = " << ax << ", " << by << endl;

    /*  calculate the total specimen volume and echo */
    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( x[i] < xmin ) xmin = x[i];
        if( x[i] > xmax ) xmax = x[i];
        if( y[i] < ymin ) ymin = y[i];
        if( y[i] > ymax ) ymax = y[i];
        if( z[i] < zmin ) zmin = z[i];
        if( z[i] > zmax ) zmax = z[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
    cout << "Total specimen range is\n" 
        << xmin << " to " << xmax << " in x\n"
            << ymin << " to " << ymax << " in y\n"
        << zmin << " to " << zmax << " in z" << endl;
		
	xmin_copy=xmin;
	xmax_copy=xmax;
	ymin_copy=ymin;
	ymax_copy=ymax;
	zmin_copy=zmin;
	zmax_copy=zmax;
	wmin_copy=wmin;
	wmax_copy=wmax;
	
	float rotat,tilt,rr,xoff,yoff,zoff,x0,y0,z0,xpos,ypos,zpos;
	double cr,sr,ct,st,crr,srr;
//	  xoff = 0.5F*( xmax + xmin );
//    yoff = 0.5F*( ymax + ymin );
//    zoff = 0.5F*( zmax + zmin );
	xoff=0.5F*ax;
	yoff=0.5F*by;
	zoff=0.5F*cz;
//	cout<<"DEBUG2"<<endl;

for(j=jmin-1;j<jmax;j++) {

    rotat = -1*phy[j] * pi/180.0;	/////////////////////
    tilt  = (180-theta[j]) * pi/180.0;	////////////////////
	rr=(180+res[j])*pi/180;
//	rr=-1*res[j]*pi/180.0;
	cr = cos( rotat );
    sr = sin( rotat );
    ct = cos( tilt );
    st = sin( tilt );
	crr= cos(rr);
	srr= sin(rr);
		for( i=0; i<natom; i++) {
            x[i] = x[i] - xoff;	// translate to center
            y[i] = y[i] - yoff;
            z[i] = z[i] - zoff;

            x0 = x[i];      //  Rotation
            y0 = y[i];
            x[i] =	cr*x0 - sr*y0;
			xpos = x[i];
            y[i] =	sr*x0 + cr*y0;

            y0 = y[i];      //  Tilt
            z0 = z[i];
            y[i] = ct*y0 + st*z0;
			ypos =  y[i];
            z[i] = -st*y0 + ct*z0;
			zpos = z[i];
			
			x0 = x[i];      //  Rotation 2
            y0 = y[i];
            x[i] =	crr*x0 - srr*y0;
			xpos = x[i];
            y[i] =	srr*x0 + crr*y0;
/*            if( (fabs( y0-y[i]) > 0.1) || (fabs(z0-z[i])>0.1 ) )

            //  find new range
            if( i==0 ) {		////////////////////	not np here
                xmin = xmax = xpos;
                ymin = ymax = ypos;
                zmin = zmax = zpos;
            } 
			else {
                if( xpos < xmin ) xmin = xpos;
                if( xpos > xmax ) xmax = xpos;
                if( ypos < ymin ) ymin = ypos;
                if( ypos > ymax ) ymax = ypos;
                if( zpos < zmin ) zmin = zpos;
                if( zpos > zmax ) zmax = zpos;
            }	*/
        }
		for(i=0;i<natom;i++){
			x[i] = x[i]+xoff;
			y[i] = y[i]+yoff;
			z[i] = z[i]+zoff;
		}

    aslice.lbeams = lbeams;
    aslice.lcross = lcross;
    aslice.lpartl = lpartl;
    aslice.lstart = lstart;
    aslice.lwobble = lwobble;

    //   set calculation parameters (some already set above)
    param[ pAX ] = ax;			// supercell size
    param[ pBY ] = by;
    param[ pNX ] = (float) nx;
    param[ pNY ] = (float) ny;
    param[pENERGY]   =  v0;
    param[ pDELTAZ ] = (float) deltaz;	// slice thickness
    param[ pOAPERT ] = aobj;
    param[ pXCTILT ] = ctiltx;		// crystal tilt
    param[ pYCTILT ] = ctilty;
    param[pCAPERT] = acmax;		// condencer angles
    param[pCAPERTMIN] = acmin;
    param[ pTEMPER ] = (float) fabs( temperature );
    param[ pNWOBBLE ] = (float) nwobble;	//  number config. to average
    param[pWAVEL] = wavlen;			//  probably recal. autoslice::calculate()

    param[pMODE] = 6;  // save mode = autoslic


    // ------- iterate the multislice algorithm proper -----------

    aslice.calculate( pix, wave0, depthpix, param, multiMode, natom, &iseed,
                Znum, x,y,z,occ,wobble, beams, hbeam, kbeam, nbout, ycross, dfdelt );
 
    pix.findRange( rmin, rmax, aimin, aimax );

    param[pRMAX]  = rmax;
    param[pIMAX]  = aimax;
    param[pRMIN]  = rmin;
    param[pIMIN]  = aimin;

    param[pDX] = dx = (float) ( ax/((float)nx) );
    param[pDY] = dy = (float) ( by/((float)ny) );

    //param[pNSLICES] = 0.0F;  /* ??? */

    for( ix=0; ix<NPARAM; ix++ ) 
		myFile.setParam( ix, param[ix] );
    myFile.resize( 2*nx, ny );
    myFile.setnpix( 2 );
    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
        myFile(ix,iy)    = pix.re(ix,iy);
        myFile(ix+nx,iy) = pix.im(ix,iy);
    }
	
	stream0 << setfill('0');
	stream0 << setw(4);
	stream0 << j+1;
	fileout="ASed_"+filein+"_"+stream0.str()+"-tif";
    i = myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax, dx, dy );
    if( i != 1 ) cout << "autoslice cannot write TIF file " << fileout << endl;
    cout << "pix range " << rmin << " to " << rmax << " real,\n" <<
            "          " << aimin << " to " << aimax << " imag" << endl;

	for(i=0;i<natom;i++){
		x[i]=x_copy[i];
		y[i]=y_copy[i];
		z[i]=z_copy[i];
		occ[i]=occ_copy[i];
		wobble[i]=wobble_copy[i];
		Znum[i]=Znum_copy[i];
	}
	xmin=xmin_copy;
	xmax=xmax_copy;
	ymin=ymin_copy;
	ymax=ymax_copy;
	zmin=zmin_copy;
	zmax=zmax_copy;
	wmin=wmin_copy;
	wmax=wmax_copy;
	stream0.str("");
	printf("finish iteration %d of %d\n\n",j+1,jmax-jmin+1);
}
    return 0;

} /* end main() */



