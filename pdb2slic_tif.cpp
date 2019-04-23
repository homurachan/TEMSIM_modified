#include <cstdio>  /* ANSI C libraries */
#include <cstring>
#include <ctime>

#include "slicelib.hpp"
#include "cfpix.hpp"       // complex image handler with FFT
#include "floatTIFF.hpp"   // file I/O routines in TIFF format
#include "autoslic.hpp"    //  the calculation engine

#define MANY_ABERR      //  define to include many aberrations


#include <cstdlib>
#include <cmath>
#include <ctype.h>   // for toupper()  */
#include <string.h>  /* for strings */
#include <time.h>

#include <string>
#include <iostream>  //  C++ stream IO
#include <sstream>   // string streams
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>
#define eulerfile_line 500000

using namespace std;



#ifdef USE_OPENMP
#include <omp.h>
/*  get wall time for benchmarking openMP */
#define walltim() ( omp_get_wtime() )
double walltimer;
#endif

const int NSMAX= 1000;   // max number of slices
const int NCMAX= 1024;   // max characters in file names
const int NZMAX= 103;    // max atomic number Z

const int TRUE=1;
const int FALSE=0;

const double CWEIGHT=16;	// molec. weight of carbon	-- modified to oxygen
const double CDENSITY=0.84;		// density in gm/cm^3 approx. for amorphous C	-- modified here to LDA
const double NAV=6.0225e23;		// Avagadro's number (#/mole)
const double RMIN=1.0;			// min separation distance of O (in Angstroms)

const int OX=0, OY=1, OZ=2;	// index for Oxygen coordinates
const int NVAL=3;

unsigned long iseed;		// random number generator seed

/*  subroutines define at the end of the file */
void newcoord( double c[], 
    const double xmax, const double ymax, const double zmax );
void insertcoord( double ctest[], double** coord, 
    const int ncoord, const int pos );
int testcoord( double ctest[], double** coord, const int ncoord, int *pos,
    const double rmin );
int fillSolid( double** coord, const int ncoord, 
    double ax, double by, double cz, double rmin );


int main()
{
    string cline, filin, filout, symb;
    int i, ic, iz,  np, nlen, ntotal, ctop,j,k,np_bak;
    long ltime;
    double xpos, ypos, zpos, occ, mytime, x0, y0, z0, cr,sr,ct,st;
    double xmin, xmax, ymin, ymax, zmin, zmax,xmin_bak, xmax_bak, ymin_bak, ymax_bak, zmin_bak, zmax_bak, rotat, tilt, pi;
	float theta[eulerfile_line+1],phy[eulerfile_line+1],res[eulerfile_line+1];
    double ax, by, cz, cthick, xoff, yoff, zoff, wobble, density;
    double **coord;			// coordinates
	int jmin,jmax;
	printf("Enter jmin jmax\n");
	scanf("%d %d",&jmin,&jmax);
	if(jmax<jmin)
	{
		printf("max smaller than min\n");
		return 0;
	}
    vector<int> znum,znum_bak;
    vector<double> xp, yp, zp, oc,xp_bak,yp_bak,zp_bak,oc_bak,wobblev;
	
	std::string outfilename[jmax+1];
	std::stringstream stream0;
	std::ifstream infile("eulerfile.txt");
    ifstream fpin;
    ofstream fpout;
    std::string str[eulerfile_line+1]; 
	int jmax0=jmax-jmin+1;
	
	cfpix pix;		// to get results of calculation
    cfpix wave0;	// initial wavefunction (if used)
    cfpix depthpix;	// to get xz cross section results
    cfpix beams;	// to get valuse of requested beams during propagation

    floatTIFF myFile;	//  file input/output
    autoslic aslice;	// has the calculation engine
    /*  the following are the chemical symbols for the periodic table */
    char symbol[] = {
                " HHeLiBe B C N O FNeNaMgAlSi P SCl"
                "Ar KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBr"
                "KrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe"
                " IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
                "YbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRn"
                "FrRaAcThPa UNpPuAmCmBkCfEsFmMdNoLr"
        };

    /*------- what this is ---------- */
    pi = 4.0 * atan( 1.0 );
    cout << "pdb2xyz version dated 6-may-2015" << endl;
    cout << "  convert PDB data file to xyz format\n" << endl;
    cout << "  (may not work on all PDB files)\n" << endl;

    //------- convert symbols to upper case for PDB comparison ---------- 
    
    nlen = (int) strlen( symbol );
    for( i=0; i<nlen; i++) symbol[i] = toupper( symbol[i] );

    //----------- get file names and open ---- 
    
    cout << "Type name of input file with PDB data:" << endl;
    cin >> filin;
    fpin.open( filin.c_str() );
    if( fpin.bad() ) { cout << "Can't open file "<< filin<<endl; exit(0); }

    cout << "Type name of output filename (without .tif) to get tif data:" << endl;
    cin >>  outfilename[0];
	float v0=300;
	int nx,ny;
    cout << "Wavefunction size in pixels, Nx,Ny:" << endl;
    cin >> nx >> ny;
	double deltaz;
	cout << "Slice thickness (in Angstroms):" << endl;
    cin >> deltaz;
    if( deltaz < 1.0 ) {
        cout << "WARNING: this slice thickness is probably too thin"
            << " for autoslice to work properly." << endl;
    }
	
    int countN=0;
	while(!infile.eof()){
		getline(infile, str[countN],'\n');
		sscanf(((str[countN]).c_str()),"%f\t%f\t%f",&theta[countN],&phy[countN],&res[countN]);
		
		//////////////////////////////////////////////////////////////
		//	  READ IN EULERANGLE, NEED TO CHANGE TO RIGHT ANGLE 	//
		//////////////////////////////////////////////////////////////
		
		
//		printf("%f %f %f %d\n",theta[countN],phy[countN],res[countN],countN);
		countN++;
    }

	for(j=1;j<jmax+1;j++){
		stream0<<j;
		outfilename[j]="ASed_"+outfilename[0]+"_"+stream0.str()+"-tif";
		stream0.str("");
//		cout << outfilename[j] << endl;
	}

    ///----------- get carbon support info ---- 
/*    cout << "Type thickness of carbon support (<0 to disable):" << endl;
    cin >> cthick;
    cthick = fabs( cthick);			*/
	cthick = 0;
    if( cthick > 0.0 ) 
        ctop = askYN( "Do you want the oxygen support on the entrance instead of exit");

    //  get CPU time for fun
    mytime = cputim();

    //--------- read data from file in complicated format -------
    //   first get total range and number of atoms 
    np = 0;	
    do {
        getline( fpin, cline );	// read a whole line

        //  select lines begin. with ATOM or HETAM  with atom coord.
        if( ( cline.find( "ATOM") == 0 ) ||
            ( cline.find( "HETATM") == 0 )  )  {
            
            //----  read x,y,z coordinates
            istringstream isbuf( cline.substr(30) );
            isbuf >> xpos >> ypos >> zpos >> occ;

            if( 0 == np ) {
                xmin = xmax = xpos;
                ymin = ymax = ypos;
                zmin = zmax = zpos;
            } else {
                if( xpos < xmin ) xmin = xpos;  // coord. range
                if( xpos > xmax ) xmax = xpos;
                if( ypos < ymin ) ymin = ypos;
                if( ypos > ymax ) ymax = ypos;
                if( zpos < zmin ) zmin = zpos;
                if( zpos > zmax ) zmax = zpos;
            }
            xp.push_back( xpos );	//  save coord.
			xp_bak.push_back( xpos );
            yp.push_back( ypos );
			yp_bak.push_back( ypos );
            zp.push_back( zpos );
			zp_bak.push_back( zpos );
            oc.push_back( occ );
			oc_bak.push_back( occ );
			wobblev.push_back(0);
            //  find atomic number
            symb= cline.substr(76,77);  // get chemical symbol
            for( i=0; i<nlen; i+=2) {
                iz = 1 + i/2;
                if( strncmp( &symbol[i], symb.c_str(), 2 ) == 0 ) break;
            }
            znum.push_back( iz );
			znum_bak.push_back( iz );
            np++;
        }

    } 
	while( cline.find( "END") == string::npos );
    fpin.close( );
    
    cout <<  "Total number of atoms = " << np << endl;
    cout <<  "  with x range " << xmin << " to " << xmax << endl;
    cout <<  "   and y range " << ymin << " to " << ymax << endl;
    cout <<  "   and z range " << zmin << " to " << zmax << endl;

	xmax_bak=xmax;
	xmin_bak=xmin;
	ymin_bak=ymin;
	ymax_bak=ymax;
	zmax_bak=zmax;
	zmin_bak=zmin;
	np_bak=np;
	
for(j=jmin-1;j<jmax;j++)	{
	
	filout=outfilename[j+1];
	fpout.open( filout.c_str() );
    if( fpout.bad() ) { cout << "Can't open file" << filout<<endl; break; }
    rotat = -1*phy[j] * pi/180.0;	/////////////////////
    tilt  = (180-theta[j]) * pi/180.0;	////////////////////
	
    if( ( fabs( rotat) > 1.0e-6 ) || ( fabs( tilt ) > 1.0e-6 ) ) {
        //-----  rotate to get better view perhaps (from slicview.cpp)
        //  Move to center of molecule and rotate
        xoff = 0.5F*( xmax + xmin );
        yoff = 0.5F*( ymax + ymin );
        zoff = 0.5F*( zmax + zmin );

        cout << "rotate about x,y,z= " << xpos << ", " << ypos << ", " << zpos << endl;
		cout << "rotat="<<theta[j]<<" tilt="<<phy[j]<<endl;
        /*  Calculate misc constants  */
        cr = cos( rotat );
        sr = sin( rotat );
        ct = cos( tilt );
        st = sin( tilt );

        for( i=0; i<np; i++) {
            xp[i] = xp[i] - xoff;	// translate to center
            yp[i] = yp[i] - yoff;
            zp[i] = zp[i] - zoff;

            x0 = xp[i];      //  Rotation
            y0 = yp[i];
            xp[i] =	cr*x0 - sr*y0;
			xpos = xp[i];
            yp[i] =	sr*x0 + cr*y0;

            y0 = yp[i];      //  Tilt
            z0 = zp[i];
            yp[i] = ct*y0 + st*z0;
			ypos =  yp[i];
            zp[i] = -st*y0 + ct*z0;
			zpos = zp[i];
            if( (fabs( y0-yp[i]) > 0.1) || (fabs(z0-zp[i])>0.1 ) )

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
            }
        }

        cout <<  "New range of atom coord. after rotation and tilt" << endl;
        cout <<  "   x range " << xmin << " to " << xmax << endl;
        cout <<  "   y range " << ymin << " to " << ymax << endl;
        cout <<  "   z range " << zmin << " to " << zmax << endl;

    }  // end rotate section

    //--------- write data to xyz file with offset -------

    //ax = 1.4*( xmax - xmin );
    //by = 1.4*( ymax - ymin );
    ax = 1.4*( xmax_bak - xmin_bak );		//////////////////
    by = 1.4*( ymax_bak - ymin_bak );		//////////////////
    if( by > ax ) ax = by;		// make it square to look right 
    else by = ax;
    if( cthick > 0.0 ) cz = (zmax_bak - zmin_bak) + cthick;
    else cz = zmax_bak - zmin_bak;
	
	double ax0=1.4*( xmax - xmin ),by0=1.4*( ymax - ymin );		///////////////////////////////////////
	if( by0 > ax0 ) ax0 = by0;		
    else by0 = ax0;
    xoff = 0.5*ax0 - 0.5*(xmax+xmin)+(ax-ax0);   // move molecule to the center 
    yoff = 0.5*ax0 - 0.5*(ymax+ymin)+(ax-ax0);		////////////////	probably right
	
    if( (cthick>0.0) && (ctop==1) ) zoff = -zmin + cthick;
    else zoff = -zmin;
    
    wobble = 0.0;   //  Debye-Waller factor 


    //----------- generate carbon support if requested -------------- 

    if( cthick > 0.0 ) {

        // ---- initialize random number generator seed ---- 
        cout << endl;
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            cout << "Type initial seed for random number generator:" << endl;
            cin >> iseed;
        } else {
            cout <<  "Random number seed initialized to " << iseed << endl;
        }
        cout << endl;

        
        cout << "generate random coord. for oxygen support"  << endl;
        density = NAV*CDENSITY*(1.0e-24)/CWEIGHT; //  # atoms/Angs^3 
        cout << "average density = "<< CDENSITY << " gm/cm^3 = "
            << density << " atoms/Angstrom^3" << endl;
        cout << "minimum allowed separation = " << RMIN << " Angstroms" << endl;

        cout << "calculate coord. in a vol. of "<< ax << " x "<< by << 
            " x "<< cthick << " Angstroms" << endl;
        ntotal = (int) ( ax*by*cthick  * density );
        cout << "Total number of oxygen atoms = " << ntotal << endl;

        coord = (double**) malloc2D( ntotal, NVAL, sizeof(double), "coord" );

        //  fill in random coord.
        ic = fillSolid( coord, ntotal, ax, by, cthick, RMIN );

        iz = 8;   //  atomix number of carbon	-- modified to oxygen
        occ = 1.0;
        wobble = 0.0;
        if( ctop == 1 ) zoff = 0.0;	//  on top  
        else zoff = zmax - zmin;	// on bottom 

    }

//////////////////////////////////////////////////////////////////////////////////////////////
//									START AUTOSLIC											//
//////////////////////////////////////////////////////////////////////////////////////////////
	float *x, *y, *z, *occ1, *wobble1;
	int NPARAM = myFile.maxParam();
	int *Znum;
    float *param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
    float *sparam = (float*) malloc1D( NPARAM, sizeof(float), "sparam" );
	int ix,iy,l;
	for( ix=0; ix<NPARAM; ix++ ) param[ix] = 0.0F;
	int ncellx = 1,ncelly = 1,ncellz = 1;
	int lpartl=0;
    float acmin = 0,acmax = 0;
	int lstart=0;
    float ctiltx = 0,ctilty = 0;
	int lbeams=0;
	int lwobble=0;
	float temperature = 0.0F;
	int lcross=0;
	int nslic0 = 0;
	float mm0 = 1.0F + v0/511.0F;
	float wavlen = (float) wavelength( v0 );
	cout << "electron wavelength = " << wavlen << " Angstroms" << endl;
	int natom=np;
	x  = (float*) malloc1D( natom, sizeof(float), "x" );
    y  = (float*) malloc1D( natom, sizeof(float), "y" );
    z  = (float*) malloc1D( natom, sizeof(float), "z" );
    occ1    = (float*) malloc1D( natom, sizeof(float), "occ" );
    wobble1 = (float*) malloc1D( natom, sizeof(float), "wobble" );
    Znum   =   (int*) malloc1D( natom, sizeof(int), "Znum" );
	for( l=0; l< natom; l++ ) {
        x[l] = float(xp[l]+xoff);
        y[l] = float(yp[l])+yoff;
        z[l] = float(zp[l]+zoff);
        occ1[l] = float(oc[l]);
		wobble1[l] = float(wobblev[l]);
        Znum[l] = float(znum[l]);
    }
	float aobj,ycross,dfdelt,rmin1, rmax, aimin, aimax,dx1,dy1;
	int nwobble=1;
	int multiMode,nbout;
	int *hbeam, *kbeam;

	cout << natom << " atomic coordinates read in"  << endl;
    aslice.lbeams = lbeams;
    aslice.lcross = lcross;
    aslice.lpartl = lpartl;
    aslice.lstart = lstart;
    aslice.lwobble = lwobble;
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

	aslice.calculate( pix, wave0, depthpix, param, multiMode, natom, &iseed,
        Znum, x,y,z,occ1,wobble1, beams, hbeam, kbeam, nbout, ycross, dfdelt );
 
    pix.findRange( rmin1, rmax, aimin, aimax );
	param[pRMAX]  = rmax;
    param[pIMAX]  = aimax;
    param[pRMIN]  = rmin1;
    param[pIMIN]  = aimin;

	param[pDX] = dx1 = (float) ( ax/((float)nx) );
    param[pDY] = dy1 = (float) ( by/((float)ny) );
    for( ix=0; ix<NPARAM; ix++ ) myFile.setParam( ix, param[ix] );
	myFile.resize( 2*nx, ny );
	myFile.setnpix( 2 );
    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
        myFile(ix,iy)    = pix.re(ix,iy);
        myFile(ix+nx,iy) = pix.im(ix,iy);
	}
    i = myFile.write( outfilename[j+1].c_str(), rmin1, rmax, aimin, aimax, dx1, dy1 );
    if( i != 1 ) cout << "autoslice cannot write TIF file " << outfilename[j+1] << endl;
    cout << "pix range " << rmin1 << " to " << rmax << " real,\n" <<
            "          " << aimin << " to " << aimax << " imag" << endl;
	free(x);
	free(y);
	free(z);
	free(occ1);
	free(Znum);
	free(param);
	free(sparam);
	x=NULL;
	y=NULL;
	z=NULL;
	occ1=NULL;
	Znum=NULL;
	param=NULL;
	sparam=NULL;
	
	memcpy(&xp[0],&xp_bak[0],xp_bak.size()*sizeof(double));
	memcpy(&yp[0],&yp_bak[0],yp_bak.size()*sizeof(double));
	memcpy(&zp[0],&zp_bak[0],zp_bak.size()*sizeof(double));
	memcpy(&oc[0],&oc_bak[0],oc_bak.size()*sizeof(double));
	memcpy(&znum[0],&znum_bak[0],znum.size()*sizeof(int));
	
	xmax=xmax_bak;
	xmin=xmin_bak;
	ymin=ymin_bak;
	ymax=ymax_bak;
	zmax=zmax_bak;
	zmin=zmin_bak;
	np=np_bak;
	
	mytime = cputim() - mytime;
    cout << "\ntotal CPU time = " << mytime << " sec." << endl;
	printf("Finish iteration %d of %d\n\n",j+1,jmax0);
}
    return( EXIT_SUCCESS );

}  // end main 


/*------------------------ fillSolid() ------------------------*/
/*
    fill a vol. with an amorphous solid with random coord.
    in the range (0,0,0) to (ax,by,cz)

  coord[][] = to get list of existing sorted by z = coord[][OZ]
                dimensions ncoord x NVAL
  ncoord = number of coordinates to generate
  rmin = minimum separation distance

  assumed globals
     NVAL

*/
int fillSolid( double** coord, const int ncoord, 
    double ax, double by, double cz, double rmin )
{
    int ic, i;
    double ctest[NVAL];

    newcoord( coord[0], ax, by, cz );   // start coord. list
    ic = 1;
    do {
        // get new coordinate 
        newcoord( ctest, ax, by, cz );
        //cout << "new coord= ", << ctest[0] << ", " << 
        //           ctest[1] << ", " <<ctest[2] << endl;

        //  add it to the list if its OK - this also sorts
        if( testcoord( ctest, coord, ic, &i, RMIN ) == TRUE ) {
            insertcoord( ctest, coord, ic, i);
            ic++;
            if( ic%1000 == 0 ) cout << "ic= \r" << ic;
        }
    }  while( ic < ncoord );

    return ic;

}  // end fillSolid()

/*-----  subroutines below are from biocells.cpp ---------------- 
   with long -> int  */

/*------------------------ insertcoord() ------------------------*/
/*
    insert a new coordinate in the list

  ctest[] = new coord to insert
  coord[][] = list of existing coordinates sorted by z = coord[][OZ]
  ncoord = number of coordinates
  pos = position to insert coord. at

  insert new coord at index pos and move all the rest down one

  assumed globals
     NVAL

*/
void insertcoord( double ctest[], double** coord, 
        const int ncoord, const int pos )
{
    long i, j;
    
    for( i=ncoord; i>pos; i--) {
        for( j=0; j<NVAL; j++)
            coord[i][j] = coord[i-1][j];
    }

    for( j=0; j<NVAL; j++)
            coord[pos][j] = ctest[j];

} // end insertcoord() 


/*------------------------ newcoord() ------------------------*/
/*
    generate a new random coordinate inside the required volume

    xmax, ymax, zmax = volume size
*/
void newcoord( double c[], 
    const double xmax, const double ymax, const double zmax )
{
    c[OX] = xmax * ranflat( &iseed );
    c[OY] = ymax * ranflat( &iseed );
    c[OZ] = zmax * ranflat( &iseed );
    return;

} // end newcoord() 

/*------------------------ testcoord() ------------------------*/
/*
    test the new coordinate to see if its too close to an existing
    coordinate (closer than RMIN)

    A straight search of the whole list is very slow (proportional to N^2)
    and was not pratical for a large set of coordinates (i.e. took
    an absurd amount of CPU time).

    This routines assumes that the list of existing coordinates is
    sorted wrt one coord (z in this case).  Once the position of the new
    coordinate is located in the list (using a binary search), then
    only the small range of coord. near this point need to be
    tested.  This makes the test dramatically faster.

  ctest[] = new coord to test
  coord[][] = list of existing coordinates sorted by z = coord[][OZ]
  ncoord = number of coordinates
  pos = returned position to insert coord.

  assumed globals
     OX, OY, OZ, RMIN

*/
int testcoord( double ctest[], double** coord, const int ncoord, int *pos,
              const double rmin )
{
    

#ifdef SLOW

    /*  this is VERY SLOW - do NOT use */
    /*  keep for comparison */

    int good;
    long i;
    double d, dx, dy, dz, rmin2;

    rmin2 = rmin*rmin;

        /*  test if this new site is occupied */
        good = TRUE;
        d = 0;
        for( i=0; i<ncoord; i++) {
            dx = ctest[OX] - coord[i][OX];
            dy = ctest[OY] - coord[i][OY];
            dz = ctest[OZ] - coord[i][OZ];
            d = dx*dx + dy*dy + dz*dz;
            if( d <= rmin2 ) {
                good = FALSE;
                break;
            }
        }
        *pos = ncoord;
        return( good );

#else
        /* this is the more -sophisticated version sorted by Z */

    long i, j, k;
    double d, dx, dy, dz, dz2, rmin2, z, range;

    /* ---- find postion of ctest[] in coord[][] using a binary search  --- 
         i should get the position to insert (between i and i+1)
         this assumes coord[][] is sorted by coord OZ
    */
//	printf("testcoord() top, ncoord= %d\n", ncoord );

    z = ctest[OZ];
    if( z <= coord[0][OZ] ) i = j = 0;
    else if( z >= coord[ncoord-1][OZ] ) i = j = ncoord;
    else { 
        i = 0;
        j = ncoord-1;
        do{ k = ( i + j ) / 2 ;
            if( z < coord[k][OZ] )  j = k;
            else if( z >=  coord[k][OZ] ) i = k;
        } while ( (j-i) > 1 );
    }

//	printf("testcoord() 4, z= %f, i= %d, j= %d\n", z, i, j );

    /* now that we have the position of the new point
       we only have to explore within RMIN of this point */

    rmin2 = rmin*rmin;
    range = 2.0*rmin2;   /* add a little safety margin */
    k = j;
    while( k >= 0 ) {
            dx = ctest[OX] - coord[k][OX];
            dy = ctest[OY] - coord[k][OY];
            dz = ctest[OZ] - coord[k][OZ];
            dz2 = dz*dz;
            d = dx*dx + dy*dy + dz2;
            if( d <= rmin2 ) return( FALSE );
            if( dz2 > range ) break;
            k--;
    }
    k = j;
    while( k < ncoord ) {
            dx = ctest[OX] - coord[k][OX];
            dy = ctest[OY] - coord[k][OY];
            dz = ctest[OZ] - coord[k][OZ];
            dz2 = dz*dz;
            d = dx*dx + dy*dy + dz2;
            if( d <= rmin2 ) return( FALSE );
            if( dz2 > range ) break;
            k++;
    }

    /* if it gets to here then this is a good point */

    *pos = j;
/*???
    for( k=0; k<ncoord; k++) {
        printf("coord[%d] = %g, %g, %g\n", k,
            coord[k][OX], coord[k][OY], coord[k][OZ] );
    }
    printf( "ctest= %g, %g, %g, insert at pos= %d\n",
        ctest[OX], ctest[OY], ctest[OZ], *pos );
    scanf( "%d", &k );
*/
    return( TRUE );

#endif

} // end testcoord() 

