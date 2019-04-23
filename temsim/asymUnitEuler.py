#!/usr/bin/env python

# $Id: asymUnitEuler.py 828 2010-06-22 13:47:21Z jiang12 $

# by Wen Jiang <jiang12@purdue.edu>, 2005-9-25
# this program will generate uniform Euler angles in an asymmetric unit using the 
# "Recursive Zonal Equal Area Sphere Partitioning" algorithm implemented in Matlab
# See reference: http://eqsp.sourceforge.net

import math, os, sys

try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(options, angstep, eulerfile) =  parse_command_line()
	
	points = asym_unit(step = angstep * math.pi/180., sym=options.sym, min_ang_dist = options.mindist*math.pi/180. )
	points2eulerfile(eulerfile, points)
	
	if options.tiltdistfile:
		plot_tilt_angle_distribution(options.tiltdistfile, points, binsize=1 * math.pi/180.)
	if options.pdbfile:
		points2pdb(options.pdbfile, points)
	
def parse_command_line():
	usage="%prog <ang step size in degree> <output euler file> [options]"
	parser = OptionParser(usage=usage, version="%prog v1.0, September 2005. By Wen Jiang <jiang12@purdue.edu>\nNOTE: requires matlab")
	parser.add_option("--sym",dest="sym",type="string",metavar="<cn|dn|icos|tet>",help="symmetry to use, default to c1", default="c1")
	parser.add_option("--min_bound_dist",dest="mindist",type="float",metavar="<degree>",help="minimal angular distance to asymmetric unit boundary to include the orientation. default to 0. can be + or -", default=0)
	parser.add_option("--pdb",dest="pdbfile",type="string",metavar="<filename>",help="to output the orientation on unit sphere in PDB file format", default="")
	parser.add_option("--tiltdist",dest="tiltdistfile",type="string",metavar="<filename>",help="to output tilt angle distribution in a text file", default="")
	
	if len(sys.argv)<3: 
		parser.print_help()
		sys.exit(-1)
	
	(options, args)=parser.parse_args()

	if len(args) !=2:
		parser.print_help()
		sys.exit(-1)
	
	angstep = float(args[0])
	eulerfile = args[1]
	
	return (options, angstep, eulerfile)

def asym_unit(step = 1.* math.pi/180., sym="c1", min_ang_dist=0):
	mindist = math.tan(min_ang_dist)
	N = int(4*math.pi/(step*step))
	
	scriptfile = "eq_point_set.m"
	fp = open(scriptfile,'w')
	fp.write(eq_point_set_matlab)
	fp.close()
	
	pointsfile = "x_points.txt"
	cmd = "points_x = eq_point_set(2,%d); points_x = points_x'; \n" % (N)
	cmd+= "save %s points_x -ascii;\n" % (pointsfile)
	cmd+= "exit\n"
	#pfp, x, x = os.popen3("matlab -nosplash -nodesktop -nojvm") 
	pfp = os.popen("matlab -nosplash -nodesktop -nojvm > /dev/null", 'w') 
	#pfp = os.popen("matlab -nosplash -nodesktop -nojvm ", 'w') 
	pfp.write(cmd)
	pfp.close()
	
	points = []
	for l in open(pointsfile,'r'):
		points.append([ float(s) for s in l.split() ])
	os.remove(pointsfile)
	os.remove(scriptfile)
	
	sym = sym.lower()
	if sym[0] == 'c':
		cyc = int(sym[1:])
		if cyc<=0: 
			print "WARNING: cyclic symmetry (%s) should be positive. set to c1" % (sym)
			cyc=1
		
		if cyc == 1:	# special for C1
			c0 = (1,0,0)
			c1 = (0,1,0)
			c2 = (-math.sqrt(2)/2,-math.sqrt(2)/2,0)
		else:
			az1=math.pi*2/cyc	
			c0 = (0,0,1)
			c1 = (1,0,0)
			c2 = (math.cos(az1), math.sin(az1), 0)
			
	elif sym[0] == 'd':
		cyc = int(sym[1:])
		if cyc<=0: 
			print "WARNING: diheral symmetry (%s) should be positive. set to d1" % (sym)
			cyc=1
		
		az1=math.pi/cyc

		c0 = (0,0,1)
		c1 = (1,0,0)
		c2 = (math.cos(az1), math.sin(az1), 0)
	elif sym == 'icos':
		ANG5to5HLF	= 	0.55357435889704525151 		# half of edge, 31.71747 degree
		ANG5to3		=	0.652358139784368185995		# 5fold to closest 3fold, 37.37737 degree
		c0 = (0,0,1)
		#c1 = (math.sin( ANG5to3 ), 0, math.cos( ANG5to3 ))
		c1 = (0, -math.sin( ANG5to3 ), math.cos( ANG5to3 ))
		#c2 = (math.sin( ANG5to5HLF )*math.cos(math.pi/5), math.sin( ANG5to5HLF )*math.sin(math.pi/5), math.cos( ANG5to5HLF ))		
		c2 = (math.sin( ANG5to5HLF )*math.cos(math.pi/5-math.pi/2), math.sin( ANG5to5HLF )*math.sin(math.pi/5-math.pi/2), math.cos( ANG5to5HLF ))
	
	elif sym == 'oct':
		omega   =  math.pi*2.0/4.0
		thetaC  =  math.acos(math.cos(omega)/(1.0-math.cos(omega)))/2.0
		c0 = (0,0,1)
		c1 = (math.sin(thetaC), 0, math.cos(thetaC))             # point to nearest vertex
		c2 = ( math.sqrt(1.0-2.0*math.cos(omega))*math.cos(omega/2.0)/(math.sqrt(3.0)*math.sin(omega/2.0)), \
		       math.sqrt(1.0-2.0*math.cos(omega))*math.sin(omega/2.0)/(math.sqrt(3.0)*math.sin(omega/2.0)), \
		                         math.cos(omega/2.0)/(math.sqrt(3.0)*math.sin(omega/2.0)))

	elif sym == 'tet':
		omega   =  math.pi*2.0/3.0
		thetaC  =  math.acos(math.cos(omega)/(1.0-math.cos(omega)))/2.0
		c0 = (0,0,1)
		c1 = (math.sin(thetaC), 0, math.cos(thetaC))
		c2 = ( math.sqrt(1.0-2.0*math.cos(omega))*math.cos(omega/2.0)/(math.sqrt(3.0)*math.sin(omega/2.0)), \
		       math.sqrt(1.0-2.0*math.cos(omega))*math.sin(omega/2.0)/(math.sqrt(3.0)*math.sin(omega/2.0)), \
		                                          math.cos(omega/2.0)/(math.sqrt(3.0)*math.sin(omega/2.0)))
	
	else:
		raise ValueError("symmetry %s is not supported" % (sym))
	
	def cross(v1, v2):
		return ([v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]])
	def norm(v):
		len2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
		return (v[0] / len2, v[1] / len2, v[2] / len2)
	def distance(p, n):
		return p[0]*n[0] + p[1]*n[1] + p[2]*n[2]
	
	n1 = norm( cross(c0, [c1[0]-c0[0], c1[1]-c0[1], c1[2]-c0[2]]) )
	n2 = norm( cross(c1, [c2[0]-c1[0], c2[1]-c1[1], c2[2]-c1[2]]) )
	n3 = norm( cross(c2, [c0[0]-c2[0], c0[1]-c2[1], c0[2]-c2[2]]) )
	
	points2 = []
	for p in points:
		if distance(p, n1)>=mindist and distance(p, n2)>=mindist and distance(p, n3)>=mindist:
			points2.append(p)
	
	return points2
	
def points2pdb(pdbfile, points):
	pdb = open(pdbfile,"w")
	for i in range(len(points)):
		p = points[i]
		pdb.write("ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%8s\n" % (i,i,p[0]*100.0,p[1]*100.0,p[2]*100.0,1.0,0.0," "))
	pdb.close()

def points2eulerfile(eulerfile, points):
	phi = 0
	eulerfp = open(eulerfile,"w")
	for i in range(len(points)):
		p = points[i]
		alt = math.acos(p[2]) * 180./math.pi
		az = math.atan2(p[1], p[0]) * 180./ math.pi
		eulerfp.write("%8g\t%8g\t%8g\n" % (alt, az, phi))
	eulerfp.close()

def plot_tilt_angle_distribution(txtfile, points, binsize=1 * math.pi/180.):
	def norm(v):
		len2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
		return (v[0] / len2, v[1] / len2, v[2] / len2)
	bins = [0] * (int(math.pi/2/binsize)+1)
	for p in points:
		theta = math.acos(norm(p)[2])
		index = int(theta/binsize)
		if index>=len(bins): index = len(bins)-1
		bins[index]+=1./len(points)
		#print "%8g %8g %8g => %8g %d" % (p[0], p[1], p[2], theta*180./math.pi, index)
	for i in range(1, len(bins)):
		bins[i] += bins[i-1]
	
	outfile = open(txtfile,"w")
	for i in range(len(bins)):
		outfile.write("%8g\t%8g\n" % (i*binsize*180/math.pi, bins[i]))
	outfile.close()
	

eq_point_set_matlab=r"""function points_x = eq_point_set(dim,N,varargin)
% EQ_POINT_SET Center points of regions of EQ partition, in Cartesian coordinates
% Reference: Recursive Zonal Equal Area Sphere Partitioning Toolbox
% URL: http://eqsp.sourceforge.net
% 
% Syntax
% points_x = eq_point_set(dim,N,options);
%
% Description
% POINTS_X = EQ_POINT_SET(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
% partition S^dim (the unit sphere in dim+1 dimensional space) into N regions
% of equal area and small diameter, and
% 2) sets POINTS_X to be an array of size (dim+1 by N), containing the center
% points of each region.
% Each column of POINTS_X represents a point of S^dim, in Cartesian coordinates.
%
% The arguments dim and N must be positive integers.
%
% POINTS_X = EQ_POINT_SET(dim,N,'offset','extra') uses experimental extra offsets
% for S^2 and S^3 to try to minimize energy.
%
% POINTS_X = EQ_POINT_SET(dim,N,extra_offset) uses experimental extra offsets if
% extra_offset is true or non-zero.
%
% Notes
% Each region is defined as a product of intervals in spherical polar
% coordinates. The center point of a region is defined via the center points
% of each interval, with the exception of spherical caps and their descendants,
% where the center point is defined using the center of the spherical cap.
%
% If dim > 3, extra offsets are not used.
% For more details on options, see help partition_options.
%
% Examples
% > points_x = eq_point_set(2,4)
% points_x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
% > size(points_x)
% ans =
%      3     4
%
% See also
% PARTITION_OPTIONS, EQ_POINT_SET_POLAR, EQ_REGIONS, S2X

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
error(nargchk(2,4,nargin));
points_x = polar2cart(eq_point_set_polar(dim,N,varargin{:}));
%
% end function


function points_s = eq_point_set_polar(dim,N,varargin)
%EQ_POINT_SET_POLAR Center points of regions of an EQ partition
%
%Syntax
% points_s = eq_point_set_polar(dim,N,options);
%
%Description
% POINTS_S = EQ_POINT_SET_POLAR(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
% partition S^dim (the unit sphere in dim+1 dimensional space) into N regions
% of equal area and small diameter, and
% 2) sets POINTS_S to be an array of size (dim by N), containing the center
% points of each region. Each column of POINTS_S represents a point of S^dim,
% in spherical polar coordinates.
%
% The arguments dim and N must be positive integers.
%
% POINTS_S = EQ_POINT_SET_POLAR(dim,N,'offset','extra') uses experimental extra
% offsets for S^2 and S^3 to try to minimize energy. If dim > 3, extra offsets
% are not used.
%
% POINTS_S = EQ_POINT_SET_POLAR(dim,N,extra_offset) uses experimental extra
% offsets if extra_offset is true or non-zero.
%
%Notes
% Each region is defined as a product of intervals in spherical polar
% coordinates. The center point of a region is defined via the center points
% of each interval, with the exception of spherical caps and their descendants,
% where the center point is defined using the center of the spherical cap.
%
% For more details on options, see HELP PARTITION_OPTIONS.
%
%Examples
% > points_s = eq_point_set_polar(2,4)
% points_s =
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
% > size(points_s)
% ans =
%      2     4
%
%See also
% PARTITION_OPTIONS, EQ_POINT_SET, EQ_REGIONS

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Function changed name from x2s2 to cart2polar2
% Optimize running time:
%   use slice assignments
%   trade space for time by using a cache
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
error(nargchk(2,4,nargin));
%
% dim is the number of dimensions
% N is the number of regions
%
% If the option 'offset' is 'extra', then use experimental extra offsets
% for S^2 and S^3 regions to try to minimize energy.
% The default is not to use extra offsets.
%
pdefault.extra_offset =  false;
if nargin < 3
    extra_offset = pdefault.extra_offset;
else
    popt = partition_options(pdefault, varargin{:});
    extra_offset = popt.extra_offset;
end
%
% Extra offset does not currently work for dim > 3,
% so quietly ignore this option in this case.
% Note that this also affects recursive calls to lower dimensions.
%
if dim > 3
    extra_offset = false;
end
%
% Check that dim and N are positive integers.
%
if ~( isnumeric(dim) && (dim >= 1) && (dim == floor(dim)) ) || ...
   ~( isnumeric(N) && (N >= 1) && (N == floor(N)) )
    fprintf('Usage: eq_point_set_polar(dim, N)\n');
    error('The arguments dim and N must be positive integers.');
end

if N == 1
    %
    % We have only one region, which must be the whole sphere.
    %
    points_s = zeros(dim,1);
    return;
end
%
% Start the partition of the sphere into N regions by partitioning
% to caps defined in the current dimension.
%
[a_cap, n_regions] = eq_caps(dim,N);
%
% a_cap is an increasing list of angles of the caps.
%
if dim == 1
    %
    % We have a circle and a_cap is an increasing list of angles of sectors,
    % with a_cap(k) being the cumulative arc length 2*pi/k.
    % The points are placed half way along each sector.
    %
    points_s = a_cap - pi/N;
    %
else
    %
    % We have a number of zones: two polar caps and a number of collars.
    % n_regions is the list of the number of regions in each zone.
    %
    n_collars = size(n_regions,2)-2;
    use_cache = dim >= 2;
    if use_cache
        cache_size = floor(n_collars/2);
        cache = cell(1,cache_size);
    end
    %
    % Start with the 'centre' point of the North polar cap.
    % This is the North pole.
    %
    points_s = zeros(dim,N);
    point_n = 2;
    %
    % Determine the 'centre' points for each collar.
    %
    if extra_offset && (dim == 3)
        R = eye(3);
    end
    if dim == 2
        offset = 0;
    end
    for collar_n = 1:n_collars
        %
        % a_top is the colatitude of the top of the current collar.
        %
        a_top = a_cap(collar_n);
        %
        % a_bot is the colatitude of the bottom of the current collar.
        %
        a_bot = a_cap(1+collar_n);
        %
        % n_in_collar is the number of regions in the current collar.
        %
        n_in_collar = n_regions(1+collar_n);
        %
        % The top and bottom of the collar are small (dim-1)-spheres,
        % which must be partitioned into n_in_collar regions.
        % Use eq_point_set_polar recursively to partition
        % the unit (dim-1)-sphere.
        % points_1 is the resulting list of points.
        %
        if use_cache
            twin_collar_n = n_collars-collar_n+1;
            if twin_collar_n <= cache_size && ...
                size(cache{twin_collar_n},2) == n_in_collar
                points_1 = cache{twin_collar_n};
            else
                points_1 = eq_point_set_polar(dim-1,n_in_collar,extra_offset);
                cache{collar_n} = points_1;
            end
        else
            points_1 = eq_point_set_polar(dim-1,n_in_collar,extra_offset);
        end
        %
        if extra_offset && (dim == 3) && (collar_n > 1)
            %
            % (Experimental)
            % Rotate 2-spheres to prevent alignment of north poles.
            %
            R = s2_offset(points_1)*R;
            points_1 = cart2polar2(R*polar2cart(points_1));
        end
        %
        % Given points_1, determine the 'centre' points for the collar.
        % Each point of points_1 is a 'centre' point on the (dim-1)-sphere.
        %
        % Angular 'centre' point;
        % The first angles are those of the current 'centre' point
        % of points_1, and the last angle in polar coordinates is the average of
        % the top and bottom angles of the collar,
        %
        a_point = (a_top+a_bot)/2;
        %
        point_1_n = 1:size(points_1,2);
        %
        if dim == 2
            %
            % The (dim-1)-sphere is a circle
            %
            points_s(1:dim-1,point_n+point_1_n-1) = mod(points_1(:,point_1_n)+2*pi*offset,2*pi);
            %
            % Given the number of sectors in the current collar and
            % in the next collar, calculate the next offset.
            % Accumulate the offset, and force it to be a number between 0 and 1.
            %
            offset = offset + circle_offset(n_in_collar,n_regions(2+collar_n),extra_offset);
            offset = offset - floor(offset);
        else
            points_s(1:dim-1,point_n+point_1_n-1) = points_1(:,point_1_n);
        end
        %
        points_s(dim, point_n+point_1_n-1) = a_point;
        point_n = point_n + size(points_1,2);
    end
    %
    % End with the 'centre' point of the bottom polar cap.
    %
    points_s(:,point_n) = zeros(dim,1);
    points_s(dim,point_n) = pi;
end
%
% end function

function [s_cap,n_regions] = eq_caps(dim,N)
%EQ_CAPS Partition a sphere into to nested spherical caps
%
%Syntax
% [s_cap,n_regions] = eq_caps(dim,N);
%
%Description
% [S_CAP,N_REGIONS] = EQ_CAPS(dim,N) does the following:
% 1) partitions the unit sphere S^dim into a list of spherical caps of
%    increasing colatitude and thus increasing area,
% 2) sets S_CAP to be an array of size (1 by N_COLLARS+2),
%    containing increasing colatitudes of caps, and
% 3) sets N_REGIONS to be an array of size (1 by N_COLLARS+2),
%    containing the intger number of regions in each corresponding zone of 
%    S^dim.
%
% The argument N is assumed to be a positive integer.
%
%Notes
% The value N_COLLARS is a positive integer function of dim and N.
%
% S_CAP[1] is C_POLAR, the colatitude of the North polar cap.
% S_CAP[N_COLLARS+1] is pi-C_POLAR.
% S_CAP[N_COLLARS+2] is pi.
%
% N_REGIONS[1] is 1.
% N_REGIONS[N_COLLARS+2] is 1.
% The sum of N_REGIONS is N.
%
%Examples
% > [s_cap,n_regions] = eq_caps(2,10)
% s_cap =
%     0.6435    1.5708    2.4981    3.1416
% n_regions =
%      1     4     4     1
%  
% > [s_cap,n_regions] = eq_caps(3,6)
% s_cap =
%     0.9845    2.1571    3.1416
% n_regions =
%      1     4     1
%
%See also
% EQ_REGIONS, EQ_POINT_SET_POLAR

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
error(nargchk(2,2,nargin));
error(nargoutchk(2,2,nargout));
%
% dim is the number of dimensions
% N is the number of regions
%
if dim == 1
    %
    % We have a circle. Return the angles of N equal sectors.
    %
    sector = 1:N;
    %
    % Make dim==1 consistent with dim>1 by
    % returning the longitude of a sector enclosing the
    % cumulative sum of arc lengths given by summing n_regions.
    %
    s_cap = sector*2*pi/N;
    n_regions = ones(size(sector));
    %
elseif N == 1
    %
    % We have only one region, which must be the whole sphere.
    %
    s_cap = [pi];
    n_regions = [1];
    %
else
    %
    % Given dim and N, determine c_polar,
    % the colatitude of the North polar spherical cap.
    %
    c_polar = polar_colat(dim,N);
    %
    % Given dim and N, determine the ideal angle for spherical collars.
    % Based on N, this ideal angle, and c_polar,
    % determine n_collars, the number of collars between the polar caps.
    %
    n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
    %
    % Given dim, N, c_polar and n_collars, determine r_regions,
    % a list of the ideal real number of regions in each collar,
    % plus the polar caps.
    % The number of elements is n_collars+2.
    % r_regions[1] is 1.
    % r_regions[n_collars+2] is 1.
    % The sum of r_regions is N.
    %
    r_regions = ideal_region_list(dim,N,c_polar,n_collars);
    %
    % Given N and r_regions, determine n_regions,
    % a list of the natural number of regions in each collar and
    % the polar caps.
    % This list is as close as possible to r_regions.
    % The number of elements is n_collars+2.
    % n_regions[1] is 1.
    % n_regions[n_collars+2] is 1.
    % The sum of n_regions is N.
    %
    n_regions = round_to_naturals(N,r_regions);
    %
    % Given dim, N, c_polar and n_regions, determine s_cap,
    % an increasing list of colatitudes of spherical caps which enclose the same area
    % as that given by the cumulative sum of regions.
    % The number of elements is n_collars+2.
    % s_cap[1] is c_polar.
    % s_cap[n_collars+1] is Pi-c_polar.
    % s_cap[n_collars+2] is Pi.
    %
    s_cap = cap_colats(dim,N,c_polar,n_regions);
    %
end
%
% end function

function c_polar = polar_colat(dim, N)
%POLAR_COLAT The colatitude of the North polar (top) spherical cap
%
% Given dim and N, determine the colatitude of the North polar spherical cap.
%
% c_polar = polar_colat(dim, N);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

enough = N > 2;
c_polar(N == 1) = pi;
c_polar(N == 2) = pi/2;
c_polar(enough) = sradius_of_cap(dim,area_of_ideal_region(dim,N(enough)));
%
% end function

function n_collars = num_collars(N,c_polar,a_ideal)
%NUM_COLLARS The number of collars between the polar caps
%
% Given N, an ideal angle, and c_polar,
% determine n_collars, the number of collars between the polar caps.
%
%  n_collars = num_collars(N,c_polar,a_ideal);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_collars = zeros(size(N));
enough = (N > 2) & (a_ideal > 0);
n_collars(enough) = max(1,round((pi-2*c_polar(enough))./a_ideal(enough)));
%
% end function


function r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%IDEAL_REGION_LIST The ideal real number of regions in each zone
%
% List the ideal real number of regions in each collar, plus the polar caps.
%
% Given dim, N, c_polar and n_collars, determine r_regions,
% a list of the ideal real number of regions in each collar,
% plus the polar caps.
% The number of elements is n_collars+2.
% r_regions[1] is 1.
% r_regions[n_collars+2] is 1.
% The sum of r_regions is N.
%
% r_regions = ideal_region_list(dim,N,c_polar,n_collars);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

r_regions = zeros(1,2+n_collars);
r_regions(1) = 1;
if n_collars > 0
    %
    % Based on n_collars and c_polar, determine a_fitting,
    % the collar angle such that n_collars collars fit between the polar caps.
    %
    a_fitting = (pi-2*c_polar)/n_collars;
    ideal_region_area = area_of_ideal_region(dim,N);    
    for collar_n = 1:n_collars
        ideal_collar_area = area_of_collar(dim, c_polar+(collar_n-1)*a_fitting, c_polar+collar_n*a_fitting);
        r_regions(1+collar_n) = ideal_collar_area / ideal_region_area;
    end
end
r_regions(2+n_collars) = 1;

% end function

function area = area_of_ideal_region(dim,N)
%AREA_OF_IDEAL_REGION Area of one region of an EQ partition
%
%Syntax
% area = area_of_ideal_region(dim,N);
%
%Description
% AREA = AREA_OF_IDEAL_REGION(dim,N) sets AREA to be the area of one of N equal
% area regions on S^dim, that is 1/N times AREA_OF_SPHERE(dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result AREA will be an array of the same size as N.
%
%Examples
% > area=area_of_ideal_region(3,1:6)
% area =
%    19.7392    9.8696    6.5797    4.9348    3.9478    3.2899
%
%See also
% AREA_OF_SPHERE

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

area = area_of_sphere(dim)./N;
%
% end function

function area = area_of_sphere(dim)
%AREA_OF_SPHERE Area of sphere
%
%Syntax
% area = area_of_sphere(dim);
%
%Description
% AREA = AREA_OF_SPHERE(dim) sets AREA to be the area of the sphere S^dim,
%
% The argument dim must be a positive integer or an array of positive integers.
% The result AREA will be an array of the same size as dim.
%
%Notes
% The area of S^dim is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% Ref: [Mue98] p39.
%
%Examples
% > area=area_of_sphere(1:7)
% area =
%     6.2832   12.5664   19.7392   26.3189   31.0063   33.0734   32.4697
%
%See also
% AREA_OF_CAP, VOLUME_OF_BALL

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

power = (dim+1)/2;
area = (2*pi.^power./gamma(power));
%
% end function

function s_cap = sradius_of_cap(dim, area)
%SRADIUS_OF_CAP Spherical radius of spherical cap of given area
%
%Syntax
% s_cap = sradius_of_cap(dim, area);
%
%Description
% S_CAP = SRADIUS_OF_CAP(dim, AREA) sets S_CAP to be the spherical radius of 
% an S^dim spherical cap of area AREA.
%
% The argument dim must be a positive integer.
% The argument AREA must be a real number or an array of real numbers.
% The result S_CAP will be an array of the same size as AREA.
%
%Notes
% S_CAP is assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% For dim <= 2, S_CAP is calculated in closed form.
% Otherwise, S_CAP is approximated using the Matlab function FZERO.
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
% > s_cap=sradius_of_cap(2,area_of_sphere(2)/2)
%
% s_cap =
%     1.5708
%
% > s_cap=sradius_of_cap(3,(0:4)*area_of_sphere(3)/4)
% s_cap =
%          0    1.1549    1.5708    1.9867    3.1416
%
%See also
% FZERO, AREA_OF_CAP

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.02 $ $Date 2005-04-25 $
% Use asin rather than acos to avoid subtraction.
% Use symmetry to avoid loss of accuracy in the Southern hemisphere.
% Remove check for when area is close to area of sphere.
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

switch dim
case 1
    s_cap = area/2;
case 2
    s_cap = 2*asin(sqrt(area/pi)/2);
otherwise
    %
    % Flatten area into a row vector.
    %
    shape = size(area);
    n = prod(shape);
    area = reshape(area,1,n);
    s_cap = zeros(size(area));
    for k = 1:n
        ak = area(k);
        as = area_of_sphere(dim);
        if ak >= as
            s_cap(k) = pi;
        else
            if (2*ak > as)
                ak = as - ak;
                flipped = true;
            else
                flipped = false;
            end
            sk = ...
                fzero(inline(sprintf('area_of_cap(%d,s)-%21.14g',dim,ak),'s'),[0,pi]);
            if flipped
                s_cap(k) = pi - sk;
            else
                s_cap(k) = sk;
            end
        end
    end
    %
    % Reshape output to same array size as original area.
    %
    s_cap = reshape(s_cap,shape);
end
% end function

function angle = ideal_collar_angle(dim,N)
%IDEAL_COLLAR_ANGLE The ideal angle for spherical collars of an EQ partition
%
%Syntax
% angle = ideal_collar_angle(dim,N);
%
%Description
% ANGLE = IDEAL_COLLAR_ANGLE(dim,N) sets ANGLE to the ideal angle for the 
% spherical collars of an EQ partition of the unit sphere S^dim into N regions.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result ANGLE will be an array of the same size as N.
%
%Notes
% The ideal collar angle is determined by the side of a dim-dimensional
% hypercube of the same volume as the area of a single region of an N region
% equal area partition of S^dim.
%
% Since the EQ partition for N < 3 has no spherical collars, 
% the recursive zonal equal area sphere partitioning algorithm does not use 
% ideal_collar_angle(dim,N) for N < 3.
%
%Examples
% > angle = ideal_collar_angle(2,10)
%  angle =
%      1.1210
%  
% > angle = ideal_collar_angle(3,1:6)
%  angle =
%      2.7026    2.1450    1.8739    1.7025    1.5805    1.4873
%
%See also
% AREA_OF_IDEAL_REGION

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

angle = area_of_ideal_region(dim,N).^(1/dim);
%
% end function

function area = area_of_collar(dim, a_top, a_bot)
%AREA_OF_COLLAR Area of spherical collar
%
%Syntax
% area = area_of_collar(dim, a_top, a_bot);
%
%Description
% AREA = AREA_OF_COLLAR(dim, A_TOP, A_BOT) sets AREA to be the area of 
% an S^dim spherical collar specified by A_TOP, A_BOT, where
% A_TOP is top (smaller) spherical radius,
% A_BOT is bottom (larger) spherical radius.
%
% The argument dim must be a positive integer.
% The arguments A_TOP and A_BOT must be real numbers or arrays of real numbers,
% with the same array size.
% The result AREA will be an array of the same size as A_TOP.
%
%Notes
% A_TOP and A_BOT are assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
% > area=area_of_collar(2,0:2,1:3)
% area =
%     2.8884    6.0095    3.6056
%
%See also
% AREA_OF_CAP, AREA_OF_SPHERE

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

area = area_of_cap(dim, a_bot) - area_of_cap(dim, a_top);
%
% end function

function area = area_of_cap(dim, s_cap)
%AREA_OF_CAP Area of spherical cap
%
%Syntax
% area = area_of_cap(dim, s_cap);
%
%Description
% AREA = AREA_OF_CAP(dim, S_CAP) sets AREA to be the area of an S^dim spherical
% cap of spherical radius S_CAP.
%
% The argument dim must be a positive integer.
% The argument S_CAP must be a real number or an array of real numbers.
% The result AREA will be an array of the same size as S_CAP.
%
%Notes
% S_CAP is assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% For dim <= 2, and for dim==3 (when pi/6 <= s_cap <= pi*5/6),
% AREA is calculated in closed form, using the analytic solution of
% the definite integral given in the reference.
% Otherwise, AREA is calculated using the Matlab function BETAINC,
% the incomplete Beta function ratio.
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
% > a=area_of_cap(2,pi/2)
% a =
%     6.2832
% > a=area_of_cap(3,0:pi/4:pi)
% a =
%          0    1.7932    9.8696   17.9460   19.7392
%
%See also
% SRADIUS_OF_CAP

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.02 $ $Date 2005-04-24 $ PL
% Use incomplete Beta function BETAINC for dim == 3,
% (when s_cap < pi/6 or s_cap > pi*5/6) and for all dim > 3.
% Use sin(s_cap).^2 in preference to (1-cos(s_cap))/2.
% $Revision 1.01 $ $Date 2005-03-16 $ PL
% Use incomplete Beta function BETAINC for dim > 8.
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

switch dim
case 1
    area = 2 * s_cap;
case 2
    area = 4*pi * sin(s_cap/2).^2;
case 3
     %
     % Flatten s_cap into a row vector.
     %
     shape = size(s_cap);
     n = prod(shape);
     s_cap = reshape(s_cap,1,n);
     area = zeros(size(s_cap));
     %
     % Near the poles, use the incomplete Beta function ratio.
     %
     pole = (s_cap < pi/6) | (s_cap > pi*5/6);
     area(pole) = area_of_sphere(dim) * betainc(sin(s_cap(pole)/2).^2,dim/2,dim/2);
     %
     % In the tropics, use closed solution to integral.
     %
     trop = s_cap(~pole);
     area(~pole) = (2*trop-sin(2*trop))*pi;

     area = reshape(area,shape);
otherwise
     area = area_of_sphere(dim) * betainc(sin(s_cap/2).^2,dim/2,dim/2);
end
%
% end function

function n_regions = round_to_naturals(N,r_regions)
%ROUND_TO_NATURALS Round off a given list of numbers of regions
%
% Given N and r_regions, determine n_regions,
% a list of the natural number of regions in each collar and the polar caps.
% This list is as close as possible to r_regions, using rounding.
% The number of elements is n_collars+2.
% n_regions[1] is 1.
% n_regions[n_collars+2] is 1.
% The sum of n_regions is N.
%
% n_regions = round_to_naturals(N,r_regions);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_regions = r_regions;
discrepancy = 0;
for zone_n = 1:size(r_regions,2)
    n_regions(zone_n) = round(r_regions(zone_n)+discrepancy);
    discrepancy = discrepancy+r_regions(zone_n)-n_regions(zone_n);
end
%
% end function

function c_caps = cap_colats(dim,N,c_polar,n_regions)
%CAP_COLATS Colatitudes of spherical caps enclosing cumulative sum of regions
%
% Given dim, N, c_polar and n_regions, determine c_caps,
% an increasing list of colatitudes of spherical caps which enclose the same area
% as that given by the cumulative sum of regions.
% The number of elements is n_collars+2.
% c_caps[1] is c_polar.
% c_caps[n_collars+1] is Pi-c_polar.
% c_caps[n_collars+2] is Pi.
%
% c_caps = cap_colats(dim,N,c_polar,n_regions);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

c_caps = zeros(size(n_regions));
c_caps(1) = c_polar;
ideal_region_area = area_of_ideal_region(dim,N);
n_collars = size(n_regions,2)-2;
subtotal_n_regions = 1;
for collar_n = 1:n_collars
    subtotal_n_regions = subtotal_n_regions+n_regions(1+collar_n);
    c_caps(collar_n+1) =sradius_of_cap(dim,subtotal_n_regions*ideal_region_area);
end
c_caps(1+n_collars+1) = pi;

% end function

function popt = partition_options(pdefault, varargin)
%PARTITION_OPTIONS Options for EQ partition
%
%Syntax
% popt = partition_options(pdefault,options);
%
%Description
% POPT = PARTITION_OPTIONS(PDEFAULT,options) collects partition options,
% specified as name, value pairs, and places these into the structure POPT.
% The structure PDEFAULT is used to define default option values.
%
% The structures pdefault and popt may contain the following fields:
% extra_offset:  boolean
%
% The following partition options are available.
%
% 'offset':      Control extra rotation offsets for S^2 and S^3 regions.
%     'extra':   Use extra rotation offsets for S^2 and S^3 regions, to try
%                to minimize energy.
%                Sets opt.extra_offset to true.
%     'normal':  Do not use extra offsets
%                Sets opt.extra_offset to false.
%
% Some shortcuts are also provided.
% POPT = PARTITION_OPTIONS(pdefault) just sets POPT to PDEFAULT.
%
% The following are equivalent to PARTITION_OPTIONS(PDEFAULT,'offset','extra'):
% PARTITION_OPTIONS(PDEFAULT,true)
% PARTITION_OPTIONS(PDEFAULT,'extra')
%
% The following are equivalent to PARTITION_OPTIONS(PDEFAULT,'offset','normal'):
% PARTITION_OPTIONS(PDEFAULT,false)
% PARTITION_OPTIONS(PDEFAULT,'normal')
%
%Examples
% > pdefault.extra_offset=false;
% > popt=partition_options(pdefault,'offset','extra')
% popt =
%     extra_offset: 1
%
% > popt=partition_options(pdefault,false)
% popt =
%     extra_offset: 0

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

popt = pdefault;
nargs = length(varargin);

if nargs == 1
    %
    % Short circuit: single argument is value of extra_offset
    %
    value = varargin{1};
    switch value
    case true
        popt.extra_offset = true;
    case false
        popt.extra_offset = false;
    case 'extra'
        popt.extra_offset = true;
    case 'normal'
        popt.extra_offset = false;
    otherwise
        value_error(value,varargin{:});
    end
    return;
end

nopts = floor(nargs/2);
opt_args = {varargin{1:2:2*nopts-1}};
for k=1:nopts
    if ~ischar([opt_args{k}])
        fprintf('Option names must be character strings\n');
        option_error(varargin{:});
    end
end    
opt_vals = {varargin{2:2:2*nopts}};

option_name = 'offset';
pos = strmatch(option_name,opt_args,'exact');
if ~isempty(pos)
    if (length(pos) == 1)
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'extra'
        popt.extra_offset = true;
    case 'normal'
        popt.extra_offset = false;
    otherwise
        value_error(value,varargin{:});
    end
end

function duplicate_error(option_name,varargin)
fprintf('Duplicate option %s\n',option_name);
option_error(varargin{:});
%
% end function

function value_error(value,varargin)
fprintf('Invalid option value ');
disp(value);
option_error(varargin{:});
%
% end function

function option_error(varargin)
fprintf('Error in options:\n');
disp(varargin);
error('Please check "help partition_options" for options');
%
% end function

function offset = circle_offset(n_top,n_bot,extra_twist)
%CIRCLE_OFFSET Try to maximize minimum distance of center points for S^2 collars
%
% Given n_top and n_bot, calculate an offset.
%
% The values n_top and n_bot represent the numbers of 
% equally spaced points on two overlapping circles.
% The offset is given in multiples of whole rotations, and
% consists of three parts;
% 1) Half the difference between a twist of one sector on each of bottom and top. 
% This brings the centre points into alignment.
% 2) A rotation which will maximize the minimum angle between
% points on the two circles.
% 3) An optional extra twist by a whole number of sectors on the second circle.
% The extra twist is added so that the location of
% the minimum angle  between circles will
% progressively twist around the sphere with each collar.
%
% offset = circle_offset(n_top,n_bot,extra_twist);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if nargin < 3
    extra_twist = false;
end
offset = (1/n_bot - 1/n_top)/2 + gcd(n_top,n_bot)/(2*n_top*n_bot);
if extra_twist
    twist = 6;
    offset = offset + twist/n_bot;
end
% end function

function x = polar2cart(s)
%POLAR2CART Convert spherical polar to Cartesian coordinates
%
%Syntax
% x = polar2cart(s);
%
%Description
% X = POLAR2CART(S) sets X to be the Cartesian coordinates of the points
% represented by the spherical polar coordinates S.
%
% S is assumed to be an array of size (dim by N) representing N points of
% S^dim in spherical polar coordinates, where dim and N are positive integers.
% N will be an array of size (dim+1 by N).
%
%Examples
% > s
% s =
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
% > x=polar2cart(s)
% x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
%See also
% CART2POLAR2

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Change name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(s,1);
n=size(s,2);
x = zeros(dim+1,n);
sinprod  = 1;
for k = dim:-1:2
    x(k+1,:) = sinprod.*cos(s(k,:));
    sinprod  = sinprod.*sin(s(k,:));
end
x(2,:)=sinprod.*sin(s(1,:));
x(1,:)=sinprod.*cos(s(1,:));
r = sqrt(sum(x.^2));
mask = (r ~= 1);
if  size(r(mask),2) > 0
    x(:,mask) = x(:,mask)./(ones(dim+1,1)*r(mask));
end
%
% end function
"""

if __name__== "__main__":
	main()
