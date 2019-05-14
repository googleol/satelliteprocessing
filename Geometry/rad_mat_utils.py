#------------------------------------------------------------------------------------
#-------- CONFIGURATION ENDS : NOT TO EDIT BELOW HERE FOR NOMINAL EDIT ----------------

import numpy as np
from math import sin,cos,radians
#import nvector as nv
import pymap3d as pm
import pyproj
from datetime import datetime,timedelta,timezone
from osgeo import gdal, osr
#from indexGeoTiff import filesinsidefolder
from os import path
import fiona
from shapely.geometry import Point,Polygon,MultiPoint,shape
#from os import walk,path


#FUNCTIONS -- :



def NORM3(v3):
    len_vec=(v3[0])**2+(v3[1])**2+(v3[2])**2
    return v3/(len_vec)**0.5
def Rotate_X(cord,t):
    #t=radians(t)
    x=cord[0]
    y=cord[1]
    z=cord[2]
    xr=x
    yr=y*cos(t)+z*sin(t)
    zr=-y*sin(t)+z*cos(t)

    return (xr,yr,zr)
def Rotate_Y(cord,t):
    #t=radians(t)
    x=cord[0]
    y=cord[1]
    z=cord[2]
    yr=y
    xr=x*cos(t)+z*sin(t)
    zr=-x*sin(t)+z*cos(t)

    return (xr,yr,zr)
def Rotate_Z(cord,t):
    #t=radians(t)
    x=cord[0]
    y=cord[1]
    z=cord[2]
    zr=z
    xr=x*cos(t)+y*sin(t)
    yr=-x*sin(t)+y*cos(t)

    return (xr,yr,zr)

def Rotate_Vector(order,vector,angles,angles_order): #len order and angles must be equal
    for rot_index in range(0,len(order)):
        if(order[rot_index]=='X'):
            vector=Rotate_X(vector,angles[angles_order[rot_index]])
        elif(order[rot_index]=='Y'):
            vector=Rotate_Y(vector,angles[angles_order[rot_index]])
        elif(order[rot_index]=='Z'):
            vector=Rotate_Z(vector,angles[angles_order[rot_index]])
    return vector

def eulerAnglesToRotationMatrix(theta) :

    R_x = np.array([[1,         0,                  0                   ],
                    [0,         cos(theta[0]), -sin(theta[0]) ],
                    [0,         sin(theta[0]), cos(theta[0])  ]
                    ])



    R_y = np.array([[cos(theta[1]),    0,      sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-sin(theta[1]),   0,      cos(theta[1])  ]
                    ])

    R_z = np.array([[cos(theta[2]),    -sin(theta[2]),    0],
                    [sin(theta[2]),    cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])


    R = np.dot(R_y, np.dot(  R_z,R_x))

    return R




def eulerAnglesToRotationMatrixZXZ(theta) :

    R_x = np.array([[1,         0,                  0                   ],
                    [0,         cos(theta[1]), -sin(theta[1]) ],
                    [0,         sin(theta[1]), cos(theta[1])  ]
                    ])


    #
    # R_y = np.array([[cos(theta[1]),    0,      sin(theta[1])  ],
    #                 [0,                     1,      0                   ],
    #                 [-sin(theta[1]),   0,      cos(theta[1])  ]
    #                 ])

    R_z = np.array([[cos(theta[2]),    -sin(theta[2]),    0],
                    [sin(theta[2]),    cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])
    R_zra = np.array([[cos(theta[0]),    -sin(theta[0]),    0],
                    [sin(theta[0]),    cos(theta[0]),     0],
                    [0,                     0,                      1]
                    ])


    R = np.dot(R_z, np.dot( R_x, R_zra ))

    return R

def ORBel_frm_SV(svr,svv):
    from astropy import units as u
    from poliastro.bodies import Earth
    from poliastro.twobody import Orbit
    from poliastro.core.elements import rv2coe

    state_vector_r=svr
    state_vector_v=svv

    # unit vector conversion


    r = np.array(state_vector_r) * u.km
    v = np.array(state_vector_v) * u.km / u.s

    ss = Orbit.from_vectors(Earth, r, v)

    #standard_grav_param_earth
    k=3.986004418*100000
    coe=rv2coe(k,r,v)
    #a,e,i,raom,apom,truean=coe[0,1,2,3,4,5]

    # print('\n')
    # print('---------- ORBIT CALCULATIONS ---------------------- \n\n')
    # print(' SEMI MAJOR AXIS IS -- : '+str(coe[0])+' Meters \n')
    # print(' ECCENTRICITY IS -- : '+str(coe[1])+' \n')
    # print(' INCLINATION IS -- : '+str(coe[2])+' radians \n') #unit confirmation?
    # print(' ARGUMENT OF RA IS -- : '+str(coe[3])+' radians \n')
    # print(' ARGUMENT OF PERIGEE IS -- : '+str(coe[4])+' radians \n')
    # print(' TRUE ANAMOLY IS -- : '+str(coe[5])+' radians \n')
    # print('------------------------------------')
    # print('\n')

    return coe
def SUBPOINT(point):
    #-----WGS 84 --------------
    dat_a=6378137 #Meters
    dat_f=1/298.257223563 #flat scale for oblateness
    dat_ep=1/((1-dat_f)*(1-dat_f))

    ray_dir=point

    xd=ray_dir[0]#-31.46071792#-
    yd=ray_dir[1]#58.59611618#-
    zd=ray_dir[2]#27.47631664#-

    xc=point[0]
    yc=point[1]
    zc=point[2]

    print(xc,yc,zc)


    # [xd,yd,zd]=pm.eci2ecef(NORM3([xd,yd,zd]), t)
    # [xc,yc,zc]=pm.eci2ecef([xc,yc,zc], t)




    qa=xd**2+yd**2+(zd**2)*dat_ep
    qb=2*(xc*xd+yc*yd+zc*zd*dat_ep)
    qc=xc**2+yc**2+((zc**2)*dat_ep)-(dat_a**2)

    #print('\n')
    #print('---------- DATUM CALCULATIONS ---------------------- \n\n')


    # print('\n-----qa  :  '+str(qa)+'---\n')
    # print('\n-----qb :  '+str(qb)+'---\n')
    # print('\n-----qc  :  '+str(qc)+'---\n')
    # print('\n-----qb**2  :  '+str(qb*qb)+'---\n')
    # print('\n-----4*qa*qc  :  '+str(4*qa*qc)+'---\n')

    #print ('---------- DATUM DETERMINANT IS --------- '+str(qb*qb-4*qa*qc)+'\n')
    ray_roots = np.roots([qa,qb,qc])
    #print (ray_roots)
    earth_point=[]
    earth_point.append([xc+ray_roots[0]*xd,yc+ray_roots[0]*yd,zc+ray_roots[0]*zd])
    earth_point.append([zc+ray_roots[1]*xd,yc+ray_roots[1]*yd,zc+ray_roots[1]*zd])

    #print(earth_point)

    #print('--NORMS '+str(np.linalg.norm(earth_point[0]))+'--:--'+str(np.linalg.norm(earth_point[1]))+'\n')

    #print('\n\n ---------- STATE VECTOR 2 NORMAL PROJECTION -------------- \n')
    return earth_point[0]

def RAY_INTERSECT_WGS84_PYGAC(point,ray_dir,t):

    centre = -point
    a__ = 6378137  # km
    # b__ = 6356.75231414 # km, GRS80
    b__ = 6356752.314245  # km, WGS84
    radius = np.array([[1 / a__, 1 / a__, 1 / b__]]).T
    shape = ray_dir.shape

    xr_ = ray_dir.reshape([3, -1])*radius
    print(xr_)
    cr_ = centre.reshape([3, -1])*radius
    print(cr_)
    ldotc = np.einsum("ij,ij->j", xr_, cr_)
    print(ldotc)
    lsq = np.einsum("ij,ij->j", xr_, xr_)
    print(lsq)
    csq = np.einsum("ij,ij->j", cr_, cr_)
    print(csq)
    #print(ldotc**2-(csq*lsq+lsq))
    #print((ldotc**2-csq*lsq+lsq)**0.5)
    d1_ = (ldotc - np.sqrt(ldotc**2-csq*lsq+lsq))/ lsq
    #print(d1_)
    #print(ray_dir * d1_.reshape(shape[1:]) - centre)
    # return the actual pixel positions
    earth_point=[]
    earth_point.append(ray_dir * d1_.reshape(shape[1:]) - centre)
    earth_point[0]=pm.eci2ecef((earth_point[0]),t)
    #earth_point[0]=6371e3*np.array(earth_point[0])
    # 6371e3 * np.vstack(ecef)
    ecx=-1*earth_point[0][0]
    ecy=-1*earth_point[0][1]
    ecz=-1*earth_point[0][2]

    llh=pm.ecef2geodetic(ecx,ecy,ecz)
    #print((ecx**2+ecy**2+ecz**2)**0.5)
    #print("LAT(deg) LON(deg) HEIGHT(m) :"+str(pm.ecef2geodetic(ecx,ecy,ecz)))
    #print("LAT(deg) LON(deg) HEIGHT(m) :"+str(ECI_TO_WGS84LLH(-1*np.array(earth_point[0]),t  )))
    #print("LAT(deg) LON(deg) HEIGHT(m) :"+str(ECEF_TO_WGS84LLH(-1*np.array(earth_point[0]))))
    #earth_point[0]=np.array([ecx,ecy,ecz))/1000
    return llh

def RAY_INTERSECT_WGS84(point,ray_dir,t):
    #-----WGS 84 --------------
    dat_a=6378137 #Meters
    dat_f=1/298.257223563 #flat scale for oblateness
    dat_ep=1/((1-dat_f)*(1-dat_f))


    xd=ray_dir[0]#-31.46071792#-
    yd=ray_dir[1]#58.59611618#-
    zd=ray_dir[2]#27.47631664#-

    xc=point[0]
    yc=point[1]
    zc=point[2]


    qa=xd**2+yd**2+(zd**2)*dat_ep
    qb=2*(xc*xd+yc*yd+zc*zd*dat_ep)
    qc=xc**2+yc**2+((zc**2)*dat_ep)-(dat_a**2)

    ray_roots = np.roots([qa,qb,qc])
    earth_point=[]
    earth_point.append([xc+ray_roots[0]*xd,yc+ray_roots[0]*yd,zc+ray_roots[0]*zd])
    earth_point.append([xc+ray_roots[1]*xd,yc+ray_roots[1]*yd,zc+ray_roots[1]*zd])

    #print(earth_point)
    #6371e3 * np.vstack(ecef)
    #print('--NORMS '+str(np.linalg.norm(earth_point[0]))+'--:--'+str(np.linalg.norm(earth_point[1]))+'\n')

    #print('\n\n ---------- STATE VECTOR 2 NORMAL PROJECTION -------------- \n')
    #print(earth_point)
    #print(6371e3*np.vstack(earth_point[0]))
    earth_point[0]=pm.eci2ecef((earth_point[1]),t)
    #earth_point[0]=6371e3*np.array(earth_point[0])
    # 6371e3 * np.vstack(ecef)
    ecx=-1*earth_point[0][0]
    ecy=-1*earth_point[0][1]
    ecz=-1*earth_point[0][2]
    #print((ecx**2+ecy**2+ecz**2)**0.5)
    #print("LAT(deg) LON(deg) HEIGHT(m) :"+str(pm.ecef2geodetic(ecx,ecy,ecz)))
    #print("LAT(deg) LON(deg) HEIGHT(m) :"+str(ECI_TO_WGS84LLH(-1*np.array(earth_point[0]),t  )))
    #print("LAT(deg) LON(deg) HEIGHT(m) :"+str(ECEF_TO_WGS84LLH(-1*np.array(earth_point[0]))))
    #earth_point[0]=np.array([ecx,ecy,ecz))/1000
    #llh=pm.ecef2geodetic(ecx,ecy,ecz)
    # ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    # lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    # lon, lat, alt = pyproj.transform(ecef, lla, ecx, ecy, ecz, radians=False)
    llh=pm.ecef2geodetic(ecx,ecy,ecz)
    #print("LAT(deg) LON(deg) HEIGHT(m) :"+str(llh))
    #print('---------- ------------------------------- -------------- \n')
    #print(llh)
    return(llh)
def MISSION_TIME_TO_DATETIME_LS(DATE,TIME):
    #FOR RS2-LISS3 DO LIGHT TIME SUBTRACTION
    ls=300000000
    height=817000#rs2-l3

    offset=4963
    res=23.5
    sectaseconds = TIME.split('.')[1]
    microseconds=(int(sectaseconds)/10)
    pro_time=TIME.split('.')[0]+'.'+str(int(microseconds))
    pt =datetime.strptime(DATE+' '+pro_time,'%Y-%m-%d %H:%M:%S.%f')
    fin_time=pt-timedelta(seconds=height/ls)
    #fin_time=pt-timedelta(seconds=(offset/res)*0.003319466666)
    return fin_time
def MISSION_TIME_TO_DATETIME(DATE,TIME):
    #FOR RS2-LISS3 DO LIGHT TIME SUBTRACTION
    ls=300000000
    height=817000#rs2-l3
    if(type(DATE)!=type('hi') or type(TIME)!=type('hi') ):
        DATE=DATE.decode('ascii')
        TIME=TIME.decode('ascii')
    #print(DATE,TIME)
    sectaseconds = TIME.split('.')[1]
    microseconds=(int(sectaseconds)/10)
    pro_time=TIME.split('.')[0]+'.'+str(int(microseconds))
    pt =datetime.strptime(DATE+' '+pro_time,'%Y-%m-%d %H:%M:%S.%f')
    fin_time=pt-timedelta(seconds=height/ls)
    return pt
def transform_wgs84_to_utm(lon, lat):
    def get_utm_zone(longitude):
        return (int(1+(longitude+180.0)/6.0))

    def is_northern(latitude):
        """
        Determines if given latitude is a northern for UTM
        """
        if (latitude < 0.0):
            return 0
        else:
            return 1
    utm_coordinate_system = osr.SpatialReference()
    utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetUTM(get_utm_zone(lon), is_northern(lat))

    wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system

    # create transform component
    wgs84_to_utm_transform = osr.CoordinateTransformation(wgs84_coordinate_system, utm_coordinate_system) # (<from>, <to>)
    return [wgs84_to_utm_transform.TransformPoint(lon, lat, 0),utm_coordinate_system.ExportToWkt()]
def Height_From_DEM(lon,lat,DEM_FILES_DIRECTORY,nom_height=0):
    #form=['tiff','tif','.img']#DEM FORMATS
    #fNames=filesinsidefolder(DEM_FILES_DIRECTORY,form)
    #print('HFD')
    indexFile=path.join(DEM_FILES_DIRECTORY,'index.shp')
    # driver= ogr.GetDriverByName('ESRI Shapefile')
    # indexLayerSource = driver.Open(indexFile, 0) # import fiona
    #shape = fiona.open("my_shapefile.shp")0 means read-only. 1 means writeable.
    # indexLayer=indexLayerSource.GetLayer()
    intersects=0
    for feature in fiona.open(indexFile):
        #first = shape_f.next()
        polygon=shape(feature['geometry'])

        point=Point(lon,lat)
        #polygon = Polygon(pts)
        #vertices=list(polygon.exterior.coords)

        if(point.intersects(polygon)):
            intersects=1
            raster_file=feature['properties']['fileName']
            #print(raster_file)
            dataset=gdal.Open(raster_file)
            cols = dataset.RasterXSize
            rows = dataset.RasterYSize

            transform = dataset.GetGeoTransform()

            xOrigin = transform[0]
            yOrigin = transform[3]
            pixelWidth = transform[1]
            pixelHeight = -transform[5]
            band = dataset.GetRasterBand(1)
            data = band.ReadAsArray(0, 0, cols, rows)
            col = int((point.x - xOrigin) / pixelWidth)
            row = int((yOrigin - point.y ) / pixelHeight)
            height=data[row][col]
            return((1,(lat,lon,height)))

    if(intersects==0):
        #print('NO INTERSECTION')
        return((0,(lat,lon,nom_height)))
def Max_Local_Height_From_DEM(lon,lat,DEM_FILES_DIRECTORY,nom_height=0,scale=(3,3)):
    #form=['tiff','tif','.img']#DEM FORMATS
    #fNames=filesinsidefolder(DEM_FILES_DIRECTORY,form)
    indexFile=path.join(DEM_FILES_DIRECTORY,'index.shp')
    # driver= ogr.GetDriverByName('ESRI Shapefile')
    # indexLayerSource = driver.Open(indexFile, 0) # import fiona
    #shape = fiona.open("my_shapefile.shp")0 means read-only. 1 means writeable.
    # indexLayer=indexLayerSource.GetLayer()
    intersects=0
    for feature in fiona.open(indexFile):
        #first = shape_f.next()
        polygon=shape(feature['geometry'])

        point=Point(lon,lat)
        #polygon = Polygon(pts)
        #vertices=list(polygon.exterior.coords)

        if(point.intersects(polygon)):
            intersects=1
            raster_file=feature['properties']['fileName']
            #print(raster_file)
            dataset=gdal.Open(raster_file)
            cols = dataset.RasterXSize
            rows = dataset.RasterYSize

            transform = dataset.GetGeoTransform()

            xOrigin = transform[0]
            yOrigin = transform[3]
            pixelWidth = transform[1]
            pixelHeight = -transform[5]
            band = dataset.GetRasterBand(1)
            data = band.ReadAsArray(0, 0, cols, rows)
            col = int((point.x - xOrigin) / pixelWidth)
            row = int((yOrigin - point.y ) / pixelHeight)
            max_height=0
            mh_index=(col,row)
            for x in range(col-scale[0],col+scale[0]):
                for y in range(row-scale[1],row+scale[1]):
                    if not(x<0 or x>=cols or y<0 or y>=rows):
                        if(data[y][x]>max_height):
                            max_height=data[y][x]
                            mh_index=(x,y)
            mh_lon=mh_index[0]*pixelWidth+xOrigin
            mh_lat=mh_index[1]*pixelHeight+yOrigin



            from scipy.ndimage import maximum_filter
            maxs = maximum_filter(data, size=(20,20))
            height=maxs[row][col]

            return((1,(mh_lat,mh_lon,max_height)))

    if(intersects==0):
        #print('NO INTERSECTION')
        return((0,(lat,lon,nom_height)))









def RAY_INTERSECT_WGS84_TERRAIN(point,ray_dir,t,DEM_FILE_LOCATION=None,nom_height=0):
    lat,lon,height=RAY_INTERSECT_WGS84(point,ray_dir,t)
    s,height_min=Height_From_DEM(lon,lat,DEM_FILE_LOCATION)
    if(s==0):
        return RAY_INTERSECT_WGS84(point,ray_dir,t)
    #print(Max_Local_Height_From_DEM(lon,lat,DEM_FILE_LOCATION,scale=(5,5)))
    s,(max_lat,max_lon,max_height)=Max_Local_Height_From_DEM(lon,lat,DEM_FILE_LOCATION,scale=(5,5))
    if(s==0):
        return RAY_INTERSECT_WGS84(point,ray_dir,t)

    eci_ter=pm.ecef2eci(pm.geodetic2ecef(lat,lon,height), t)
    #print([(max_lon,max_lat),max_height, t])
    eci_ter_lmax=pm.ecef2eci(pm.geodetic2ecef(max_lat,max_lon,max_height),t)

    #print(mid_llh[2],height_mid


    xd=ray_dir[0]#-31.46071792#-
    yd=ray_dir[1]#58.59611618#-
    zd=ray_dir[2]#27.47631664#-

    xc=-point[0]
    yc=-point[1]
    zc=-point[2]

    xn=eci_ter_lmax[0]
    yn=eci_ter_lmax[1]
    zn=eci_ter_lmax[2]

    #print(xn,yn,zn)

    t_par=(((xn)**2+(yn)**2+(zn)**2)-(xc*xn)-(yc*yn)-(zc*zn))/(xd*xn+yd*yn+zd*zn)
    rez=np.array([xc+t_par*xd,yc+t_par*yd,zc+t_par*zd])
    #print(pm.eci2geodetic((rez),t))

    #print(eci_ter)
    height_delta=9999999999
    eci_ter_mid=(eci_ter+rez)/2
    height_delta_old=9898989898
    while(abs(height_delta)>5 and height_delta_old!=height_delta):
        #print(str(height_delta)+': height delta')
        mid_llh=pm.eci2geodetic(eci_ter_mid,t)#latlonheight
        s,height_mid=Height_From_DEM(mid_llh[1],mid_llh[0],DEM_FILE_LOCATION)
        if(s==0):
            return RAY_INTERSECT_WGS84(point,ray_dir,t)
        #print(mid_llh[2],height_mid[2],mid_llh[2]-height_mid[2])
        height_delta_old=height_delta
        height_delta=mid_llh[2]-height_mid[2]

        if(height_delta<0):
            eci_ter_mid_new=(rez+eci_ter_mid)/2
            eci_ter=eci_ter_mid
            eci_ter_mid=eci_ter_mid_new
        else:
            eci_ter_mid_new=(eci_ter_mid+eci_ter)/2
            rez=eci_ter_mid
            eci_ter_mid=eci_ter_mid_new
    #print(mid_llh)
    return(mid_llh)
#def Interpolate_State_Vector(svb,sva,t):
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total:
        print()
