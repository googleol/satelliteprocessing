from loguru import logger
import numpy as np
from math import sin,cos,radians,atan,degrees,acos
import pymap3d as pm
from .rad_mat_utils import *
import h5py
import gdal,ogr,osr
import matplotlib.pyplot as plt
import json
from pprint import pprint
from datetime import datetime,timedelta,timezone
import bisect
from .gdal_warp_georef import warp_with_gcps
import utm





class CCD_TO_GROUND_TERRAIN:

    #CCD_TO_GROUND REQUIRES MISSION CONFIGURATION , H5 INPUT , START_TIME ,SAT,SEN
    def __init__(self, mission_config,h5r_file_location,SAT_ID,SEN_ID,SCENE_START_DATE,SCENE_START_TIME):
        self.mission_config=mission_config
        self.h5r_file_location=h5r_file_location
        self.SAT_ID=SAT_ID
        self.SEN_ID=SEN_ID
        self.SCENE_START_DATE=SCENE_START_DATE
        self.SCENE_START_TIME=SCENE_START_TIME
        #self.DETECTOR_PERCENTILE_FOR_INTERPOLATION=DETECTOR_PERCENTILE_FOR_INTERPOLATION
        #self.LINE_PERCENTILE_FOR_INTERPOLATION=LINE_PERCENTILE_FOR_INTERPOLATION

    # -- GEOMETRY CALCULATION --
    def geom_for_pass(self,DETECTOR_PERCENTILE_FOR_INTERPOLATION,LINE_PERCENTILE_FOR_INTERPOLATION,DEM_FILE_LOCATION=None,GCP_FILE_LOCATION=None,MRC2PLD_ROLL_BIAS=0,MRC2PLD_PITCH_BIAS=0):
        logger.info(" Terrain Correction Geometry Calculation Started  ")
        try:
#--------------------------- LOAD MISSION CONFIG TDI--------------------
            with open(self.mission_config, "r") as mcf:
                mission_config=json.load(mcf)
            #pprint(mission_config)
            h5r_file_location=self.h5r_file_location
            #out_shape_location=grid_vector_location


            GCP_FILE_LOCATION=GCP_FILE_LOCATION
            #GENERATE_GRID=0
            f = h5py.File(h5r_file_location, 'r')
            svs=f['Geo_Location']['Ephemeries_Info']
            #BAND=f['ImageData']['B2']#which band
            SAT_ID=self.SAT_ID
            SEN_ID=self.SEN_ID
            #d=2/0
            NUMBER_OF_BANDS=mission_config[SAT_ID][SEN_ID]["NUMBER_OF_BANDS"]
            NUMBER_OF_DETECTORS=mission_config[SAT_ID][SEN_ID]["NUMBER_OF_DETECTORS"]
            MID_PRINCIPLE_POINT=mission_config[SAT_ID][SEN_ID]["MID_PRINCIPLE_POINT"]

            MRC2PLD_RPY=mission_config[SAT_ID][SEN_ID]["MRC2PLD_RPY"]
            LINE_INTEGRATION_TIME=mission_config[SAT_ID][SEN_ID]["LINE_INTEGRATION_TIME"]

            MRC2PLD_PITCH_BIAS=MRC2PLD_PITCH_BIAS
            MRC2PLD_ROLL_BIAS=MRC2PLD_ROLL_BIAS


            PLD_ROLL=[]
            PLD_PITCH=[]
            PLD_YAW=[]

            SCENE_START_DATE=self.SCENE_START_DATE
            SCENE_START_TIME=self.SCENE_START_TIME


            NUMBER_OF_SCANLINES=f['ImageData']['B2'].shape[0]
            #d=2/0
            DETECTOR_PERCENTILE_FOR_INTERPOLATION=DETECTOR_PERCENTILE_FOR_INTERPOLATION
            LINE_PERCENTILE_FOR_INTERPOLATION=LINE_PERCENTILE_FOR_INTERPOLATION
        except Exception as e:
             logger.error("PARAMETER INITIALIZATION FAILED. MAY ME THE FILE IS CORRUPT OR THE DATA IS NOT IN EXPECTED FORMAT. FIND THE ERROR CODE BELOW")
             logger.exception(e)
             return None



        #pprint(mission_config[SAT_ID][SEN_ID]["PLD_PITCH_SME"])
        #------------------------------------ PAYLOAD_ROLL_PITCH_YAW_INTERPOLATION FOR DETECTORS --------------------------

        logger.info('INTERPOLATING DETECTOR ROLL PITCH YAW VALUES WRT MISSION REFERENCE CUBE ')
        try:
        #ROLL INTERPOLATION

            PLD_ROLL_VALS=mission_config[SAT_ID][SEN_ID]["PLD_ROLL_VALS"]
            pld_roll_val_dets=PLD_ROLL_VALS[0]

            for band in range(1,NUMBER_OF_BANDS+1):
                detector_roll=[]
                slope=(PLD_ROLL_VALS[band][2]-PLD_ROLL_VALS[band][0])/(PLD_ROLL_VALS[0][2]-PLD_ROLL_VALS[0][0])

                for detector in range(1,NUMBER_OF_DETECTORS+1):
                    roll=slope*(detector-PLD_ROLL_VALS[0][0])+PLD_ROLL_VALS[band][0]
                    detector_roll.append(roll)
                PLD_ROLL.append(detector_roll)
            #PITCH_INTERPOLATION

            PLD_PITCH_VALS=mission_config[SAT_ID][SEN_ID]["PLD_PITCH_VALS"]
            #pld_roll_val_dets=PLD_ROLL_VALS[0]

            for band in range(1,NUMBER_OF_BANDS+1):
                detector_pitch=[]
                slope=(PLD_PITCH_VALS[band][2]-PLD_PITCH_VALS[band][0])/(PLD_PITCH_VALS[0][2]-PLD_PITCH_VALS[0][0])

                for detector in range(1,NUMBER_OF_DETECTORS+1):
                    pitch=slope*(detector-PLD_PITCH_VALS[0][0])+PLD_PITCH_VALS[band][0]
                    detector_pitch.append(pitch)
                PLD_PITCH.append(detector_pitch)

            #print(PLD_PITCH[0][2999])

            PLD_YAW_VALS=mission_config[SAT_ID][SEN_ID]["PLD_YAW_VALS"]
            #pld_roll_val_dets=PLD_ROLL_VALS[0]

            for band in range(1,NUMBER_OF_BANDS+1):
                detector_yaw=[]
                slope=(PLD_YAW_VALS[band][2]-PLD_YAW_VALS[band][0])/(PLD_YAW_VALS[0][2]-PLD_YAW_VALS[0][0])

                for detector in range(1,NUMBER_OF_DETECTORS+1):
                    yaw=slope*(detector-PLD_PITCH_VALS[0][0])+PLD_PITCH_VALS[band][0]
                    detector_yaw.append(yaw)
                PLD_YAW.append(detector_yaw)
                #d=2/0
        except Exception as e:
            logger.error("PAYLOAD-MRC INTERPOLATION FAILED - IMPROPER CONFIG OR INFO, PLEASE CHECK BELOW ERROR FOR DETAILS")
            logger.exception(e)
            return None



        #--------------------- PAYLOAD INTERPOLATION FINISH--------------------------------
        pos_list=[]
        time_list=[]
        svc=1
        logger.info('PROCESSING STATE VECTORS FOR INDEXING')
        try:
            for sv in svs:
                date=sv[14].decode('ascii')
                time=sv[15].decode('ascii')
                time_list.append(MISSION_TIME_TO_DATETIME(date,time).replace(tzinfo=timezone.utc).timestamp())
            logger.success('Completed.')
        except Exception as e:
            logger.error("State vector - time indexing failed with below error")
            logger.exception(e)
            return None
        logger.info(' INTERPOLATING STATE VECTORS FOR SCENE CAPTURE TIMES')
        try:
            scene_start_dt=MISSION_TIME_TO_DATETIME_LS(SCENE_START_DATE,SCENE_START_TIME)
            gcps=[]
            #nom_h=300
            scan_line_time=[]
            #printProgressBar(0, l, prefix = 'Progress:', suffix = 'Complete', length = NUMBER_OF_SCANLINES)
            logger.info('CALCULATING TERRAIN GROUND CO-ORDINATES FROM IMAGE COORDINATES')
            for line in range(0,NUMBER_OF_SCANLINES):
                #print('line processing '+str(line))
                if(line%LINE_PERCENTILE_FOR_INTERPOLATION==0):
                    time=scene_start_dt+line*timedelta(seconds=LINE_INTEGRATION_TIME/1000)
                    scan_line_time.append(time)
                    timestamp=time.replace(tzinfo=timezone.utc).timestamp()
                    navigation_lindex=bisect.bisect_left(time_list,timestamp )
                    #print(navigation_lindex)
                    #navigation_rindex=navigation_lindex+1

                    for detector in range(0,NUMBER_OF_DETECTORS):

                        if(detector%DETECTOR_PERCENTILE_FOR_INTERPOLATION==0 or detector==NUMBER_OF_DETECTORS-1):
                            #time_list.append(MISSION_TIME_TO_DATETIME(date,time).replace(tzinfo=timezone.utc).timestamp())
                            #date_str = date+' '+time
                            svbd=svs[navigation_lindex]
                            svad=svs[navigation_lindex+1]

                            svb=[]
                            sva=[]
                            for x in range(0,16):
                                if(x<10):
                                    svb.append(float(svbd[x]))
                                    #print(svb)
                                    sva.append(float(svad[x]))
                                else:
                                    svb.append(svbd[x])
                                    sva.append(svad[x])

                            svb_n=np.array(svb[:10])
                            sva_n=np.array(sva[:10])
                            #print(svb)
                            svbt=MISSION_TIME_TO_DATETIME(svb[14],svb[15]).replace(tzinfo=timezone.utc).timestamp()
                            svat=MISSION_TIME_TO_DATETIME(sva[14],sva[15]).replace(tzinfo=timezone.utc).timestamp()
                            sv=[]
                            for x in range(0,10):
                                #print(float(sva[x]),svb[x])
                                sv.append(((sva_n[x]-svb_n[x])/(svat-svbt))*(svbt-timestamp)+svb_n[x])
                            #sv=((sva-svb)/(svat-svbt))*(svbt-timestamp)+svb
                            qroll=sv[0]
                            qpitch=sv[1]
                            qyaw=sv[2]
                            X=sv[4] #must be Kilo meter right? -RK
                            Xv=sv[5]
                            Y=sv[6]
                            Yv=sv[7]
                            Z=sv[8]
                            Zv=sv[9]
                            # q1=sv[10]
                            # q2=sv[11]
                            # q3=sv[12]
                            # q4=sv[13]
                            #date=sv[14].decode('ascii')
                            #time=sv[15].decode('ascii')
                            t = time

                            state_vector_r=(X,Y,Z)
                            state_vector_v=(Xv,Yv,Zv)

                            #coe=ORBel_frm_SV(state_vector_r,state_vector_v)

                            state_vector_r=(1000*np.array(state_vector_r))
                            state_vector_v=(1000*np.array(state_vector_v))

                            #print("PER_TEST "+str(X*Xv+Y*Yv+Z*Zv))
                            h_cap=(np.cross((state_vector_r),(state_vector_v)))
                            t_cap=NORM3(np.cross(h_cap,state_vector_r))


                            yaw_cap=-NORM3(state_vector_r)
                            roll_cap=t_cap
                            pitch_cap=-NORM3(h_cap)

                            ECI_ORB_DC=np.array([

                            roll_cap,pitch_cap,yaw_cap,

                            ]
                            )

                            sat_rot_theta=[radians(qroll),radians(qpitch),radians(qyaw),]

                            mission_rot_theta=[radians(MRC2PLD_RPY[0]+MRC2PLD_ROLL_BIAS),radians(MRC2PLD_RPY[1]+MRC2PLD_PITCH_BIAS),radians(MRC2PLD_RPY[2]),]
                            trace_list=[]

                            payload_rot_theta=[radians(PLD_ROLL[0][detector]),radians(PLD_PITCH[0][detector]),radians(PLD_YAW[0][detector])]

                            SAT_ROT_MATRIX=eulerAnglesToRotationMatrix(sat_rot_theta)
                            MIS_ROT_MATRIX=eulerAnglesToRotationMatrix(mission_rot_theta)
                            PLD_ROT_MATRIX=eulerAnglesToRotationMatrix(payload_rot_theta)

                            ECI_ORB_INV_DC=np.linalg.inv(ECI_ORB_DC)

                            ROTATED_STATE_VECTOR=np.matmul(ECI_ORB_DC,state_vector_r)



                            ROTATED_STATE_VECTOR=np.matmul(PLD_ROT_MATRIX,ROTATED_STATE_VECTOR)
                            ROTATED_STATE_VECTOR=np.matmul(MIS_ROT_MATRIX,ROTATED_STATE_VECTOR)
                            ROTATED_STATE_VECTOR=np.matmul(SAT_ROT_MATRIX,ROTATED_STATE_VECTOR)


                            ROTATED_STATE_VECTOR=np.matmul(ECI_ORB_INV_DC,ROTATED_STATE_VECTOR)




                            #print('---------- SATELLITE GROUND TRACE VECTOR NORMAL PROJECTION -------------- \n')
                            llh=RAY_INTERSECT_WGS84_TERRAIN(-state_vector_r,ROTATED_STATE_VECTOR,t,DEM_FILE_LOCATION)

                            gcps.append([llh[1],llh[0],llh[2],detector,line])
                            #if(line!=0):
                            #printProgressBar(line,NUMBER_OF_SCANLINES, prefix = 'Progress:', suffix = 'Complete', length = 50)
                            #c=c+1
            logger.success('Terrain co-ordinates calculation completed successfully')
            if(GCP_FILE_LOCATION):

                gcp_json={"GCP_LIST":gcps}
                with open(GCP_FILE_LOCATION, "w") as gcpf:
                    logger.info('Writing Obtained Terrain corrected GCPS to :'+GCP_FILE_LOCATION)
                    #print('WRITING GCPS TO GCP FILE')
                    json.dump(gcp_json,gcpf)

            logger.success('GCPs written to files Succesfully')
            self.gcps=gcps
        except Exception as e:
            logger.error("Ground coordiante calculation failed with below error")
            logger.exception(e)
            return None
        #logger.success('Succesfully calculated GCPS for the pass')
        return 1
        #return gcps

    def gcp_grid(self,gcps,gcp_grid_location):
        try:
            logger.info('MAKING GCP GRIDS FROM GCPS')
            multiline = ogr.Geometry(ogr.wkbMultiLineString)
            for gcp_index in range(1,len(gcps)):
                line1 = ogr.Geometry(ogr.wkbLineString)
                line1.AddPoint(gcps[gcp_index-1][0], gcps[gcp_index-1][1])
                line1.AddPoint(gcps[gcp_index][0], gcps[gcp_index][1])
                multiline.AddGeometry(line1)
            outDriver = ogr.GetDriverByName('ESRI Shapefile')

            # Create the output
            spatialReference = osr.SpatialReference()
            spatialReference.ImportFromEPSG(4326)

            outDataSource = outDriver.CreateDataSource(gcp_grid_location)
            outLayer = outDataSource.CreateLayer('gcp_grid_location',spatialReference, geom_type=ogr.wkbMultiLineString )
            # Get the output Layer's Feature Definition
            featureDefn = outLayer.GetLayerDefn()
            # creae a new feature
            outFeature = ogr.Feature(featureDefn)
            # Set new geometry
            outFeature.SetGeometry(multiline)
            # Add new feature to output Layer
            outLayer.CreateFeature(outFeature)
            # dereference the feature
            outFeature = None
            # Save and close DataSources
            outDataSource = None
            logger.success('GCP GRID PREPARED Succesfully')
            return 1
        except Exception as e:
            logger.error("GEOM NOT CALCULATED. RETURNED NONE WITH BELOW ERROR")
            logger.exception(e)
            return None


    def Warp_Raster(self,gcps,BAND_ID,OUTPUT_LOCATION):


                #line1=[]
        try:
            logger.info('WARPING THE IMAGE WITH THE GCPS')
            f = h5py.File(self.h5r_file_location, 'r')
            BAND=f['ImageData'][BAND_ID]
            #input_path = BAND_LOCATION
            output_path = OUTPUT_LOCATION
            # GCP input
            # xyz = [...]
            # row_col = [...]
            band_data = gdal.GetDriverByName('GTiff').Create(output_path, BAND.shape[1], BAND.shape[0], 1)
            utm_c=transform_wgs84_to_utm(gcps[0][0], gcps[0][1])
            band_srs = osr.SpatialReference()
            band_data.SetProjection(utm_c[1])
            #print([ utm_c[0][0], 23.5, 0,utm_c[0][1], 0, 23.5 ])
            #band_data.SetGeoTransform( [ utm_c[0][0], 23.5, 0,utm_c[0][1], 0, -23.5 ] )
            #ulx,xres,xskew,tly,skewy,res
            band_data_band = band_data.GetRasterBand(1)
            #outData = numpy.zeros((rows,cols), numpy.int16)
            BAND=np.array(BAND)
            band_data_band.WriteArray(BAND, 0, 0)
            band_data.FlushCache()
            gcps_fin = []
            for gcp in gcps:
                gcps_fin.append(gdal.GCP(gcp[0],gcp[1], gcp[2], gcp[3], gcp[4]))
            #band_data=[]
            warp_with_gcps(band_data, output_path, gcps_fin, gcp_epsg=4326, output_epsg=4326)
            logger.success('WARPING IMAGE COMPLETED USING GCPS. OUTPUT AT '+output_path)
            return 1
        except Exception as e:
            logger.error("WARPING NOT DONE. RETURNED NONE WITH BELOW ERROR")
            logger.exception(e)
            return None



if(__name__=="__main__"):

    #-------------CONFIG----------------------------
    # with open("../MISSION_CONFIG_TDI.json", "r") as mcf:
    #     mission_config=json.load(mcf)
    #pprint(mission_config)
    mission_config="../MISSION_CONFIG_TDI.json"
    h5r_file_location='/home/radhakrishna/transfer/dgr/1846810121/BAND.h5'
    SAT_ID="RS2"
    SEN_ID="LISS-III"
    SCENE_START_DATE="2018-08-09"
    SCENE_START_TIME="04:27:26.442982"#21"04:27:06.162866"#11
    GCP_FILE_LOCATION ="/home/radhakrishna/transfer/dgr/1846810121/GCP_PHYSICAL_SENSOR_P125_R043.json"
    GCP_grid_LOCATION ="/home/radhakrishna/transfer/dgr/vectors/21/GCP_PHYSICAL_SENSOR_P125_R043.shp"
    OUTPUT_LOCATION='/home/radhakrishna/transfer/dgr/1846810121/FIN_OUT_BAND2_SENSOR_P125_R043_SVI.tif'
    DEM_FILE_LOCATION='/home/radhakrishna/transfer/dgr/DEM'
    #PROGRAM STARTS
    CCD_TO_GROUND_TR=CCD_TO_GROUND_TERRAIN(mission_config,h5r_file_location,SAT_ID,SEN_ID,SCENE_START_DATE,SCENE_START_TIME)
    rez=CCD_TO_GROUND_TR.geom_for_pass(2000,1000,DEM_FILE_LOCATION=DEM_FILE_LOCATION,GCP_FILE_LOCATION=GCP_FILE_LOCATION,MRC2PLD_ROLL_BIAS=-0.043,MRC2PLD_PITCH_BIAS=-0.125)
    if(rez==1):
        CCD_TO_GROUND_TR.gcp_grid(CCD_TO_GROUND_TR.gcps,GCP_grid_LOCATION)
        CCD_TO_GROUND_TR.Warp_Raster(CCD_TO_GROUND_TR.gcps,"B2",OUTPUT_LOCATION)
