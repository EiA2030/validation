import ee
import datetime
import time
import sys
from google.cloud import storage

out = str(sys.argv[5])

ee.Initialize()

xmin = float(sys.argv[1])
xmax = float(sys.argv[3])
ymin = float(sys.argv[2])
ymax = float(sys.argv[4])

bb = ee.Geometry.BBox(xmin, ymin, xmax, ymax)

# Clip images function
def clipbb(image):
  return image.clip(bb)

sdate = datetime.date.today()-datetime.timedelta(days = 3650)
edate = datetime.date.today()
sdate = ee.Date.fromYMD(int(sdate.strftime("%Y")),
                        int(sdate.strftime("%-m")),
                        int(sdate.strftime("%-d")))
edate = ee.Date.fromYMD(int((edate).strftime("%Y")),
                        int((edate).strftime("%-m")),
                        int((edate).strftime("%-d")))

def means(col, band, sdate, edate, bb):
  img = ee.ImageCollection(col).filterDate(sdate, edate).filterBounds(bb).select(band).map(clipbb).mean().rename(band)
  return img

def exprt(img, filename, bb):
  task = ee.batch.Export.image.toCloudStorage(**{
    'image': img,
    'scale': 1000,
    'region': bb,
    'crs': 'EPSG:4326',
    'maxPixels' : 1e10,
    'fileNamePrefix' : str(filename),
    'fileFormat': 'GeoTIFF'
  })
  task.start()
  while task.active():
    print('Polling for task (id: {}).'.format(task.id))
    time.sleep(5)
  blob = bucket.blob(str(filename) + ".tif")
  blob.download_to_filename(out + "/" + str(filename) + ".tif")
  
bio01 = ee.Image('WORLDCLIM/V1/BIO').select('bio01').expression('b("bio01")*0.1')  
bio12 = ee.Image('WORLDCLIM/V1/BIO').select('bio12')
srtm = ee.Image('CGIAR/SRTM90_V4').select('elevation').reproject(bio01.projection())
soc = ee.Image('projects/soilgrids-isric/soc_mean').select('soc_5-15cm_mean').expression('b("soc_5-15cm_mean")*10').reproject(bio01.projection())
ph = ee.Image('projects/soilgrids-isric/phh2o_mean').select('phh2o_5-15cm_mean').expression('b("phh2o_5-15cm_mean")*10').reproject(bio01.projection())
cec = ee.Image('projects/soilgrids-isric/cec_mean').select('cec_5-15cm_mean').expression('b("cec_5-15cm_mean")*10').reproject(bio01.projection())
clay = ee.Image('projects/soilgrids-isric/clay_mean').select('clay_5-15cm_mean').expression('b("clay_5-15cm_mean")*10').reproject(bio01.projection())
sand = ee.Image('projects/soilgrids-isric/sand_mean').select('sand_5-15cm_mean').expression('b("sand_5-15cm_mean")*10').reproject(bio01.projection())
ndvi = means(col = "MODIS/061/MOD13A2", band = "NDVI", sdate = sdate, edate = edate, bb = bb).expression('b("NDVI")*0.0001').reproject(bio01.projection())

exprt(img = bio01, filename = 'bio01',  bb = bb)
exprt(img = bio12, filename = 'bio12',  bb = bb)
exprt(img = srtm, filename = 'srtm',  bb = bb)
exprt(img = soc, filename = 'soc',  bb = bb)
exprt(img = ph, filename = 'ph',  bb = bb)
exprt(img = cec, filename = 'cec',  bb = bb)
exprt(img = clay, filename = 'clay',  bb = bb)
exprt(img = sand, filename = 'sand',  bb = bb)
exprt(img = ndvi, filename = 'ndvi',  bb = bb)
