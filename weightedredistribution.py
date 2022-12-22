# -*- coding: utf-8 -*-
"""
Created on Fri May  6 20:35:30 2022

@author: 10853401
"""




# set start time

'''ALL CODE MUST BE INSIDE HERE'''


from numpy import  mean, zeros
from numpy.random import uniform

from math import hypot, pi, sqrt, ceil, floor
from geopandas import read_file

from shapely.geometry import Point
from matplotlib.lines import Line2D

from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from rasterio import open as rio_open
from rasterio.plot import show as rio_show
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.pyplot import subplots, savefig
from matplotlib.patches import Patch 
from time import time
from fiona.errors import DriverError



start_time = time()	# NO CODE ABOVE HERE

################################################################################FUNCTION DEFINITIONS #############################
def distance(a, b):
    
   # This is a function created to compute distance using pythogoras theorem. This will be used in the Weighted Distribution Function below        
       
    return hypot(a[0] - b[0], a[1] - b[1])

def WeightedRedistribution(w,s,pointData,weightingSurface,administrativeAreas):
    
    '''
    *   using a function to implement weighted redistribution where:
    *   w=[user defined] desired influence of weighting surface
    *    s=[user defined] desired level of spatial ambiguity
    *    pointData=the input point dataset to be redistributed
    *    weightingSurface=[user defined] raster data
    *    administrativeAreas=[user defined] polygons for relevant administrative areas at each “level”
    '''
    
    #looping through each administrative area   
    for id, admin in administrativeAreas.iterrows():
                     
       #getting the geometry bounds of each admin area to use in uniform function
       minx, miny, maxx, maxy = admin.geometry.bounds
       
       #extracting all points within the administrative area    
       extractedpoints = pointData.loc[pointData.within(administrativeAreas.loc[id, 'geometry'])]

       # loop through the points in the administrative area
       for id, point in extractedpoints.iterrows():
          # creating a variable random_points to hold w number of random points inside the admin
          maxpoints=0
          #creating a variable maximum to show the maximum value of w, or population density.          
          maximum = 0
          
          #running a loop to get w random points within admin         
          while maxpoints<w:
             #using the uniform random function to generate random points within 
             #admin bounds obtained above 
             
             #obtaining x of random point
             xs = uniform(low=minx, high=maxx, size=1)
             
             #obtaining y of random point
             ys = uniform(low=miny, high=maxy, size=1)
                          
             #initialising it to point rp using shapely Point                 
             rp=Point(xs,ys)
             
             #checking if the random point is within geometry              
             if rp.within(admin.geometry):
                 
                 #getting the x and y of the point in image space               
                 rp1=weightingSurface.index(xs,ys)
                 
                 #assigning the weighting surface value at point rp1 to variable value
                 weightingvalue = weightingSurface.read(1)[rp1]
                 
                 #incrementing random points
                 maxpoints = maxpoints + 1
                 
                 # making sure we get the seed point with the highest population density value 
                 if weightingvalue > maximum:
                     seedpoint = rp1
                     maximum=weightingvalue

          #getting the image space coordinates of the seed point and assignging it to r0 and c0
          r0, c0 =seedpoint
          
          #popping r0 and c0from a list to get their locaton in image space
          x0 = r0.pop()
          y0 = c0.pop()
   
          
          #calculate radius using bounds of geom area and spatial ambiguity 's'
          radius=sqrt(admin.geometry.area*s/pi) 
          
          # convert the radius to pixels by dividing it by the resolution of the dem data
          radius_px=int(radius/weightingSurface.res[0])
          
          #subtracting radius in pixels from seed point values to get circumefrence of circle

          i0=x0-radius_px # row 
          j0=y0-radius_px # col
          
          #calculate v using circle boundaries and make output map
          
          #looping through diameter of circle with radius px
          for i in range(2*radius_px):
              #looping through diameter of circle with radius px

              for j in range(2*radius_px):
              # making sure that we are within the boundary of the dataset
                  if i0+i>0 and i0+i<weightingSurface.height and j0+j>0 and j0+j<weightingSurface.width:
                      
                      #initialising variable T to hold the location between seed point and circle circumference
                      T=(i0+i,j0+j)
                      
                      # check if the cell being drawn is in the circle by measuring the distance between the seedpoint and T is within
                      #the circle radius
                      if distance((x0,y0),T) < radius_px:
                              
                              # calculate v
                              v = 1-sqrt((x0-T[0])*(x0-T[0])+(y0-T[1])*(y0-T[1]))/radius_px

                              
                              x=int(i0+i)
                              y=int(j0+j)
                              # add v to output
                              output[(x,y)] =output[(x,y)]+v
       
    return output

####################################################################IMPORTING DATA AND RUNNING THE FUNCTION ##################################################
# selecting british national grid - ESPG 27700
UKcrs = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs '


#reading Greater Manchester Shapefile

try:
    GreaterMan = read_file("./wr/gm-districts.shp").to_crs(UKcrs)
    
#If the file does not read, exit    
except DriverError:
    print("DriverError:Invalid Filepath. Please check your filepath and try again ")


#reading tweet Shapefile
    
try:
    tweetpoints=read_file("./wr/level3-tweets-subset.shp").to_crs(UKcrs)
    
#If the file does not read, exit       
except DriverError:
    print("DriverError:Invalid Filepath. Please check your filepath and try again ")


#reading weighting layer (population density map)


with rio_open("./wr/100m_pop_2019.tif") as d:
    #reading the band
    band_1=d.read(1)
    
    # create a an empty raster with the same size as the population density map
    output = zeros((d.height, d.width))
    
    # running the weighted distribution algorithm
    WeightedRedistribution(20, 0.1, tweetpoints, d, GreaterMan)
    


    # creating Map of Tweets
    fig, my_ax = subplots(1, 1, figsize=(16, 10))
    my_ax.set(title="Greater Manchester Tweet Distribution")
    
    # plot tweets
    tweetpoints.plot(
        ax = my_ax,
        color = '#000000',
        legend=True,
        edgecolor = '#000000',
        linewidth = 0.5,
        label='Tweets'
        ) 
    # plot Greater Manchester
    GreaterMan.plot(
        ax = my_ax,
        color = '#f0e0e0',
        alpha = 0.4,
        edgecolor = '#660000',
        linewidth = 0.5,
        ) 

 
    #adding legend    
    my_ax.legend(handles=[
            Patch(facecolor='#f0e0e0', edgecolor='#660000', label='Districts'),
            Line2D([0], [0], marker='o', color='w', label='Tweets',
            markerfacecolor='k', markersize=10),
    
        ],loc='lower right',  fontsize='x-large')
           
        
    # add north arrow
    x, y, arrow_length = 0.97, 0.99, 0.1
    my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
    arrowprops=dict(facecolor='white', width=5, headwidth=15),
    ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)
    
    # add scalebar
    my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower left"))  
    
    #save the figure
    savefig('./out/assessment2tweets.png', bbox_inches='tight')    
    
    
    
    # creating weighted redistribution image
    fig, my_ax = subplots(1, 1, figsize=(16, 10))
    my_ax.set(title="Weighted Redistribution of Greater Manchester Tweets")
    
    # plot the new raster layer 'output'
    rio_show(
        output,
        ax=my_ax,
        transform = d.transform,
		cmap = 'jet',

        )
    # plot the Great Manchester vector layer
    GreaterMan.plot(
    ax = my_ax,
    color = 'none',
    edgecolor = '#FFFFFF',
    linewidth = 1,

    ) 
    

    # plotting the scalebar
    colorbar = fig.colorbar(ScalarMappable(norm=Normalize(vmin=floor(band_1.min()), vmax=ceil(band_1.max())), cmap='jet'),ticks=[band_1.min(),band_1.max()], ax=my_ax)
    colorbar.ax.set_yticklabels(['Low', 'High'])  
   
    # add north arrow
    x, y, arrow_length = 0.97, 0.99, 0.1
    my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length), color = 'w',
    arrowprops=dict(facecolor='white', width=5, headwidth=15),
    ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)
    # add scalebar
    my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower left"))
    my_ax.legend(handles=[
        Patch(facecolor='none', edgecolor='white', label='Districts'), ],loc='lower right',  fontsize='x-large')

    # save the figure
    savefig('./out/assessment2.png', bbox_inches='tight')
  
    

	
# report runtime
print(f"completed in: {time() - start_time} seconds")	# NO CODE BELOW HERE    