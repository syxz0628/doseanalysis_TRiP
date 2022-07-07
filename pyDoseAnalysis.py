import os
import numpy as np 
from packaging import version
import matplotlib
import matplotlib.pyplot as plt
import nrrd
from scipy.interpolate import interp1d
# Handle mouse clicks on the plot:
class LineSlice:
    '''Allow user to drag a line on a pcolor/pcolormesh plot, and plot the Z values from that line on a separate axis.

    Example
    -------
    fig, (ax1, ax2) = plt.subplots( nrows=2 )    # one figure, two axes
    img = ax1.pcolormesh( x, y, Z )     # pcolormesh on the 1st axis
    lntr = LineSlice( img, ax2 )        # Connect the handler, plot LineSlice onto 2nd axis

    Arguments
    ---------
    img: the pcolormesh plot to extract data from and that the User's clicks will be recorded for.
    ax2: the axis on which to plot the data values from the dragged line.


    '''
    def __init__(self, img, ax, interactive, **kwargs):
        '''
        img: the pcolormesh instance to get data from/that user should click on
        ax: the axis to plot the line slice on
        '''
        self.img = img
        self.ax = ax
	
        if version.parse(matplotlib.__version__ )>version.parse("3.5.0"): 
            self.data = img.get_array().reshape(img.get_coordinates().shape[0]-1, img.get_coordinates().shape[1]-1)
        else: 
            self.data = img.get_array().reshape(img._meshWidth, img._meshHeight)

        #print(np.where((self.data-data)!=0))

        # register the event handlers:
        if interactive:
            self.cidclick = img.figure.canvas.mpl_connect('button_press_event', self.on_click)
            self.cidrelease = img.figure.canvas.mpl_connect('button_release_event', self.on_click)

        self.markers, self.arrow = None, None   # the lineslice indicators on the pcolormesh plot
        self.line = None    # the lineslice values plotted in a line
        self.out_array = []
        self.hlut = []
        self.mode='profile_only'

    #end __init__

    def on_click(self, event):
        '''Matplotlib will run this function whenever the user triggers an event on our figure'''
        if event.inaxes != self.img.axes: return      # exit if clicks weren't within the `img` axes
        if self.img.figure.canvas.manager.toolbar.mode == '':   # exit if pyplot toolbar (zooming etc.) is active
            if event.name == 'button_press_event':
                self.p1 = (event.xdata, event.ydata)    # save 1st point
            elif event.name == 'button_release_event':
                self.p2 = (event.xdata, event.ydata)    # save 2nd point
                self.drawLineSlice()    # draw the Line Slice position & data
                if self.mode=='wepl':
                    if np.any(self.hlut):
                        self.siddon_sparse_2d()
                    else: 
                        print('Need to load HLUT to perform WEPL estimation!')
                        raise AttributeError
        else: 
            #print(f'event {self.img.figure.canvas.manager.toolbar.mode}')
            return
    #end __call__
    def set_HLUT(self,hlut): 
        hlut_fine = interp1d(hlut[:,0],hlut[:,1])
        x_fine = np.linspace(np.min(hlut[:,0]), np.max(hlut[:,0]),num=int(np.max(hlut[:,0])-np.min(hlut[:,0])), endpoint=True)

        self.hlut = np.c_[x_fine,hlut_fine(x_fine)]

    def set_mode(self,mode):
        self.mode = mode

    def set_pixel_dims(self,dimx,dimy):
        self.dimx = dimx
        self.dimy = dimy

    def drawLineSlice( self ):
        ''' Draw the region along which the Line Slice will be extracted, onto the original self.img pcolormesh plot.  Also update the self.axis plot to show the line slice data.'''
        '''Uses code from these hints:
        http://stackoverflow.com/questions/7878398/how-to-extract-an-arbitrary-line-of-values-from-a-numpy-array
        http://stackoverflow.com/questions/34840366/matplotlib-pcolor-get-array-returns-flattened-array-how-to-get-2d-data-ba
        '''
        x0,y0 = self.p1[0], self.p1[1]  # get user's selected coordinates
        x1,y1 = self.p2[0], self.p2[1]
        print(f'Line profile from (x,y)={self.p1} to (x,y)={self.p2}')

        length = int( np.hypot(x1-x0, y1-y0) )
        x, y = np.linspace(x0, x1, length),   np.linspace(y0, y1, length)

        # Extract the values along the line with nearest-neighbor pixel value:
        # get temp. data from the pcolor plot
        #zi = self.data[x.astype(np.int), y.astype(np.int)]
        # Extract the values along the line, using cubic interpolation:
        import scipy.ndimage
        zi = scipy.ndimage.map_coordinates(self.data.astype(float), np.vstack((y,x)))
        if self.mode == 'wepl':
            if np.any(self.hlut):
                for i,val in enumerate(zi):
                    rsp_index = np.argmin(np.abs(self.hlut[:,0]-val))
                    rsp = self.hlut[rsp_index,1]
                    zi[i] = rsp
            else: 
                print('Need to load HLUT in WEPL mode!')
                raise AttributeError

        # if plots exist, delete them:
        if self.markers != None:
            if isinstance(self.markers, list):
                self.markers[0].remove()
            else:
                self.markers.remove()
        if self.arrow is not None:
            self.arrow.remove()

        # plot the endpoints
        self.markers = self.img.axes.plot([x0, x1], [y0, y1], 'wo')   
        # plot an arrow:
        self.arrow = self.img.axes.annotate("",
                    xy=(x0, y0),    # start point
                    xycoords='data',
                    xytext=(x1, y1),    # end point
                    textcoords='data',
                    arrowprops=dict(
                        arrowstyle="<-",
                        connectionstyle="arc3", 
                        color='white',
                        alpha=0.7,
                        linewidth=3
                        ),

                    )

        # plot the data along the line on provided `ax`:
        if self.line != None:
            self.line[0].remove()   # delete the plot
        self.line = self.ax.plot(zi)
        self.ax.set_ylim([zi.min(), zi.max()])
        self.img.figure.canvas.draw()
        self.img.figure.canvas.flush_events()
        self.out_array = zi
    #end drawLineSlice()

    #Siddon ray tracing 
    """
    Code from public Github :https://github.com/ZhouYzzz/MedicalReconstruction/blob/master/siddon.py
    """
    
    def alpha_range(self,x0, x1, x_min, x_max):
        """
        Helper function to calculate the 1D alpha range used in Siddon algorithm
        :param x0: start point
        :param x1: end point
        :param x_min: border start range
        :param x_max: border end range
        :return: alpha_min, alpha_max
        """
        if x0 == x1:
            raise ValueError('x1 and x2 should be different, get {} and {}'.format(x0, x1))
        alpha_x1 = (x_min - x0) / (x1 - x0)
        alpha_x2 = (x_max - x0) / (x1 - x0)
        alpha_min = max(0, min(alpha_x1, alpha_x2))
        alpha_max = min(1, max(alpha_x1, alpha_x2))
        return alpha_min, alpha_max
    
    
    def alpha_interceptions(self,x0, x1, x_min, x_max, alpha_min, alpha_max, dim, d):
        """
        Helper function to calculate all 1D intercept points
        :param x0: start point
        :param x1: end point
        :param x_min: border start point
        :param x_max: border end point
        :param alpha_min: start alpha
        :param alpha_max: end alpha
        :param dim: number of cells
        :param d: cell size
        :return: a list of intercepted alphas
        """
        full_alpha_interceptions = (x_min + np.arange(dim + 1) * d - x0) / (x1 - x0)
        # print(full_alpha_interceptions)
        i_min = dim - (x_max - (alpha_min if x1 > x0 else alpha_max) * (x1 - x0) - x0) / d
        i_max = 1 + (x0 + (alpha_max if x1 > x0 else alpha_min) * (x1 - x0) - x_min) / d
        # print('i min {}, max {}'.format(i_min, i_max))
        alpha = full_alpha_interceptions[np.arange(i_min, i_max, dtype=np.int64)]
        # print(alpha)
        return alpha
    
    
    def siddon_sparse_2d(self):
        """
        2D siddon algorithm
        :param x0: start x
        :param y0: start y
        :param x1: end x
        :param y1: end y
        :param dims: 2d number of cells
        :param voxel_sizes: cellsize
        :return: sparse matrix elements (i, v)
        """
        voxel_sizes = [self.dimx,self.dimy]
        x0 = self.p1[0]*voxel_sizes[0]
        y0 = self.p1[1]*voxel_sizes[1]
        x1 = self.p2[0]*voxel_sizes[0]
        y1 = self.p2[1]*voxel_sizes[1]

        #print(voxel_sizes)
        dims = self.data.shape
        #print(dims)

        #if x0 == x1 or y0 == y1:
        #    return siddon_special_case_sparse_2d(x0, y0, x1, y1, dims=dims, voxel_sizes=voxel_sizes)
        # calculate the valid alpha range
        a_min_x, a_max_x = self.alpha_range(x0, x1, 0, dims[0] * voxel_sizes[0])
        a_min_y, a_max_y = self.alpha_range(y0, y1, 0, dims[1] * voxel_sizes[1])
        a_min = max(a_min_x, a_min_y)
        a_max = min(a_max_x, a_max_y)
        #print('a min {}, max {}'.format(a_min, a_max))
        if a_min >= a_max:
            return 
        # calculate the total length of the ray
        length = np.sqrt(((x0 - x1))** 2 + ((y0 - y1))** 2)
        #print(length)

        # calculate interceptions (combine x and y interceptions)
        a_inter_x = self.alpha_interceptions(x0, x1, 0, dims[0] * voxel_sizes[0], a_min, a_max, dims[0], voxel_sizes[0])
        a_inter_y = self.alpha_interceptions(y0, y1, 0, dims[1] * voxel_sizes[1], a_min, a_max, dims[1], voxel_sizes[1])
        alpha = np.concatenate(([a_min], a_inter_x, a_inter_y, [a_max]))
        alpha = np.clip(alpha, a_min, a_max)
        alpha = np.unique(alpha)
        # calculate the sparse indices and the values
        WED = 0.
        for v in range(len(alpha) - 1):
            ix = (x0 + np.mean([alpha[v], alpha[v+1]]) * (x1 - x0))/voxel_sizes[0]
            jy = (y0 + np.mean([alpha[v], alpha[v+1]]) * (y1 - y0))/voxel_sizes[1]
            # print('ix {}, jy {}'.format(ix, jy))
            ix = int(ix)
            jy = int(jy)
            assert ix < dims[1] and jy < dims[0]
            lij = length * (alpha[v+1] - alpha[v])
            #print(f'V {lij}')
            rsp_index = np.argmin(np.abs(self.hlut[:,0] - self.data[jy,ix]))
            #print(f' RSP index {rsp_index}')
            rsp = self.hlut[rsp_index,1]
            #print(f'HU {self.data[ix,jy]}; RSP {rsp}')
            wepl = lij*rsp
            WED += wepl
        print(f'The Waterquivalent path length is: {WED}.')
        return 
#end class LineTrace



from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, from_levels_and_colors

class DoseAnalysis: 
    def __init__(self, **kwargs): 
        print("New object")
        self.dose_data = []
        self.wet_data = []
        self.CT = []
        self.target_dose = kwargs.get('dose', 3.)
        self.voi = []
        self.profile_data = []
        self.profile_p1 = []
        self.profile_p2 = []

    def NewColorMap(self):
        # create colormap based on rainbow (less bright than GSI)
        # ---------------
        ## we want the CMAP to be simply "red" between -5% and +7% around target dose
        red_min = int(0.95*256/1.1)
        red_max = int(1.07*256/1.1)
        pink_min = int(red_max)
        pink_max = 256

        rainbow = cm.get_cmap('rainbow', 256)
        bottom = rainbow(np.linspace(0, 1, 256)) #red_min))
        new_colors = rainbow(np.linspace(0,1,256))
        #new_colors = np.vstack((bottom,top))

        pink = np.array([248/256, 24/256, 148/256, 1])
        red = [1,0,0,1] #bottom[-1,:]
        new_colors[red_min:red_max, :] = red
        new_colors[red_max:,:] = pink
        #Convert to cmap
        newcmp = ListedColormap(new_colors)
        return newcmp

    def GenerateGSIColorMap(self): 
        darkBlue = 47.0
        lightBlue = 95.0
        darkGreen = 142.0
        lightGreen = 189.0
        yellow = 225.0
        red = 248.0


        colors = []
        for i in range(0,256):
            if i>=0 and i<darkBlue:
              colors.append((0, 0 + i/darkBlue, 1))
            elif i>=darkBlue and i<lightBlue:
              colors.append((0, 1-0.5*(i-darkBlue)/(lightBlue-darkBlue), 1 - (i-darkBlue)/(lightBlue-darkBlue)))
            elif i>=lightBlue and i<darkGreen:
              colors.append((0, 0.5 + 0.5*(i-lightBlue)/(darkGreen-lightBlue), 0))
            elif i>=darkGreen and i<lightGreen:
              colors.append((0 + (i-darkGreen)/(lightGreen-darkGreen), 1, 0))
            elif i>=lightGreen and i<yellow:
              colors.append((1, 1, 0))
            elif i<=red:
              colors.append((1, 0, 0))
            else:
              colors.append((1, 0, 1))


        n_bins = 256
        cmap_name = 'GSIColorMap_OverdoseOn'
        gsicolors = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

        return gsicolors

    def DiscreteColorMap(self, levels, colors, extend): 
        levels = np.array(levels).astype(float)/100.*self.target_dose

        cmap = from_levels_and_colors(levels, colors, extend=extend)
        return cmap[0]


    def readDose(self, filepath, **kwargs):
        directory = filepath[:filepath.rfind("/")+1]
 
        temppath = "MakeContour_tempFileWoSpaceDir.nrrd"
        with open(filepath) as inFile:
            with open(temppath, "w") as outFile: 
                for line in inFile: 
                    if line.find("space directions")!=-1: 
                        continue
                    if line.find("data file")!=-1: 
                        v = line.split(' ')
                        new = directory+v[-1]
                        line = v[0]+' '+v[1]+' '+new
                    outFile.write(line)
     
        data,header = nrrd.read(temppath)
        os.system("rm "+temppath)
        data = np.array(data)
        data = np.resize(data,(header['sizes'][0], header['sizes'][1], header['sizes'][2]))
        if kwargs.get('overwrite', False) and np.any(self.dose_data): ## do not override if there is nothing to override
            dose_file = kwargs.get('dose_file', -1)
            self.dose_data[dose_file] = data
        else:
            print(f'The dose cube has dimensions {data.shape}')
            self.dose_data.append(data)

    def readHLUT(self, filepath,**kwargs):
        data = np.loadtxt(filepath,comments=['#','!'])
        self.hlut = np.c_[data[:,0], data[:,3]]
        if(kwargs.get('plot',False)):
            plt.plot(self.hlut[:,0],self.hlut[:,1])
            plt.xlabel('Hounsfield units')
            plt.ylabel('Relative stopping power')
            plt.show()

    def readCT(self, filepath): 
        self.CT,header = nrrd.read(filepath)
        self.CT = np.array(self.CT)
        self.CT = np.resize(self.CT, (header['sizes'][0], header['sizes'][1], header['sizes'][2]))
        self.ct_sizes = header['space directions']
        print(f'The CT has dimensions {self.CT.shape}')

    def readWETCube(self, filepath):
        directory = filepath[:filepath.rfind("/")+1]

        temppath = directory+"MakeContour_tempFileWoSpaceDir.nrrd"
        with open(filepath) as inFile:
            with open(temppath, "w") as outFile: 
                for line in inFile: 
                    if line.find("space directions")==-1: outFile.write(line)

    
        data,header = nrrd.read(temppath) 
        os.system("rm "+temppath)
        data = np.array(data)
        data = np.resize(data,(header['sizes'][0], header['sizes'][1], header['sizes'][2]))
        self.wet_data = data

    def readVoi(self, filepath):
        directory = filepath[:filepath.rfind("/")+1]

        temppath = directory+"MakeContour_tempFileWoSpaceDir.nrrd"
        with open(filepath) as inFile:
            with open(temppath, "w") as outFile: 
                for line in inFile: 
                    if line.find("space directions")!=-1: continue
                    if line.find("data file")!=-1: 
                        v = line.split(' ')
                        new = directory+v[-1]
                        line = v[0]+' '+v[1]+' '+new
                    outFile.write(line)
    
        data,header = nrrd.read(temppath) 
        os.system("rm "+temppath)
        data = np.array(data)
        data = np.resize(data,(header['sizes'][0], header['sizes'][1], header['sizes'][2]))
        self.voi.append(data)
    
    def PlotOnCT(self,**kwargs): 
        if not np.any(self.CT):
            print("No CT available")
            raise ValueError
        if not np.any(self.dose_data):
            print("No Dose data available")
            raise ValueError
        dose_vmax = kwargs.get('max_dose', self.target_dose)
        name = kwargs.get('name', "Test")
        showPlot = kwargs.get('show',True)
        color_map = kwargs.get('cmap', 'GSI')
        discrete_colors = kwargs.get('colors',[])
        discrete_levels = kwargs.get('levels',[])
        ct_min = kwargs.get('CT_min', -1000)
        ct_max = kwargs.get('CT_max', 1000)
        dose_alpha = kwargs.get('dose_alpha', 0.7)
        dose_threshold = kwargs.get('dose_threshold',0.1)
        slice_x = kwargs.get('slice_x', slice(0,-1))
        slice_y = kwargs.get('slice_y', slice(0,-1))
        slice_z = kwargs.get('slice_z', int(len(self.CT[0,0,:])/2))
        if not np.isscalar(slice_x) and not np.isscalar(slice_y) and not np.isscalar(slice_z): 
            print('Exactly one of the slices must be scalar for a 2D plot')
            raise ValueError
        elif np.isscalar(slice_x) and (np.isscalar(slice_z) or np.isscalar(slice_y)): 
            print('Exactly one of the slices must be scalar for a 2D plot')
            raise ValueError
        elif np.isscalar(slice_y) and (np.isscalar(slice_z) or np.isscalar(slice_x)): 
            print('Exactly one of the slices must be scalar for a 2D plot') 
            raise ValueError
        elif np.isscalar(slice_z) and (np.isscalar(slice_x) or np.isscalar(slice_y)): 
            print('Exactly one of the slices must be scalar for a 2D plot')
            raise ValueError


        save_format = kwargs.get('save_format','png')
        rotate = kwargs.get('rotate90',1)
        dose_files = kwargs.get('dose_files', [0]) 

        cmap = cm.get_cmap('rainbow', 256) 
        if color_map =='GSI': cmap = self.GenerateGSIColorMap()
        elif color_map == 'discrete':
            #discrete_colors = np.array(discrete_colors)
            #discrete_levels = np.array(discrete_levels)
            if not discrete_colors: 
                print('You need to provide a list of colors to use a discrete color map')
                raise ValueError
            elif not discrete_levels: 
                print('You need to provide a list of levels to use a discrete color map')
                raise ValueError
            else: 
                if len(discrete_colors)==len(discrete_levels)-1: 
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'neither')
                elif len(discrete_colors)==len(discrete_levels):
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'max')
                elif len(discrete_colors)==len(discrete_levels)+1: 
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'both')
                else: 
                    print('The length of the color list must be equal to +-1 the length of the level list.')
                    raise ValueError
        else: 
            try:
               cmap = cm.get_cmap(color_map)
            except: 
                print('The given colormap name could not be found! Please try either "discrete" or "GSI" or a named matplotlib colormap')
                raise ValueError

        im = plt.imshow(np.rot90(self.CT[slice_x,slice_y,slice_z], rotate), origin='lower', cmap='Greys_r', vmin=ct_min, vmax = ct_max)

        if np.any(self.voi):
            
            voi_color = kwargs.get('voi_colorlist',['k']*len(self.voi))
            for voi, color in zip(self.voi,voi_color):
                co = plt.contour(np.rot90(voi[slice_x,slice_y,slice_z], rotate), origin='lower',colors=color,linewidths=[2],levels=[0.1])

        for i in dose_files: 
             dose_masked = np.ma.masked_where(self.dose_data[i]<dose_threshold, self.dose_data[i], copy = True)
             do = plt.imshow(np.rot90(dose_masked[slice_x,slice_y,slice_z]), origin='lower', cmap = cmap,alpha=dose_alpha, vmin=0., vmax = dose_vmax*1.1)
        plt.colorbar(do, label="Dose [Gy]")
        plt.axis('off')
        plt.savefig(name+'.'+save_format, bbox_inches='tight', transparent=False, format=save_format,dpi=300)
        if showPlot:plt.show()
    
    def PlotDoseDifferenceOnCT(self,**kwargs): 
        if not np.any(self.CT):
            print("No CT available.")
            raise ValueError
        if not np.any(self.dose_data):
            print("No Dose data available.")
            raise ValueError
        if not len(self.dose_data)>=2:
            print("Need at least two dose cubes to take difference.")
            raise ValueError

        dose_vmax = kwargs.get('max', 100.) # percent
        dose_vmin = kwargs.get('min', -100.)
        name = kwargs.get('name', "Test")
        showPlot = kwargs.get('show',True)
        color_map = kwargs.get('cmap', 'GSI')
        discrete_colors = kwargs.get('colors',[])
        discrete_levels = kwargs.get('levels',[])
        ct_min = kwargs.get('CT_min', -1000)
        ct_max = kwargs.get('CT_max', 1000)
        dose_alpha = kwargs.get('dose_alpha', 0.7)
        slice_x = kwargs.get('slice_x', slice(0,-1))
        slice_y = kwargs.get('slice_y', slice(0,-1))
        slice_z = kwargs.get('slice_z', int(len(self.CT[0,0,:])/2))
        save_format = kwargs.get('save_format','png')
        dose_indices = kwargs.get('dose_files', [0,1])
        threshold = kwargs.get('threshold', 1) # percent
        rotate = kwargs.get('rotate90', 1)

        cmap = cm.get_cmap('rainbow', 256)

        if color_map =='GSI': cmap = self.GenerateGSIColorMap()
        elif color_map == 'discrete':
            #discrete_colors = np.array(discrete_colors)
            #discrete_levels = np.array(discrete_levels)
            if not discrete_colors:
                print('You need to provide a list of colors to use a discrete color map')
                raise ValueError
            elif not discrete_levels:
                print('You need to provide a list of levels to use a discrete color map')
                raise ValueError
            else:
                if len(discrete_colors)==len(discrete_levels)-1:
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'neither')
                elif len(discrete_colors)==len(discrete_levels):
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'max')
                elif len(discrete_colors)==len(discrete_levels)+1:
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'both')
                else:
                    print('The length of the color list must be equal to +-1 the length of the level list.')
                    raise ValueError
        else:
            try:
               cmap = cm.get_cmap(color_map)
            except:
                print('The given colormap name could not be found! Please try either "discrete" or "GSI" or a named matplotlib colormap')
                raise ValueError

        im = plt.imshow(np.rot90(self.CT[slice_x,slice_y,slice_z],rotate),origin='lower', cmap='Greys_r', vmin=ct_min, vmax = ct_max)

        dose_diff = 100*(self.dose_data[dose_indices[0]] - self.dose_data[dose_indices[1]])/self.dose_data[dose_indices[0]] #relative to input 0

        if np.any(self.voi):    
            voi_color = kwargs.get('voi_colorlist',['k']*len(self.voi))
            for voi, color in zip(self.voi,voi_color):
                co = plt.contour(np.rot90(voi[slice_x,slice_y,slice_z], rotate), origin='lower',colors=color,linewidths=[2],levels=[0.1])


        dose_diff_masked = np.ma.masked_where(np.abs(dose_diff)<threshold, dose_diff, copy = True) ##Mask anything smaller than 0.1
        do = plt.imshow(np.rot90(dose_diff_masked[slice_x,slice_y,slice_z],rotate),origin='lower', cmap = cmap ,alpha=dose_alpha, vmin=dose_vmin, vmax = dose_vmax)
        plt.colorbar(do, label="Dose difference [%]")
        plt.axis('off')
        plt.savefig(name+'.'+save_format, bbox_inches='tight', transparent=False, format=save_format)
        if showPlot: plt.show()

    def GetProfile(self,**kwargs):
        dose_vmax = kwargs.get('max', 1.1*self.target_dose) # percent
        dose_vmin = kwargs.get('min', 0.)
        showPlot = kwargs.get('show',True)
        color_map = kwargs.get('cmap', 'GSI')
        discrete_colors = kwargs.get('colors',[])
        discrete_levels = kwargs.get('levels',[])
        ct_min = kwargs.get('CT_min', -1000)
        ct_max = kwargs.get('CT_max', 1000)
        dose_alpha = kwargs.get('dose_alpha', 0.7)
        slice_x = kwargs.get('slice_x', slice(0,-1))
        slice_y = kwargs.get('slice_y', slice(0,-1))
        slice_z = kwargs.get('slice_z', int(len(self.dose_data[0][0,0,:])/2))
        data_axis = [0,1]
        if np.isscalar(slice_x): 
            data_axis = [1,2]
        elif np.isscalar(slice_y):
            data_axis = [0,2]
        elif np.isscalar(slice_z):
            data_axis = [0,1]
        else: 
            print('At least one slice must be scalar valued for a 2D plot')
            raise ValueError
        rotate = kwargs.get('rotate90',1) 
        dose_indices = kwargs.get('dose_files', [0,1])
        analyze_CT = kwargs.get('analyze_ct', False)
        from_points = kwargs.get('from_points', False)
        p1 = kwargs.get('p1', [])
        p2 = kwargs.get('p2', [])
        interactive = True
        if from_points: 
            if not p1 or not p2: 
                print('Please provide start and end coordinates, if from_points is set to True')
                raise AttributeError
            interactive = False

        cmap = cm.get_cmap('rainbow', 256)


        if color_map =='GSI': cmap = self.GenerateGSIColorMap()
        elif color_map == 'discrete':
            if not discrete_colors:
                print('You need to provide a list of colors to use a discrete color map')
                raise ValueError
            elif not discrete_levels:
                print('You need to provide a list of levels to use a discrete color map')
                raise ValueError
            else:
                if len(discrete_colors)==len(discrete_levels)-1:
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'neither')
                elif len(discrete_colors)==len(discrete_levels):
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'max')
                elif len(discrete_colors)==len(discrete_levels)+1:
                    cmap = self.DiscreteColorMap(discrete_levels,discrete_colors,'both')
                else:
                    print('The length of the color list must be equal to +-1 the length of the level list.')
                    raise ValueError
        else:
            try:
               cmap = cm.get_cmap(color_map)
            except:
                print('The given colormap name could not be found! Please try either "discrete" or "GSI" or a named matplotlib colormap')
                raise ValueError

        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
        img = None

        ## Plot CT (it its not used for analysis)
        ## If we want to analyze it, we plot it as pcolormesh
        ## else, we plot it as imshow
        if not analyze_CT:
            if kwargs.get('plot_ct', False):
                if not np.any(self.CT): 
                    print("A CT scan needs to be available, if option Plot_CT = True")
                    raise ValueError

                ct_plot = np.rot90(self.CT[slice_x,slice_y,slice_z])
                im = ax1.imshow(ct_plot, origin = 'lower', cmap='Greys_r', vmin=ct_min, vmax = ct_max)
        else: 
            if not np.any(self.CT):
                print('A CT scan needs to be available, if option analyze_CT = True')
                raise ValueError

            ct_plot = np.rot90(self.CT[slice_x,slice_y,slice_z])
            img = ax1.pcolormesh( np.arange( len(ct_plot[0,:]) ), np.arange(len(ct_plot[:,0])), ct_plot, shading ='auto', cmap = 'Greys_r', vmin=ct_min, vmax = ct_max)
            self.profile_content = 'CT number'
     
        if kwargs.get('plot_dose', True) or not analyze_CT: 
            if not np.any(self.dose_data):
                print("No Dose data available.")
                raise ValueError
            dose_plot = self.dose_data[dose_indices[0]]
            if kwargs.get('plot_difference', False):
                if not len(self.dose_data)>=2:
                    print('For a difference plot, at least two dose cubes are needed')
                    raise ValueError
                if not len(dose_indices)==2:
                    print('Please provide the numbers (starting from zero) of exactly two of the loaded doses, you want to take the difference between.')
                    raise ValueError
                dose_plot = 100*(self.dose_data[dose_indices[0]] - self.dose_data[dose_indices[1]])/self.dose_data[dose_indices[0]] #relative to input 0
                if not analyze_CT: self.profile_content = 'Difference [%]'
        
            dose_plot = np.rot90(dose_plot[slice_x,slice_y,slice_z], rotate)
            
            # plot the data        
            #img = ax1.imshow(np.rot90(np.flip(dose_plot[slice_x,slice_y,slice_z],1)), cmap = cmap ,alpha=dose_alpha, vmin=dose_vmin, vmax = dose_vmax)
            if not analyze_CT:
                img = ax1.pcolormesh( np.arange( len(dose_plot[0,:]) ), np.arange(len(dose_plot[:,0])), dose_plot, shading ='auto', cmap = cmap,alpha=dose_alpha, vmin=dose_vmin, vmax = dose_vmax) 
                self.profile_content = 'Dose [Gy]'
            else:  
                do = ax1.imshow(dose_plot, origin='lower', cmap = cmap ,alpha=dose_alpha, vmin=dose_vmin, vmax=dose_vmax) 

        fig.colorbar(img, ax=ax1) 

        # register the event handler for interactive plot: 
        LnTr = LineSlice(img, ax2, interactive=interactive)    # args: the pcolor plot (img) & the axis to plot the values on (ax2)
        LnTr.set_mode(kwargs.get('profile_mode', 'profile_only'))
        if analyze_CT and kwargs.get('profile_mode', 'profile_only')=='wepl':
            if not np.any(self.hlut):
                print('A HLUT needs to be available for analyzing the WEPL.')
                raise AttributeError
            self.profile_content = 'Relative stopping power'
            LnTr.set_HLUT(self.hlut)
            if rotate ==1 or rotate ==3: data_axis = data_axis[::-1]
            LnTr.set_pixel_dims(self.ct_sizes[data_axis[0],data_axis[0]],self.ct_sizes[data_axis[1],data_axis[1]])
        if not interactive: 
            LnTr.p1 = p1
            LnTr.p2 = p2
            LnTr.drawLineSlice()

        plt.show()
        if np.any(LnTr.out_array): 
            self.profile_data.append(LnTr.out_array)
            self.profile_p1.append(LnTr.p1)
            self.profile_p2.append(LnTr.p2)
            print(f'The final coordinates of the line profile in the plot were: start (x,y)={LnTr.p1}, end (x,y)={LnTr.p2}')
            print('A new object "profile_data" was created.\nYou can save the data to txt by using the function SaveProfile.\nYou can independently plot it with PlotProfileData.')
            print('You can of course handle the data directly, accessingit with <yourObject>.profile_data')

    def GetProfileFromPoints(self,p1, p2):
        img = ax1.pcolormesh( np.arange( len(dose_plot[0,:]) ), np.arange(len(dose_plot[:,0])), dose_plot, shading ='auto', cmap = cmap,alpha=dose_alpha, vmin=dose_vmin, vmax = dose_vmax)
        fig.colorbar(img, ax=ax1) 
        # register the event handler:
        LnTr = LineSlice(img, ax2, interactive=False)    # args: the pcolor plot (img) & the axis to plot the values on (ax2)
        LnTr.p1 = p1
        LnTr.p2 = p2
        LnTr.drawLineSlice()
        plt.show()
    def PlotProfile(self, **kwargs):
        if not np.any(self.profile_data): 
            print('This funtion is only available if you already created a profile with the GetProfile function.')
            raise ValueError
        profile_indices = kwargs.pop('dose_indices', [0])
        for p in profile_indices: 
            plt.plot(self.profile_data[p], **kwargs)
        plt.xlabel('Pixel')
        plt.ylabel(self.profile_content)
        plt.show()

    def SaveProfileData(self, **kwargs):
        if not np.any(self.profile_data): 
            print('This funtion is only available if you already created a profile with the GetProfile function.')
            raise ValueError

        profile_indices = kwargs.get('profile_indices', [0])
        name = kwargs.get('name', 'Profile_data')
        for p in profile_indices: 
            header = f'These are profile data generated with the GetProfile function in MakeDoseOverlay.py\n The profile starts at point (x,y)={self.profile_p1[p]} and ends at (x,y)={self.profile_p2[p]} in image coordinates.\n'
            np.savetxt(name+'_'+str(p)+'.txt',self.profile_data[p], header=header)

         

import sys
def main():
    dosePlot = DoseAnalysis() ##New object 
    dosePlot.readDose(sys.argv[1])
    dosePlot.readCT(sys.argv[2])
    dosePlot.readVoi(sys.argv[3])
    try:
        dosePlot.PlotOnCT(show=True, slice_z=40,name=sys.argv[4], cmap = 'discrete', levels=[0,10,30,50,70,90,110], colors=["blue", "azure","lightgreen","green", "yellow", "red", "magenta"])
    except ValueError: 
        print("ValueError in plot on CT, data read correctly?")
    
    return 0

if __name__=="__main__":
    main()
