from neuron import h
import neuron
import collections
import numpy as np
import pyqtgraph as pg
import os
import re
import os.path

class HocReader(object):
    """
    Provides useful methods for reading hoc structures.
    
    Parameters
    -----------
        hoc: a hoc object or a "xxx.hoc" file name.
    """
    
    def __init__(self, hoc):
        """
        Parameters
        ----------
        hoc : :obj: `hoc` or str
            Either a hoc object that hs been already created, or a string that defines a hoc file name.
        
        """
        self.file_loaded = False
        if isinstance(hoc, basestring):
            fullfile = os.path.join(os.getcwd(), hoc)
            if not os.path.isfile(fullfile):
                raise Exception("File not found: %s" % (fullfile))
            success = neuron.h.load_file(1, fullfile)
            if success == 0: # indicates failure to read the file
                raise NameError("Found file, but NEURON load failed: %s" % (fullfile))
            self.file_loaded = True
            self.h = h # save a copy of the hoc object itself.
        else:
            self.h = hoc # just use the passed argument
            self.file_loaded = True

        # geometry containers
        self.edges = None
        self.vertexes = None
        
        # all sections in the hoc  {sec_name: hoc Section}
        self.sections = collections.OrderedDict()
        # {sec_name: index} indicates index into self.sections.values() where
        # Section can be found. 
        self.sec_index = {}
        # {sec_name: [mechanism, ...]}
        self.mechanisms = collections.OrderedDict()
        # {sec_group_name: set(sec_name, ...)}
        self.sec_groups = {} 
        
        # topology  {section: (parent, [children,...])}
        self.topology = {}
        
        # populate self.sections, self.sec_index, and self.mechanisms
        self.read_section_info()
        
        # auto-generate section groups based on either hoc section lists, or
        # on section name prefixes.
        sec_lists = self.get_section_lists()
        sec_prefixes = self.get_section_prefixes()

        
        # Add groupings by section list if possible:
        if len(sec_lists) > 1:
            self.add_groups_by_section_list(sec_lists)
            
        # Otherwise, try section prefixes
        elif len(sec_prefixes) > 1:
            for group, sections in sec_prefixes.items():
                self.add_section_group(group, sections)

        # generate topology
        self._generate_topology()

    def get_section(self, sec_name):
        """
        Get the section associated with the section name
        
        Parameters
        ----------
        sec_name : str
            The name of the section object.
        
        Returns
        -------
        The hoc Section object with the given name.
        
        """
        try:
            return self.sections[sec_name]
        except KeyError:
            raise KeyError("No section named '%s'" % sec_name)

    def get_section_prefixes(self):
        """
        Go through all the sections and generate a dictionary mapping their
        name prefixes to the list of sections with that prefix. 
        
        For example, with sections names axon[0], axon[1], ais[0], and soma[0],
        we would generate the following structure:
        
            {'axon': ['axon[0]', 'axon[1]'],
             'ais':  ['ais[0]'],
             'soma': ['soma[0]']}
        """
        prefixes = {}
        regex = re.compile('(?P<prefix>\w+)\[(\d*)\]')
        for sec_name in self.sections:
            g = regex.match(sec_name)
            if g is None:
                continue
            prefix = g.group('prefix')
            prefixes.setdefault(prefix, []).append(sec_name)
        return prefixes

    def get_mechanisms(self, section):
        """
        Get a set of all of the mechanisms inserted into a given section

        Parameters
        ----------
        section : :obj: `NEURON section`
            The NEURON section object.

        Returns
        -------
        A list of mechanism names

        """
        return self.mechanisms[section]

    def get_density(self, section, mechanism):
        """
        Get density mechanism that may be found the section.
        mechanism is a list ['name', 'gbarname']. This is needed because
        some mechanisms do not adhere to any convention and may have a different
        kind of 'gbarname' than 'gbar<name>_mechname'
        returns the average of the conductance density, as that may range across different
        values in a section (e.g., can vary by segments)
        Parameters
        ----------
        section : :obj: `NEURON section`
            The NEURON section object.
        mechanism : list
            mechanism is a list ['name', 'gbarname']. It is used to 
            retrieve the mechanism density from HOC as 
            `segment.name.gbarname`.

        Returns
        -------
            Mean conductance of the selected mechanism in the section, averaged across all segments of the section.
        """
        
        gmech = []
        for seg in section:
            try:
                x =  getattr(seg,  mechanism[0])
                mecbar = '%s_%s' % (mechanism[1], mechanism[0])
                if mecbar in dir(x):
                    gmech.append(getattr(x, mechanism[1]))
                else:
                    print 'hoc_reader:get_density did not find the mechanism in dir x', dir(x)
            except NameError:
                return(0.)
            except:
                print 'hoc_reader:get_density failed to evaluate the mechanisms... '
                raise

#        print gmech
        if len(gmech) == 0:
            gmech = 0.
        return np.mean(gmech)

    def get_sec_info(self, section):
        """
        Get the info of the given section
        modified from: neuronvisio
        
        Parameters
        ----------
        section : :obj: `NEURON section`
            The NEURON section object.
        
        Returns
        -------
        str
            containing the information, with html formatting.
        
        """
        info = "<b>Section Name:</b> %s<br/>" %section.name()
        info += "<b>Length [um]:</b> %f<br/>" % section.L
        info += "<b>Diameter [um]:</b> %f<br/>" % section.diam
        info += "<b>Membrane Capacitance:</b> %f<br/>" % section.cm
        info += "<b>Axial Resistance :</b> %f<br/>" % section.Ra
        info += "<b>Number of Segments:</b> %f<br/>" % section.nseg
        mechs = []
        for seg in section:
            for mech in seg:
                mechs.append(mech.name())
        mechs = set(mechs) # Excluding the repeating ones

        mech_info = "<b>Mechanisms in the section</b><ul>"
        for mech_name in mechs:
            s = "<li> %s </li>" % mech_name
            mech_info += s
        mech_info += "</ul>"
        info += mech_info
        return info

    def read_section_info(self):
        """
        Read all the information about the sections in the current hoc file
        Stores the result in the mechanisms class variable.
        """
        # Collect list of all sections and their mechanism names.
        self.sections = collections.OrderedDict()
        self.mechanisms = collections.OrderedDict()
        for i, sec in enumerate(self.h.allsec()):
            self.sections[sec.name()] = sec
            self.sec_index[sec.name()] = i
            mechs = set()
            for seg in sec:
                for mech in seg:
                    mechs.add(mech.name())
            self.mechanisms[sec.name()] = mechs

    def hoc_namespace(self):
        """
        Get a dict of the HOC namespace {'variable_name': hoc_object}.
        NOTE: this method requires NEURON >= 7.3
        """
        names = {}
        for hvar in dir(self.h): # look through the whole list, no other way
            try:
                # some variables can't be pointed to...
                if hvar in ['nseg', 'diam_changed', 'nrn_shape_changed_', 
                            'secondorder', 'stoprun']: 
                    continue
                u = getattr(self.h, hvar)
                names[hvar] = u
            except:
                continue
        return names
    
    def find_hoc_hname(self, regex):
        """
        Find hoc names matching a pattern
        
        Parameters
        ----------
        regex : str
            Regular expression (Python Re module) to search for.
         
        Returns
        -------
        list 
            The names of HOC objects whose *hname* matches regex.        
        """
        objs = []
        ns = self.hoc_namespace()
        for n, v in ns.items():
            try:
                hname = v.hname()
                if re.match(regex, hname):
                    objs.append(n)
            except:
                continue
        return objs

    def add_section_group(self, name, sections, overwrite=False):
        """
        Declare a grouping of sections (or section names). Sections may be
        grouped by any arbitrary criteria (cell, anatomical type, etc).
        
        Parameters
        ----------
        name : str
                name of the section group
        sections: list
            section names or hoc Section objects.
        
        """
        if name in self.sec_groups and not overwrite:
            raise Exception("Group name %s is already used (use overwrite=True)." % name)
        
        group = set()
        for sec in sections:
            if not isinstance(sec, basestring):
                sec = sec.name()
            group.add(sec)
        self.sec_groups[name] = group

    def get_section_group(self, name):
        """
        Get a section group by name
        Parameters
        ----------
        name : str
            name of the group (dendrite, for example)
        
        Returns
        -------
        The set of section names in the group *name*.
        
        """
        return self.sec_groups[name]
    
    def get_section_lists(self):
        """
        Search through all of the hoc variables to find those that are "SectionLists"
        """
        return self.find_hoc_hname(regex=r'SectionList\[')
        #ns = self.hoc_namespace()
        #return [name for name in ns if ns[name].hname().startswith('SectionList[')]
        
    def add_groups_by_section_list(self, names):
        """
        Add a new section groups from the hoc variables indicated in *names*.
        
        Parameters
        -----------
        names : list 
            List containing variable names as strings. Each name must refer to a list of 
                   Sections in hoc. If a dict is supplied instead, then it
                   maps {hoc_list_name: section_group_name}.
        
        Side effects (modifies)
        -----------------------
           calls add_section_group
        
        """
        # if a list is supplied, then the names of groups to create are 
        # exactly the same as the names of hoc lists.
        if not isinstance(names, dict):
            names = {name:name for name in names}
        for hoc_name, group_name in names.items():
            var = getattr(self.h, hoc_name)
            self.add_section_group(group_name, list(var))

    def get_geometry(self):
        """
        modified from:neuronvisio
        Generate structures that describe the geometry of the sections and their segments (all segments are returned)
        
        Returns
        -------
        vertexes : record array containing {pos: (x,y,z), dia, sec_id}
                           for each segment.
        edges :    array of pairs indicating the indexes of connected
                           vertexes.
        Side effects
        ------------
            Modifies vertexes and edges.
        
        """
        
        # return cached geometry if this method has already run.
        if self.vertexes is not None:
            return self.vertexes, self.edges
        
        self.h.define_shape()

        # map segments (lines) to the section that contains them
        self.segment_to_section = {}
        
        vertexes = []
        connections = []
        
        for secid, sec in enumerate(self.sections.values()):
            x_sec, y_sec, z_sec, d_sec = self.retrieve_coordinate(sec)

            for i,xi in enumerate(x_sec):
                vertexes.append(((x_sec[i], y_sec[i], z_sec[i]), d_sec[i], secid))
                indx_geom_seg = len(vertexes) - 1
                if len(vertexes) > 1 and i > 0:
                    connections.append([indx_geom_seg, indx_geom_seg-1])


        self.edges = np.array(connections)
        self.vertexes = np.empty(len(vertexes), dtype=[
            ('pos', float, 3),
            ('dia', float),
            ('sec_index', int)])
        self.vertexes[:] = vertexes
        return self.vertexes, self.edges

    def retrieve_coordinate(self, section):
        """Retrieve the coordinates of a section avoiding duplicates
        
        Parameters
        ----------
        section : :obj: `NEURON section`
            The NEURON section object.
        
        Returns
        -------
        array
            arrays of x, y, z, d for all the segments of the specified section.
        
        """

        section.push()
        x, y, z, d = [],[],[],[]

        tot_points = 0
        connect_next = False
        for i in range(int(self.h.n3d())):
            present = False
            x_i = self.h.x3d(i)
            y_i = self.h.y3d(i)
            z_i = self.h.z3d(i)
            d_i = self.h.diam3d(i)
            # Avoiding duplicates in the sec
            if x_i in x:
                ind = len(x) - 1 - x[::-1].index(x_i) # Getting the index of last value
                if y_i == y[ind]:
                    if z_i == z[ind]:
                        present = True

            if not present:
                k =(x_i, y_i, z_i)
                x.append(x_i)
                y.append(y_i)
                z.append(z_i)
                d.append(d_i)
        self.h.pop_section()
        return (np.array(x),np.array(y),np.array(z),np.array(d))

    def _generate_topology(self):
        for name, sec in self.sections.items():
            sref = self.h.SectionRef(sec=sec)
            parent = sref.parent().sec.name() if sref.has_parent() else None
            if name not in self.topology:
                self.topology[name] = [None, []]
            self.topology[name][0] = parent
            if parent is not None:
                if parent not in self.topology:
                    self.topology[parent] = [None, []]
                self.topology[parent][1].append(name)
            
    def get_branch(self, root):
        """
        Return all sections in a branch, starting with root.
        
        Parameters
        ----------
        root : :obj: `NEURON section`
            The NEURON section object.
        
        """
        branch = [root]
        childs = [root]
        while len(childs) > 0:
            new_childs = []
            for ch in childs:
                new_childs.extend(self.topology[ch][1])
            childs = new_childs
            branch.extend(childs)
        return branch

    def translate_branch(self, root, dx):
        """
        Move the branch beginning at *root* by *dx*.
        
        Parameters
        ----------
        root : :obj: `NEURON section`
            The NEURON section object.
        dx : array
           Which must be an array of length 3 defining the translation.
        
        """
        self.get_geometry()
        dx[np.newaxis, :]
        for name in self.get_branch(root):
            sid = self.sec_index[name]
            mask = self.vertexes['sec_index'] == sid
            self.vertexes['pos'][mask] += dx

    def make_volume_data(self, resolution=1.0, max_size=500e6):
        """
        Using the current state of vertexes, edges, generates a scalar field
        useful for building isosurface or volumetric renderings.
        
        Parameters
        ----------
        resolution: float, default=1.0 microns
            width (um) of a single voxel in the scalar field.
        max_size: int
            maximum allowed scalar field size (bytes).
        Returns
        -------
            * 3D scalar field indicating distance from nearest membrane,
            * 3D field indicating section IDs of nearest membrane,
            * QTransform that maps from 3D array indexes to original vertex 
                coordinates.
        """
        vertexes, lines = self.get_geometry()
        
        maxdia = vertexes['dia'].max() # maximum diameter (defines shape of kernel)
        kernel_size = int(maxdia/resolution) + 3 # width of kernel
        
        
        # read vertex data
        pos = vertexes['pos']
        d = vertexes['dia']
        sec_id = vertexes['sec_index']
        
        # decide on dimensions of scalar field
        mx = pos.max(axis=0)
        mn = pos.min(axis=0)
        diff = mx - mn
        shape = tuple((diff / resolution + kernel_size).astype(int))

        # prepare blank scalar field for drawing
        size = np.dtype(np.float32).itemsize * shape[0] * shape[1] * shape[2]
        if size > max_size:
            raise Exception("Scalar field would be larger than max_size (%dMB > %dMB), resolution is%f" % (size/1e6, max_size/1e6, resolution))
        scfield = np.zeros(shape, dtype=np.float32)
        scfield[:] = -1000
        
        # array for holding IDs of sections that contribute to each area
        idfield = np.empty(shape, dtype=int)
        idfield[:] = -1
        
        # map vertex locations to voxels
        vox_pos = pos.copy()
        vox_pos -= mn.reshape((1,3))
        vox_pos *= 1./resolution

        # Define kernel used to draw scalar field along dendrites
        def cone(i,j,k):
            # value decreases linearly with distance from center of kernel.
            w = kernel_size / 2
            return w - ((i-w)**2 + (j-w)**2 + (k-w)**2)**0.5
        kernel = resolution * np.fromfunction(cone, (kernel_size,)*3)
        kernel -= kernel.max()

        def array_intersection(arr1, arr2, pos):
            """
            Return slices used to access the overlapping area between two 
            arrays that are offset such that the origin of *arr2* is a *pos* 
            relative to *arr1*.            
            """
            s1 = [0]*3
            s2 = [0]*3
            t1 = [0]*3
            t2 = [0]*3
            pos = map(int, pos)
            for axis in range(3):
                s1[axis] = max(0, -pos[axis])
                s2[axis] = min(arr2.shape[axis], arr1.shape[axis]-pos[axis])
                t1[axis] = max(0, pos[axis])
                t2[axis] = min(arr1.shape[axis], pos[axis]+arr2.shape[axis])
            slice1 = (slice(t1[0],t2[0]), slice(t1[1],t2[1]), slice(t1[2],t2[2]))
            slice2 = (slice(s1[0],s2[0]), slice(s1[1],s2[1]), slice(s1[2],s2[2]))
            return slice1, slice2            

        # Draw lines into volume using *kernel* as the brush
        vox_pos[:,0] = np.clip(vox_pos[:,0], 0, scfield.shape[0]-1)
        vox_pos[:,1] = np.clip(vox_pos[:,1], 0, scfield.shape[1]-1)
        vox_pos[:,2] = np.clip(vox_pos[:,2], 0, scfield.shape[2]-1)
        for c in range(lines.shape[0]):
            i = lines[c, 0]
            j = lines[c, 1]
            p1 = vox_pos[i].copy()
            p2 = vox_pos[j].copy()
            diff = p2-p1
            axis = np.argmax(np.abs(diff))
            dia = d[i]
            nvoxels = abs(int(diff[axis]))+1
            for k in range(nvoxels):
                kern = kernel + (dia/2.0)
                sl1, sl2 = array_intersection(scfield, kern, p1) # find the overlapping area between the field and the kernel
                idfield[sl1] = np.where(scfield[sl1] > kern[sl2], idfield[sl1], sec_id[i])
                scfield[sl1] = np.where(scfield[sl1] > kern[sl2], scfield[sl1], kern[sl2])
                dia += (d[j]-d[i]) / nvoxels
                p1 += diff / nvoxels
                
        # return transform relating volume data to original vertex data
        transform = pg.Transform3D()
        w = resolution * kernel_size / 2 # offset introduced due to kernel
        transform.translate(*(mn-w))
        transform.scale(resolution, resolution, resolution)
        transform.translate(1, 1, 1)
        return scfield, idfield, transform
