# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.
__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2012-02-17"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ..dol import *

from numpy import linspace

import os
from glob import glob
import re
import subprocess

from ..core.utils import cbcflow_warning, cbcflow_print


def clean_hdf5_file(filename, datasets):
    delete_from_hdf5_file = '''
    namespace dolfin {
        #include <hdf5.h>  
        void delete_from_hdf5_file(std::string filename, std::string dataset)
        {
            const hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
            // Open file existing file for append
            hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
            
            H5Ldelete(file_id, dataset.c_str(), H5P_DEFAULT);
            
            herr_t status = H5Fclose(file_id);
        }
    }
    '''

    cpp_module = compile_extension_module(delete_from_hdf5_file)
    
    hdf5_file = "%s/%s/%s.hdf5" %(self.casedir, fieldname, fieldname)
    if not os.path.isfile(hdf5_file):
        cbcflow_warning("Unable to find file %s. Likely discrepancy between metadata.txt and data." %hdf5_file)
        return
    
    if len(datasets) > 0:
        for dataset in datasets:
            cpp_module.delete_from_hdf5_file(filename, dataset)
        
        tmp_hdf5_file = "%s/%s/%s_tmp.hdf5" %(self.casedir, fieldname, fieldname)
            
        subprocess.call("h5repack %s %s" %(hdf5_file, tmp_hdf5_file), shell=True)
        os.remove(hdf5_file)
        os.rename(tmp_hdf5_file, hdf5_file)
    

class Restart(object):
    def __init__(self, problem, postprocessor, restart_time, restart_timestep):
        self.problem = problem
        self.postprocessor = postprocessor
                
        assert not (restart_time > 0 and restart_timestep > 0), "Restart failed. Cannot specify both restart timestep and restart time."
        assert self.postprocessor != None, "Restart failed. Must specify postprocessor with associated casedir."
        
        # If time and timestep not given, restart from last existing timestep
        if restart_timestep < 0:
            restart_timestep = 1e8
               
        self.casedir = self.postprocessor.params.casedir
        
        u_metadata = self._read_metadata("Velocity")
        assert u_metadata != None, "Unable to find Velocity metadata."
        p_metadata = self._read_metadata("Pressure")
        
        
        u_time, u_timestep, u_datatype, u_name = self._find_restart_data(u_metadata, restart_time, restart_timestep)
        
        
        if u_datatype == "dataset":
            u_filepath = self.casedir + "/Velocity/Velocity.hdf5"
        elif u_datatype == "filename":
            u_filepath = self.casedir + "/Velocity/" + u_name
        
        try:
            p_time, p_timestep, p_datatype, p_name = self._find_restart_data(p_metadata, u_time, u_timestep)
        except:
            p_time, p_timestep, p_datatype, p_name, p_filepath = None, None, None, None, None
        
        if p_datatype == "dataset":
            p_filepath = self.casedir + "/Pressure/Pressure.hdf5"
        elif p_datatype == "filename":
            p_filepath = self.casedir + "/Pressure/" + p_name
            
        if u_timestep != p_timestep or not os.path.isfile(p_filepath):
            cbcflow_warning("Unable to find matching restart info for pressure. Using Constant(0).")
        
        # Replace problems initial_condition function    
        self._replace_ic(problem, u_datatype, u_name, u_filepath, p_datatype, p_name, p_filepath)

        # Replace update problem parameters
        problem.params.start_timestep = u_timestep
        problem.params.T0 = u_time
        
        # Correct postprocessing fields
        self._correct_postprocessing(postprocessor, u_time, u_timestep)
    
    def _replace_ic(self, problem, u_datatype, u_name, u_filepath, p_datatype, p_name, p_filepath):
        def initial_conditions(spaces, controls):
            """ Function to replace existing problem.initial_conditions """
            # Velocity initial conditions
            V = spaces.V            
            if u_datatype == "filename":
                u = Function(V, u_filepath)
            elif u_datatype == "dataset":
                u = Function(V)
                f = HDF5File(u_filepath, 'r')
                f.read(u, u_name)
            
            # TODO: Expand icu to contain two timesteps (required for some solvers)
            icu = [u[0], u[1], u[2]]
            
            # Pressure initial conditions
            Q = spaces.Q
            if p_datatype == "filename":
                p = Function(Q, p_filepath)
            elif p_datatype == "dataset":
                p = Function(Q)
                f = HDF5File(p_filepath, 'r')
                f.read(p, p_name)
            else:
                p = Constant(0)
            
            icp = p
            
            return icu, icp
        
        # Replace problem.initial_conditions, and update parameters
        problem.initial_conditions = initial_conditions
        
    def _read_metadata(self, fieldname):
        filepath = self.casedir+"/"+fieldname+"/metadata.txt"
        if not os.path.isfile(filepath):
            return None
        
        metadata_f = open(filepath, 'r')
        metadata_lines = metadata_f.readlines()
        metadata_f.close()
        return metadata_lines
    
    def _find_restart_data(self, metadata, restart_time, restart_timestep):
        metadata_split = [line.split('\t') for line in metadata if "save_count" in line]
        if len(metadata_split) == 0:
            return None, None, None, None, None
        
        save_counts = []
        timesteps = []
        times = []
        saved_data = []
        for line in metadata_split:
            save_counts.append(int(line[0].split('=')[1]))
            timesteps.append(int(line[1].split('=')[1]))
            times.append(float(line[2].split('=')[1]))
        
            saved_data.append(line[3:])
            
            
        if restart_time > 0:
            timediffs = [abs(t-restart_time) for t in times]
            i = timediffs.index(min(timediffs))
        else:
            timestepdiffs = [abs(ts-restart_timestep) for ts in timesteps]
            i = timestepdiffs.index(min(timestepdiffs))
            
        restart_time = times[i]
        restart_timestep = timesteps[i]
        
        j = ["dataset" in sd or "filename" in sd for sd in saved_data[i]].index(True)
        datatype, name = saved_data[i][j].replace("'", '').replace('"', '').replace('\n', '').split('=')
        
        return restart_time, restart_timestep, datatype, name
        
    def _correct_postprocessing(self, postprocessor, restart_time, restart_timestep):
        for fieldname in postprocessor._fields:
            metadata = self._read_metadata(fieldname)
            del_metadata = []
            new_metadata = []
            if metadata == None:
                break
            
            for line in metadata:
                if "save_count" in line:
                    spl_line = line.split('\t')
                    t = float(spl_line[2].split('=')[1])
                    if t > restart_time:
                        del_metadata.append(line)
                        continue
                    else:
                        save_count = int(spl_line[0].split('=')[1])
                
                new_metadata.append(line)
            
            postprocessor._save_counts[fieldname] = save_count + 1
            
            self._rewrite_metadata(fieldname, new_metadata)
            self._clean_data(fieldname, del_metadata)
            self._rewrite_metadata(fieldname, new_metadata)

    def _clean_metadata(self, metadata, restart_timestep):
        
        
        del_metadata = []
        return del_metadata
    
    def _clean_data(self, fieldname, del_metadata):
        try:
            self._clean_hdf5(fieldname, del_metadata)
        except:
            cbcflow_warning("Cleaning hdf5-file for restart failed. hdf5-file might be corrupted.")
        
        try:
            self._clean_xml(fieldname, del_metadata)
        except:
            cbcflow_warning("Cleaning xml-files for restart failed. xml-files might be corrupted.")
        
        try:
            self._clean_xmlgz(fieldname, del_metadata)
        except:
            cbcflow_warning("Cleaning xml.gz-files for restart failed. xml-files might be corrupted.")
        
        try:
            self._clean_xdmf(fieldname, del_metadata)
        except:
            cbcflow_warning("Cleaning xdmf-file for restart failed. xdmf-file might be corrupted.")
        
        self._clean_pvd(fieldname, del_metadata)
    
    def _clean_hdf5(self, fieldname, del_metadata):
        del_datasets = []
        for line in del_metadata:
            if "dataset=" in line:
                del_datasets.append(re.search("dataset='(.\S+)'", line).group(1))
                       
        if len(del_datasets) == 0:
            return
        
        delete_from_hdf5_file = '''
        namespace dolfin {
            #include <hdf5.h>  
            void delete_from_hdf5_file(std::string filename, std::string dataset)
            {
                const hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
                // Open file existing file for append
                hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
                
                H5Ldelete(file_id, dataset.c_str(), H5P_DEFAULT);
                
                herr_t status = H5Fclose(file_id);
            }
        }
        '''
    
        cpp_module = compile_extension_module(delete_from_hdf5_file)
        
        hdf5_file = "%s/%s/%s.hdf5" %(self.casedir, fieldname, fieldname)
        if not os.path.isfile(hdf5_file):
            cbcflow_warning("Unable to find file %s. Likely discrepancy between metadata.txt and data." %hdf5_file)
            return
        
        for dataset in del_datasets:
            cpp_module.delete_from_hdf5_file(hdf5_file, dataset)

        tmp_hdf5_file = "%s/%s/%s_tmp.hdf5" %(self.casedir, fieldname, fieldname)
            
        subprocess.call("h5repack %s %s" %(hdf5_file, tmp_hdf5_file), shell=True)
        os.remove(hdf5_file)
        os.rename(tmp_hdf5_file, hdf5_file)
        
    def _clean_xml(self, fieldname, del_metadata):
        del_files = []
        for line in del_metadata:
            p = re.compile("filename='(.\S+.xml)'")
            if re.search(p, line):
                del_files.append(re.search(p, line).group(1))

        for fil in del_files:
            path = "%s/%s/%s" %(self.casedir, fieldname, fil)
            if os.path.isfile(path):
                os.remove(path)
    
    def _clean_xmlgz(self, fieldname, del_metadata):
        del_files = []
        for line in del_metadata:
            p = re.compile("filename='(.\S+.xml.gz)'")
            if re.search(p, line):
                del_files.append(re.search(p, line).group(1))

        for fil in del_files:
            path = "%s/%s/%s" %(self.casedir, fieldname, fil)
            if os.path.isfile(path):
                os.remove(path)
    
    def _clean_xdmf(self, fieldname, del_metadata):
        basename = "%s/%s/%s" %(self.casedir, fieldname, fieldname)
        i = 0
        while True:
            h5_filename = basename+"_RS"+str(i)+".h5"
            if not os.path.isfile(h5_filename):
                break
            i = i + 1
        
        xdmf_filename = basename+"_RS"+str(i)+".xdmf"
        
        os.rename(basename+".h5", h5_filename)
        os.rename(basename+".xdmf", xdmf_filename)
        
        f = open(xdmf_filename, 'r').read()
        
        new_f = open(xdmf_filename, 'w')
        new_f.write(f.replace(os.path.split(basename)[1]+".h5", os.path.split(h5_filename)[1]))
        new_f.close()
    
    def _clean_pvd(self, fieldname, del_metadata):
        if os.path.isfile(self.casedir+"/"+fieldname+"/"+fieldname+".pvd"):
            cbcflow_warning("No functionality for cleaning pvd-files for restart. Will overwrite.")
    
    def _rewrite_metadata(self, fieldname, metadata):
        f = open(self.casedir+"/"+fieldname+"/metadata.txt", 'w')
        for line in metadata:
            f.write(line)
        f.write("### RESTART ###\n")
        f.close()

    
