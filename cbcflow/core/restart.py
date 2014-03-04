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

from cbcflow.dol import compile_extension_module, Function, HDF5File

import shelve
import os
import subprocess
from commands import getstatusoutput

from cbcflow.utils.common import cbcflow_warning

fetchable_formats = ["hdf5", "xml"]

def find_common_savetimesteps(play_log, fields):
    common_keys = []
    for ts, data in play_log.items():
        if "fields" not in data:
            continue
        
        present = {}
        for f in fields:
            present[f] = False
            for v in data["fields"].values():
                if f == v['type'] and any([saveformat in fetchable_formats for saveformat in v["save_as"]]):
                    present[f] = True

        if all(present.values()):
            common_keys.append(int(ts))
    
    return list(sorted(common_keys))


class Restart(object):
    """
    Class to specify restart data. Searches within the given postprocessors
    casedir for solution data to find a suitable timestep to restart from,
    based on the parameters restart_time or restart_timestep.
    
    Overwrites problem.initial_conditions with solution loaded from file and sets
    problem.params.T0. Handles postprocessing to avoid overwrite conflicts with
    existing data.
    """
    def __init__(self, problem, postprocessor, restart_time, restart_timestep):
        self.problem = problem
        self.postprocessor = postprocessor
                
        assert not (restart_time > 0 and restart_timestep > 0), "Restart failed. Cannot specify both restart timestep and restart time."
        assert self.postprocessor != None, "Restart failed. Must specify postprocessor with associated casedir."
        
        # If time and timestep not given, restart from last existing timestep
        if restart_timestep < 0:
            restart_timestep = 1e8
               
        self.casedir = self.postprocessor.get_casedir()
        
        play_log = shelve.open(os.path.join(self.casedir, "play.db"), 'c')
        timesteps = find_common_savetimesteps(play_log, ["Velocity", "Pressure"])
        
        assert len(timesteps) > 0, "Unable to find any timesteps to start from."
        
        times = [play_log[str(ts)]["t"] for ts in timesteps]

        if restart_time > 0:
            timediffs = [abs(t-restart_time) for t in times]
            i = timediffs.index(min(timediffs))
        else:
            timestepdiffs = [abs(ts-restart_timestep) for ts in timesteps]
            i = timestepdiffs.index(min(timestepdiffs))
        
        # Replace problems initial_condition function    
        self._replace_ic(problem, timesteps[i], play_log[str(timesteps[i])])
        
        
        # Replace update problem parameters
        problem.params.start_timestep = timesteps[i]
        problem.params.T0 = times[i]
        
        # Correct postprocessing fields
        self._correct_postprocessing(play_log, timesteps[i])
    
    def _replace_ic(self, problem, timestep, play_log_item):
        "Replace problem.initial_conditions with data acquired from the play log."
        
        # Find velocity metadata
        keys = {}
        for k,v in play_log_item['fields'].items():
            keys[v['type']] = k
        
        u_metadata = shelve.open(os.path.join(self.casedir, keys["Velocity"], "metadata.db"), 'r')[str(timestep)]
        p_metadata = shelve.open(os.path.join(self.casedir, keys["Pressure"], "metadata.db"), 'r')[str(timestep)]
        
        def initial_conditions(spaces, controls):
            """ Function to replace existing problem.initial_conditions """
            # Velocity initial conditions
            V = spaces.V
            if 'hdf5' in u_metadata.keys():
                hdf5filename = os.path.join(self.casedir, keys["Velocity"], keys["Velocity"]+".hdf5")
                hdf5file = HDF5File(hdf5filename, 'r')
                u = Function(V)
                hdf5file.read(u, u_metadata['hdf5']['dataset'])
                del hdf5file
            elif "xml" in u_metadata.keys():
                xmlfilename = os.path.join(self.casedir, keys["Velocity"], keys["Velocity"]+str(timestep)+".xml")
                u = Function(V, xmlfilename)
            elif "xml.gz" in u_metadata.keys():
                xmlfilename = os.path.join(self.casedir, keys["Velocity"], keys["Velocity"]+str(timestep)+".xml.gz")
                u = Function(V, xmlfilename)

            # TODO: Allow for two initial conditions
            icu = [ui for ui in u]
            
            # Pressure initial conditions
            Q = spaces.Q
            if 'hdf5' in p_metadata.keys():
                hdf5filename = os.path.join(self.casedir, keys["Pressure"], keys["Pressure"]+".hdf5")
                hdf5file = HDF5File(hdf5filename, 'r')
                p = Function(Q)
                hdf5file.read(p, p_metadata['hdf5']['dataset'])
            elif "xml" in p_metadata.keys():
                xmlfilename = os.path.join(self.casedir, keys["Pressure"], keys["Pressure"]+str(timestep)+".xml")
                u = Function(Q, xmlfilename)
            elif "xml.gz" in p_metadata.keys():
                xmlfilename = os.path.join(self.casedir, keys["Pressure"], keys["Pressure"]+str(timestep)+".xml.gz")
                u = Function(Q, xmlfilename)
            
            icp = p
            
            return icu, icp
        
        # Replace problem.initial_conditions, and update parameters
        problem.initial_conditions = initial_conditions

    def _correct_postprocessing(self, play_log, restart_timestep):
        "Removes data from casedir found at timestep>restart_timestep."
        play_log_to_remove = {}
        
        for k,v in play_log.items():
            if int(k) > restart_timestep:
                play_log_to_remove[k] = play_log.pop(k)
        
        all_fields_to_clean = []
        
                
        for k,v in play_log_to_remove.items():
            if not "fields" in v:
                continue
            else:
                all_fields_to_clean += v["fields"].keys()
        all_fields_to_clean = list(set(all_fields_to_clean))
        for fieldname in all_fields_to_clean:
            self._clean_field(fieldname, restart_timestep)
    
    def _clean_field(self, fieldname, restart_timestep):
        "Deletes data from field found at timestep>restart_timestep."
        metadata = shelve.open(os.path.join(self.postprocessor.get_savedir(fieldname), 'metadata.db'), 'w')

        metadata_to_remove = {}
        for k in metadata.keys():
            try:
                k = int(k)
            except:
                continue
            if k > restart_timestep:
                metadata_to_remove[str(k)] = metadata.pop(str(k))
        
        # Remove files and data for all save formats
        self._clean_hdf5(fieldname, metadata_to_remove)
        self._clean_files(fieldname, metadata_to_remove)
        self._clean_txt(fieldname, metadata_to_remove)
        self._clean_shelve(fieldname, metadata_to_remove)
        self._clean_xdmf(fieldname, metadata_to_remove)
        self._clean_pvd(fieldname, metadata_to_remove)
    
    def _clean_hdf5(self, fieldname, del_metadata):        
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
        
        hdf5filename = os.path.join(self.postprocessor.get_savedir(fieldname), fieldname+'.hdf5')
        if not os.path.isfile(hdf5filename):
            return
        for k, v in del_metadata.items():
            if not 'hdf5' in v:
                continue
            else:
                cpp_module.delete_from_hdf5_file(hdf5filename, v['hdf5']['dataset'])
        hdf5tmpfilename = os.path.join(self.postprocessor.get_savedir(fieldname), fieldname+'_tmp.hdf5')
        
        status, result = getstatusoutput("h5repack")
        if status != 0:
            cbcflow_warning("Unable to run h5repack. Will not repack hdf5-files before replay, which may cause bloated hdf5-files.")
        else:
            subprocess.call("h5repack %s %s" %(hdf5filename, hdf5tmpfilename), shell=True)
            os.remove(hdf5filename)
            os.rename(hdf5tmpfilename, hdf5filename)
        
        
    def _clean_files(self, fieldname, del_metadata):
        for k, v in del_metadata.items():
            if not 'filename' in v:
                continue
            else:
                fullpath = os.path.join(self.postprocesor.get_savedir(fieldname), v['filename'])
                os.remove(fullpath)

    def _clean_txt(self, fieldname, del_metadata):
        txtfilename = os.path.join(self.postprocessor.get_savedir(fieldname), fieldname+".txt")
        if not os.path.isfile(txtfilename):
            return
        
        txtfile = open(txtfilename, 'r')
        txtfilelines = txtfile.readlines()
        txtfile.close()
        
        num_lines_to_strp = ['txt' in v for v in del_metadata.values()].count(True)
        
        txtfile = open(txtfilename, 'w')
        [txtfile.write(l) for l in txtfilelines[:-num_lines_to_strp]]
        txtfile.close()
        
    def _clean_shelve(self, fieldname, del_metadata):
        shelvefilename = os.path.join(self.postprocessor.get_savedir(fieldname), fieldname+".db")
        if not os.path.isfile(shelvefilename):
            return
        
        shelvefile = shelve.open(shelvefilename, 'c')
        for k,v in del_metadata.items():
            if 'shelve' in v:
                shelvefile.pop(str(k))
        shelvefile.close()
    
    def _clean_xdmf(self, fieldname, del_metadata):
        basename = os.path.join(self.postprocessor.get_savedir(fieldname), fieldname)
        if not os.path.isfile(basename+".xdmf"):
            return
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
    