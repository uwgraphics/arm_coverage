import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import time
import subprocess
import pyvista as pv
from utils import *
import time
from watchdog.observers import Observer
import static_variables
from solver_base import SolverBase

class JointGTSP(SolverBase):
    def __init__(self, file_name, tolerances):
        SolverBase.__init__(self, file_name, tolerances)
        self.motions = Motions()
        self.routines = []

        f = open("./GLKH/joint_gtsp.tour", 'a+') 
        f.close()
        self.file_handler = self.TourFileChangeHandler(self, "./GLKH/joint_gtsp.tour", self.file_callback)

    def file_callback(self, routine):
        print("Got a motion")
        if routine not in self.routines:
            self.routines.append(routine)
            j_motion = [self.ik_solutions[i] for i in routine]
            c_path = [self.poses[self.ik_i_to_pose_i[i]] for i in routine]
            self.motions.c_paths.append(c_path)
            self.motions.motions.append(j_motion)
            self.motions.times.append(time.time())
            self.motions.n += 1

    def run_glkh(self):
        p = subprocess.Popen(["./GLKH", "./joint_gtsp.par"], cwd='./GLKH', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        observer = Observer()
        observer.schedule(self.file_handler, "./GLKH/joint_gtsp.tour", recursive=False)
        observer.start()

        out, err = p.communicate()

        # save out to log file
        out = out.decode('utf-8')
        with open("GLKH/joint_gtsp.log", mode='w') as file:
            file.write(out)

        observer.stop()
        observer.join()


    def solve(self, time_limit, num_samples):
        self.motions.start_time = time.time()
        self.ik_solutions, self.ik_i_to_pose_i, self.pose_i_to_ik_i = self.sample_iks(num_samples, self.poses)
        self.motions.sampling_ik_time = time.time() - self.motions.start_time

        self.joint_distance_matrix = self.compute_joint_distance_matrix(self.ik_solutions, self.is_neighbour, self.poses, self.ik_i_to_pose_i)
        self.save_to_gtsplib("joint_gtsp", self.ik_solutions, self.joint_distance_matrix, self.pose_i_to_ik_i, time_limit)

        self.run_glkh()

        self.motion_performance()
        self.save_motion("gtsp")
        
