import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import time
import subprocess
from utils import *
import time
from watchdog.observers import Observer
import static_variables
from solver_base import SolverBase


class CartesianTSP(SolverBase):
    def __init__(self, file_name, tolerances):
        SolverBase.__init__(self, file_name, tolerances)
        self.routines = []
        self.motions = Motions()

        f = open("./GLKH/cartesian_tsp.tour", 'a+') 
        f.close()

        self.file_handler = self.TourFileChangeHandler(self, "./GLKH/cartesian_tsp.tour",  self.file_callback)

    def iklink(self, tsp_routine):

        start = time.time()

        assert len(tsp_routine) == len(self.poses)

        pre_num_reconfigs = np.zeros(len(self.pose_i_to_ik_i[tsp_routine[0]]))
        pre_ds = np.zeros(len(self.pose_i_to_ik_i[tsp_routine[0]]))
        pre_idx = np.zeros(len(self.ik_solutions))

        for x in range(1, len(tsp_routine)):
            ds = []
            num_reconfigs = []

            ik_indices1 = self.pose_i_to_ik_i[tsp_routine[x]]
            ik_indices2 = self.pose_i_to_ik_i[tsp_routine[x-1]]

            for y1 in range(len(ik_indices1)):
                ds.append(1e9)
                num_reconfigs.append(1e9)
                ik1 = np.array(self.ik_solutions[ik_indices1[y1]])
                for y2 in range(len(ik_indices2)):
                    ik2 = np.array(self.ik_solutions[ik_indices2[y2]])
                    d = np.linalg.norm(ik1 - ik2)
                    num_reconfig = 1 if self.needs_reconfig(ik2, ik1, self.poses[tsp_routine[x-1]], self.poses[tsp_routine[x]]) else 0

                    d += pre_ds[y2]
                    num_reconfig += pre_num_reconfigs[y2]
                    
                    if num_reconfig < num_reconfigs[-1] or (num_reconfig == num_reconfigs[-1] and d < ds[-1]):
                        ds[-1] = d
                        num_reconfigs[-1] = num_reconfig
                        pre_idx[ik_indices1[y1]] = int(ik_indices2[y2])

            pre_ds = ds
            pre_num_reconfigs = num_reconfigs

        min_d = 1e9
        min_num_reconfig = 1e9
        min_idx = -1

        last_ik_indices = self.pose_i_to_ik_i[tsp_routine[-1]]

        assert len(pre_ds) == len(pre_num_reconfigs)
        assert len(pre_ds) == len(last_ik_indices)

        for y in range(len(pre_ds)):
            if pre_ds[y] < min_d or (pre_ds[y] == min_d and pre_num_reconfigs[y] < min_num_reconfig):
                min_d = pre_ds[y]
                min_num_reconfig = pre_num_reconfigs[y]
                min_idx = last_ik_indices[y]

        # print("min num reconfig: ", min_num_reconfig, "Min d: ", min_d, "Min idx: ", min_idx)
        motion = [self.ik_solutions[min_idx]]
        i = min_idx
        while len(motion) < len(self.poses):
            i = int(pre_idx[i])
            motion.append(self.ik_solutions[i])

        # reverse routine
        motion = motion[::-1]

        return motion

    def file_callback(self, routine):
        print("Got a cartesian path")
        if routine not in self.routines:
            self.routines.append(routine)
            motion = self.iklink(routine)
            c_path = [self.poses[i] for i in routine]
            self.motions.c_paths.append(c_path)
            self.motions.motions.append(motion)
            self.motions.times.append(time.time())
            self.motions.n += 1

    def save_to_tsplib(self, name, distance_matrix, time_limit):
        with open("./GLKH/{}.par".format(name), mode='w') as file:
            file.write("PROBLEM_FILE = {}.tsp\n".format(name))
            file.write("OUTPUT_TOUR_FILE = {}.tour\n".format(name))
            file.write("RUNS = 1\n")
            file.write("PRECISION = 1\n")
            file.write("TRACE_LEVEL = 1\n")
            file.write("MAX_TRIALS = 1000000\n")
            file.write("CANDIDATE_SET_TYPE = POPMUSIC\n")
            file.write("TIME_LIMIT = {}\n".format(time_limit))

        with open("./GLKH/{}.tsp".format(name), mode='w') as file:
            file.write("NAME: {}\n".format(name))
            file.write("TYPE: TSP\n")
            file.write("COMMENT: Nothing\n")

            n = len(self.poses)+1
            file.write("DIMENSION: {}\n".format(n))
            file.write("EDGE_WEIGHT_TYPE: EXPLICIT\n")
            file.write("EDGE_WEIGHT_FORMAT: UPPER_ROW\n")
            file.write("EDGE_WEIGHT_SECTION\n")
            for i in range(n):
                row = ""
                for j in range(i+1, n):
                    row += "{} ".format(int(distance_matrix[i, j]))
                row += "\n"
                file.write(row)
            file.write("EOF\n")

    def run_lkh(self):
        p = subprocess.Popen(["./LKH", "./cartesian_tsp.par"], cwd='./GLKH', stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        observer = Observer()
        observer.schedule(self.file_handler, "./GLKH/cartesian_tsp.tour", recursive=False)
        observer.start()

        out, err = p.communicate()

        # save out to log file
        out = out.decode('utf-8')
        with open("GLKH/cartesian_tsp.log", mode='w') as file:
            file.write(out)

        observer.stop()
        observer.join()

        # time.sleep(1)

    def solve(self, time_limit, num_samples):
        self.motions.start_time = time.time()
        self.ik_solutions, self.ik_i_to_pose_i, self.pose_i_to_ik_i = self.sample_iks(num_samples, self.poses)
        self.motions.sampling_ik_time = time.time() - self.motions.start_time

        self.compute_cartesian_distance_matrix()

        self.save_to_tsplib("cartesian_tsp", self.cartesian_distance_matrix, time_limit)

        self.run_lkh()

        self.motion_performance()
        self.save_motion("tsp_iklink")
