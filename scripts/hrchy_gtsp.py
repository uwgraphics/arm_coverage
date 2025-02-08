import sys
import os
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import time
from datetime import datetime
import csv
import subprocess
import pyvista as pv
from sklearn.cluster import DBSCAN, AffinityPropagation
from utils import *
import time
from watchdog.observers import Observer
import static_variables
from solver_base import SolverBase

class HrchyGTSP(SolverBase):
    def __init__(self, file_name, tolerances):
        SolverBase.__init__(self, file_name, tolerances)
        self.got_new_iks = False

        self.ik_solutions = []
        self.ik_i_to_pose_i = []
        self.pose_i_to_ik_i = []
        for i in range(len(self.poses)):
            self.pose_i_to_ik_i.append([])

        self.ik_solutions_h = []
        self.propergated_ik_h = []
        self.nearby_iks = []
        self.ik_i_h_to_pose_i = []
        self.ik_i_h_to_pose_i_h = []
        self.pose_i_h_to_ik_i_h = []

        self.h_routines = []
        self.propergated_index = 0
        self.motions = Motions()

        f = open("./GLKH/joint_gtsp_high_level.tour", 'a+') 
        f.close()

        f = open("./GLKH/joint_gtsp_low_level.tour", 'a+') 
        f.close()

        self.h_file_handler = self.TourFileChangeHandler(self, "./GLKH/joint_gtsp_high_level.tour", self.high_level_file_callback)
        self.low_level_file_handler = self.TourFileChangeHandler(self, "./GLKH/joint_gtsp_low_level.tour", self.low_level_file_callback)
        
        # cluster points using AffinityPropagation
        points = np.array(self.poses)[:, :3]

        cloud = pv.PolyData(points)
        surf = cloud.delaunay_2d()

        self.pose_distance_matrix = np.zeros((len(self.poses), len(self.poses)))
        for i in range(len(self.poses)):
            for j in range(i+1, len(self.poses)):
                d = self.distance_between_poses(self.poses[i], self.poses[j])
                # d = np.linalg.norm(points[i] - points[j])

                # check if d is NaN
                if np.isnan(d):
                    print("Pose i: ", self.poses[i])
                    print("Pose j: ", self.poses[j])
                    raise Exception("Distance between poses is NaN")

                self.pose_distance_matrix[i, j] = d
                self.pose_distance_matrix[j, i] = d

        # make sure distance matrix doesn't contain NaN
        assert np.isnan(self.pose_distance_matrix).any() == False, "Distance matrix contains NaN"

        self.poses_h = []
        # clustering = AffinityPropagation().fit(points)
        # clustering = AffinityPropagation( preference= - 2).fit(distance_matrix)
        clustering = AffinityPropagation().fit(self.pose_distance_matrix)
        self.clustered_pose_centers_indices = clustering.cluster_centers_indices_

        for i in self.clustered_pose_centers_indices:
            self.poses_h.append(self.poses[i])

        self.clustered_pose_nearby_indices = []
        for i in range(len(self.clustered_pose_centers_indices)):
            self.clustered_pose_nearby_indices.append([])
            for j in range(len(clustering.labels_)):
                if clustering.labels_[j] == i and j not in self.clustered_pose_centers_indices:
                    self.clustered_pose_nearby_indices[-1].append(j)

        cloud = pv.PolyData([points[i] for i in self.clustered_pose_centers_indices])
        surf_h = cloud.delaunay_2d()
        # surf_h.plot(show_edges=True)

        plotter = pv.Plotter() 
        # translate surf_h to avoid overlapping
        plotter.add_mesh(surf, show_edges=True, color=[0.1, 0.6, 0.6])
        plotter.add_mesh(cloud, point_size=10, color='red')

        # plotter.show()

        self.is_neighbour_h = np.zeros((len(self.clustered_pose_centers_indices), len(self.clustered_pose_centers_indices)), dtype=bool)

        i = 0
        while i < len(surf_h.faces):
            n = surf_h.faces[i]
            assert n == 3
            a = surf_h.faces[i+1]
            b = surf_h.faces[i+2]
            c = surf_h.faces[i+3]
            self.is_neighbour_h[a, b] = True
            self.is_neighbour_h[b, a] = True
            self.is_neighbour_h[b, c] = True
            self.is_neighbour_h[c, b] = True
            self.is_neighbour_h[a, c] = True
            self.is_neighbour_h[c, a] = True
            i += 4

        # print(self.is_neighbour_h)

    def save_h_routines(self):
        routine = self.h_routines[-1]
        pwd = os.getcwd()

        joint_names = [self.robot_name + "-" + jn for jn in self.robot.rk._joint_names]
        fieldnames = ['time'] + joint_names + ['target-POS_X', 'target-POS_Y', 'target-POS_Z', 'target-ROT_X', 'target-ROT_Y', 'target-ROT_Z', 'target-ROT_W']

        file_name = self.robot_name + "_h_"+ datetime.now().strftime("%Y%m%d_%H%M%S") + ".csv"
        with open("{}/motions/{}".format(pwd, file_name), mode='w') as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            time = 0.0
            for i in range(len(routine)):
                
                row = [
                    time,
                    *self.ik_solutions_h[routine[i]],
                    *self.poses_h[self.ik_i_h_to_pose_i_h[routine[i]]]
                ]

                writer.writerow(dict(zip(fieldnames, row)))
                time += 1.
            
            print("Motion saved to: ", file_name)


    def high_level_file_callback(self, routine):
        print("Got new high-level routine.")
        if routine not in self.h_routines:
            self.h_routines.append(routine)
            # self.save_h_routines()
            self.got_new_iks = True
        
    def low_level_file_callback(self, routine):
        print("Got new low-level routine: ", routine)
        j_motion = [self.ik_solutions[i] for i in routine]
        c_path = [self.poses[self.ik_i_to_pose_i[i]] for i in routine]
        self.motions.c_paths.append(c_path)
        self.motions.motions.append(j_motion)
        self.motions.times.append(time.time())
        self.motions.n += 1

    def run_glkh_h(self, time_limit):
        p = subprocess.Popen(["./GLKH", "./joint_gtsp_high_level.par"], cwd='./GLKH', stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        observer = Observer()
        observer.schedule(self.h_file_handler, "./GLKH/joint_gtsp_high_level.tour", recursive=False)
        observer.start()

        no_solution = True

        while time.time() - self.motions.start_time < time_limit or no_solution:
            if self.got_new_iks:
                self.got_new_iks = False
                if self.ik_propergate():
                    self.joint_distance_matrix = self.compute_joint_distance_matrix(self.ik_solutions, self.is_neighbour, self.poses, self.ik_i_to_pose_i)

                    self.save_to_gtsplib("joint_gtsp_low_level", self.ik_solutions, self.joint_distance_matrix, self.pose_i_to_ik_i, time_limit/5)
                    
                    self.run_glkh_low_level()

                    no_solution = False
           
        out, err = p.communicate()

        # save out to log file
        out = out.decode('utf-8')
        with open("GLKH/joint_gtsp_high_level.log", mode='w') as file:
            file.write(out)
        
        observer.stop()
        observer.join()

    def run_glkh_low_level(self):
        p = subprocess.Popen(["./GLKH", "./joint_gtsp_low_level.par"], cwd='./GLKH', stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        observer = Observer()
        observer.schedule(self.low_level_file_handler, "./GLKH/joint_gtsp_low_level.tour", recursive=False)
        observer.start()

        out, err = p.communicate()

        # save out to log file
        out = out.decode('utf-8')
        with open("GLKH/joint_gtsp_low_level.log", mode='w') as file:
            file.write(out)

        observer.stop()
        observer.join()

    def ik_propergate(self):

        l = len(self.ik_solutions)

        while self.propergated_index < len(self.h_routines):
            routine = self.h_routines[self.propergated_index]
            self.propergated_index += 1

            pre_i = len(self.ik_solutions)
            for i in routine:
                
                if self.propergated_ik_h[i]:
                    continue
                
                self.propergated_ik_h[i] = True

                pose_i_h = self.ik_i_h_to_pose_i_h[i]
                pose_i = self.ik_i_h_to_pose_i[i]
                ik = self.ik_solutions_h[i]
                self.ik_solutions.append(ik)
                self.ik_i_to_pose_i.append(pose_i)
                self.pose_i_to_ik_i[pose_i].append(len(self.ik_solutions)-1)
                for j in self.clustered_pose_nearby_indices[pose_i_h]:
                    results = self.nearby_iks[i][j]
                    # assert len(results) > 0, "Failed to find IK for pose: {}, ik: {}".format(pose, ik)
                    self.ik_solutions.append(results)
                    self.ik_i_to_pose_i.append(j)
                    self.pose_i_to_ik_i[j].append(len(self.ik_solutions)-1)

        return l != len(self.ik_solutions)

    def sample_iks_h(self, num_iks):
        # print(self.clustered_pose_centers_indices)

        for i in range(len(self.clustered_pose_centers_indices)):
            center_pose_i = self.clustered_pose_centers_indices[i]
            pose = self.poses[center_pose_i]
            self.poses_h.append(pose)

            iks = []
            for j in range(num_iks):
                ik = self.robot.reach_with_relaxed_ik(pose, False,  self.tolerances)
                counter = 0
                while len(ik) == 0:
                    ik = self.robot.reach_with_relaxed_ik(pose, False, self.tolerances)
                    counter += 1
                    if counter > 1000:
                        raise Exception("Failed to find IK for pose: {}".format(pose))
                iks.append(ik)

            ik_array = np.array(iks)
            db = DBSCAN(eps=0.1, min_samples=2).fit(ik_array)
            labels = db.labels_
            centers = []
            used_labels = []
            for j in range(len(labels)):
                if labels[j] == -1 or labels[j] not in used_labels:
                    centers.append(ik_array[j])
                    used_labels.append(labels[j])   

            self.pose_i_h_to_ik_i_h.append([])
            for ik in centers:
                can_reach_to_all_nearby = True
                nearby_iks = {}
                for k in self.clustered_pose_nearby_indices[i]:
                    pose = self.poses[k]
                    results = self.robot.reach_with_relaxed_ik(pose, True,  self.tolerances, start_config=ik, max_iter=100)
                    nearby_iks[k] = results
                    if len(results) == 0:
                        can_reach_to_all_nearby = False
                        break

                if can_reach_to_all_nearby:
                    self.ik_solutions_h.append(ik)
                    self.propergated_ik_h.append(False)
                    self.nearby_iks.append(nearby_iks)
                    self.ik_i_h_to_pose_i_h.append(i)
                    self.ik_i_h_to_pose_i.append(center_pose_i)
                    self.pose_i_h_to_ik_i_h[i].append(len(self.ik_solutions_h)-1)

            assert len(self.pose_i_h_to_ik_i_h[i]) > 0, "Failed to find IK for pose: {}".format(pose)
        
    def solve(self, time_limit, num_samples):
        self.motions.start_time = time.time()
        self.sample_iks_h(num_samples)
        self.motions.sampling_ik_time = time.time() - self.motions.start_time
        
        self.joint_distance_matrix_h = self.compute_joint_distance_matrix(self.ik_solutions_h, self.is_neighbour_h, self.poses_h, self.ik_i_h_to_pose_i_h, distance_factor=3.0)

        self.save_to_gtsplib("joint_gtsp_high_level", self.ik_solutions_h, self.joint_distance_matrix_h, self.pose_i_h_to_ik_i_h, time_limit)

        self.run_glkh_h(time_limit)

        self.motion_performance()
        self.save_motion("hrchy_gtsp")
