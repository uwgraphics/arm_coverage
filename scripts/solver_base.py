from robot import Robot
import sys
import os
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import time
from datetime import datetime
import csv
import transformations as T
import pyvista as pv
from sklearn.cluster import DBSCAN
from utils import *
import time
from watchdog.events import FileSystemEventHandler

class SolverBase:
    def __init__(self, file_name, tolerances):
        robot_name = file_name.split("/")[-1].split("_")[0]
        self.robot_name = robot_name
        self.robot = Robot(robot_name)
        self.ik_solutions = []
        self.ik_i_to_pose_i = []
        self.pose_i_to_ik_i = []
        self.tolerances = tolerances

        self.ik_solutions_h = []
        self.ik_i_h_to_pose_i_h = []
        self.pose_i_to_ik_i_h = []

        with open(file_name, mode='r') as file:
            reader = csv.DictReader(file)
            poses = []
            for row in reader:
                poses.append([float(row['target-POS_X']), float(row['target-POS_Y']), float(row['target-POS_Z']), float(row['target-ROT_X']), float(row['target-ROT_Y']), float(row['target-ROT_Z']), float(row['target-ROT_W'])])
        self.poses = poses

        points = np.array(self.poses)[:, :3]
        cloud = pv.PolyData(points)
        
        surf = cloud.delaunay_2d()
        # surf = cloud.delaunay_3d()

        self.is_neighbour = np.zeros((len(self.poses), len(self.poses)), dtype=bool)

        i = 0
        while i < len(surf.faces):
            n = surf.faces[i]
            assert n == 3
            a = surf.faces[i+1]
            b = surf.faces[i+2]
            c = surf.faces[i+3]
            self.is_neighbour[a, b] = True
            self.is_neighbour[b, a] = True
            self.is_neighbour[b, c] = True
            self.is_neighbour[c, b] = True
            self.is_neighbour[a, c] = True
            self.is_neighbour[c, a] = True
            i += 4

        current_dir = os.getcwd()
        self.GLKH_dir = os.path.join(current_dir, "GLKH")
        assert os.path.isdir(self.GLKH_dir), f"Folder {self.GLKH_dir} does not exist"

    class TourFileChangeHandler(FileSystemEventHandler):
        def __init__(self, parent, file_path, callback):
            self.parent = parent
            self.file_path = file_path
            self.callback = callback

        def on_modified(self, event):
            if event.src_path == self.file_path:
                with open(event.src_path, mode='r') as file:
                    lines = file.readlines()
                    tour = []
                    start = False
                    for line in lines:
                        if line.startswith("TOUR_SECTION"):
                            start = True
                            continue
                        if start:
                            if line.startswith("-1"):
                                break
                            tour.append(int(line))

                gap_idx = 0

                dummy_node = np.max(tour)

                for i in range(0, len(tour)):
                    if tour[i] == dummy_node:
                        gap_idx = i+1
                        if gap_idx == len(tour):
                            gap_idx = 0
                        break
                
                if gap_idx == 0:
                    routine = tour[:-1]
                else:
                    routine = tour[gap_idx:] + tour[:gap_idx-1]
            
                self.callback([x-1 for x in routine])

    def distance_between_points(self, a, b):
        return np.linalg.norm(np.array(a) - np.array(b))

    def motion_performance(self):
        motions = self.motions

        ds = []
        num_reconfigs = []

        for i in range(motions.n):
            motion = motions.motions[i]
            c_path = motions.c_paths[i]

            d = 0
            num_reconfig = 0
            for j in range(len(motion)-1):
                ik1 = motion[j]
                ik2 = motion[j+1]
                pose1 = c_path[j]
                pose2 = c_path[j+1]
                if self.needs_reconfig(ik1, ik2, pose1, pose2):
                    num_reconfig += 1
                else:
                    d += self.distance_between_points(ik1, ik2)

            ds.append(d)
            num_reconfigs.append(num_reconfig)
                
        ds.append(ds[-1])
        num_reconfigs.append(num_reconfigs[-1])

        best_num_reconfig = num_reconfigs[0]
        best_d = ds[0]
        best_idx = 0

        for i in range(1, len(ds)):
            if num_reconfigs[i] < best_num_reconfig or (num_reconfigs[i] == best_num_reconfig and ds[i] < best_d):
                best_num_reconfig = num_reconfigs[i]
                best_d = ds[i]
                best_idx = i

        motions.ds = ds
        motions.num_reconfigs = num_reconfigs
        motions.best_idx = best_idx
        motions.best_motion = motions.motions[best_idx]
        motions.best_c_path = motions.c_paths[best_idx]
    
        print("Number of reconfigs:", motions.num_reconfigs[best_idx])
        print("Joint movement:", motions.ds[best_idx])

    def sample_iks(self, num_samples, poses):
        print("Sampling IK: self.tolerances: ", self.tolerances)
        self.num_ik_per_ee_pose = num_samples
        ik_solutions = []
        ik_i_to_pose_i = []
        pose_i_to_ik_i = []
        for (pose_i, pose) in enumerate(poses):
            iks = []
            for _ in range(num_samples):
                # print("Sampling for pose: ", pose)
                ik = self.robot.reach_with_relaxed_ik(pose, False, tolerances = self.tolerances)
                counter = 0
                while len(ik) == 0:
                    ik = self.robot.reach_with_relaxed_ik(pose, False, tolerances = self.tolerances)
                    counter += 1
                    if counter > 100:
                        raise Exception("Failed to find IK for pose: ", pose)
                iks.append(ik)
            
            # cluster using DBSCAN
            ik_array = np.array(iks)
            db = DBSCAN(eps=0.1, min_samples=2).fit(ik_array)
            labels = db.labels_

            centers = []
            uesd_labels = []
            for i in range(len(labels)):
                if labels[i] == -1 or labels[i] not in uesd_labels:
                    centers.append(ik_array[i])
                    uesd_labels.append(labels[i])
                
            pose_i_to_ik_i.append([])
            for center in centers:
                ik_solutions.append(center)
                ik_i_to_pose_i.append(pose_i)
                pose_i_to_ik_i[pose_i].append(int(len(ik_solutions) - 1))

        return ik_solutions, ik_i_to_pose_i, pose_i_to_ik_i

    def interpolate_ik(self, l_ik, r_ik, r):
        return l_ik + (r_ik - l_ik) * r

    def point_to_line(self, p, l1, l2):
        return np.linalg.norm(np.cross(np.array(l2) - np.array(l1), np.array(l1) - np.array(p))) / np.linalg.norm(np.array(l2) - np.array(l1))
        
    def max_joint_move(self, l_ik, r_ik):
        return np.max(np.abs(np.array(l_ik) - np.array(r_ik)))

    def needs_reconfig(self, l_ik, r_ik, l_ee_pose, r_ee_pose, factor = 1.0):
        posi_d = np.linalg.norm(np.array(l_ee_pose[0:3]) - np.array(r_ee_pose[0:3]))
        q1 = [l_ee_pose[6], l_ee_pose[3], l_ee_pose[4], l_ee_pose[5]]
        q2 = [r_ee_pose[6], r_ee_pose[3], r_ee_pose[4], r_ee_pose[5]]
        angle_d = T.angle_between_quaternion_z_axes(q1, q2)

        if posi_d > 0.2 * factor or angle_d > 0.5 * factor:
            return True

        if self.max_joint_move(l_ik, r_ik) > 1.0 * factor:
            return True

        l1 = l_ee_pose[:3]
        l2 = r_ee_pose[:3]
        q1 = l_ee_pose[3:]
        q2 = r_ee_pose[3:]

        n = 2
        for i in range(1, n):
            r = i / n
            ik = self.interpolate_ik(np.array(l_ik), np.array(r_ik), r)

            ee_pose = self.robot.rk.forward_position_kinematics(ik)

            p = ee_pose[:3]
            q = ee_pose[3:]

            d = self.point_to_line(p, l1, l2)

            if d > 0.1 * factor:
                return True

        return False

    def compute_joint_distance(self, l_ik, r_ik, l_ee_pose, r_ee_pose,  distance_factor = 1.0):
        if self.needs_reconfig(l_ik, r_ik, l_ee_pose, r_ee_pose, factor = distance_factor):
            return SEMI_LARGE_WEIGHT
        else:
            return int(np.linalg.norm(np.array(l_ik) - np.array(r_ik)) * 1e4)

    def distance_between_poses(self, l_ee_pose, r_ee_pose):
        posi_d = np.linalg.norm(np.array(l_ee_pose[0:3]) - np.array(r_ee_pose[0:3]))
        q1 = [l_ee_pose[6], l_ee_pose[3], l_ee_pose[4], l_ee_pose[5]]
        q2 = [r_ee_pose[6], r_ee_pose[3], r_ee_pose[4], r_ee_pose[5]]
        angle_d = T.angle_between_quaternion_z_axes(q1, q2)
        return posi_d + 0.1 * angle_d

    def compute_cartesian_distance_matrix(self):
        n = len(self.poses)
        distance_matrix = np.zeros((n+1, n+1), dtype=np.int32)
        for i in range(n):
            for j in range(i+1, n):
                if i == j:
                    distance_matrix[i, j] = 0
                elif self.is_neighbour[i, j]:
                    distance_matrix[i, j] = int(self.distance_between_poses(self.poses[i], self.poses[j]) * 1e4)
                else:
                    distance_matrix[i, j] = LARGE_WEIGHT

        distance_matrix[n, n] = LARGE_WEIGHT

        self.cartesian_distance_matrix = distance_matrix

    def compute_joint_distance_matrix(self, ik_solutions, is_neighbour, poses, ik_i_to_pose_i, distance_factor = 1.0):
        n = len(ik_solutions)
        distance_matrix = np.zeros((n+1, n+1), dtype=np.int32)

        for i in range(n):
            start = time.time()
            l_pose_i = ik_i_to_pose_i[i]
            l_ik = ik_solutions[i]
            l_ee_pose = poses[l_pose_i]

            for j in range(i+1, n):
                r_pose_i = ik_i_to_pose_i[j]

                if i == j:
                    distance_matrix[i, j] = 0
                elif is_neighbour[l_pose_i, r_pose_i]:
                    r_ik = ik_solutions[j]
                    r_ee_pose = poses[r_pose_i] 
                    distance_matrix[i, j] = self.compute_joint_distance(l_ik, r_ik, l_ee_pose, r_ee_pose,  distance_factor = distance_factor)
                else:
                    distance_matrix[i, j] = LARGE_WEIGHT

        distance_matrix[n, n] = LARGE_WEIGHT

        return distance_matrix

    def gtsp_routine_distance(self, routine):
        distance = 0.
        for i in range(len(routine)-1):
            distance += np.linalg.norm(np.array(self.ik_solutions[routine[i]]) - np.array(self.ik_solutions[routine[i+1]]))
        return distance

    def cal_total_distance(self, routine):
        num_points, = routine.shape
        return sum([self.joint_distance_matrix[routine[i], routine[(i + 1)]] for i in range(num_points-1)])

    def save_motion(self, tag):

        motions = self.motions

        motion = motions.best_motion
        c_path = motions.best_c_path

        assert len(motion) == len(c_path) 

        file_name = self.robot_name + "_" + tag + "_"+ datetime.now().strftime("%Y%m%d_%H%M%S") + ".csv"

        joint_names = [self.robot_name + "-" + jn for jn in self.robot.rk._joint_names]
        fieldnames = ['time'] + joint_names + ['target-POS_X', 'target-POS_Y', 'target-POS_Z', 'target-ROT_X', 'target-ROT_Y', 'target-ROT_Z', 'target-ROT_W']

        # get pwd
        pwd = os.getcwd()
        with open("{}/motions/{}".format(pwd, file_name), mode='w') as file:

            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            time = 0.0
            for i in range(len(motion)):
                
                row = [
                    time,
                    *motion[i],
                    *c_path[i]
                ]

                writer.writerow(dict(zip(fieldnames, row)))
                time += 1.
            
            print("Motion saved to: ", file_name)

    def save_ik_samplings(self):
        file_name = self.robot_name + "_" + datetime.now().strftime("%Y%m%d_%H%M%S") + ".csv"

        joint_names = [self.robot_name + "-" + jn for jn in self.robot.rk._joint_names]
        fieldnames = ['time'] + joint_names + ['target-POS_X', 'target-POS_Y', 'target-POS_Z', 'target-ROT_X', 'target-ROT_Y', 'target-ROT_Z', 'target-ROT_W']

        # get pwd
        pwd = os.getcwd()
        with open("{}/motions/{}".format(pwd, file_name), mode='w') as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            time = 0.0
            for i in range(len(self.ik_solutions)):
                
                row = [
                    time,
                    *self.ik_solutions[i],
                    *self.poses[self.ik_i_to_pose_i[i]]
                ]

                writer.writerow(dict(zip(fieldnames, row)))
                time += 1. / 30.
            
            print("Saved to: ", file_name)

    def save_to_gtsplib(self, name, ik_solutions, distance_matrix, pose_i_to_ik_i, time_limit=60):
        with open("./GLKH/{}.par".format(name), mode='w') as file:
            file.write("PROBLEM_FILE = {}/{}.gtsp\n".format(self.GLKH_dir, name))
            file.write("ASCENT_CANDIDATES = 500\n")
            file.write("CANDIDATE_SET_TYPE = POPMUSIC\n")
            file.write("POPMUSIC_SAMPLE_SIZE = 100\n")
            file.write("POPMUSIC_MAX_NEIGHBORS = 30\n")
            file.write("POPMUSIC_TRIALS = 0\n")
            file.write("INITIAL_PERIOD = 1000\n")
            file.write("MAX_CANDIDATES = 30\n")
            file.write("MAX_TRIALS = 1000000\n")
            file.write("OUTPUT_TOUR_FILE = {}.tour\n".format(name))
            file.write("RUNS = 1\n")
            file.write("POPULATION_SIZE = 1\n")
            file.write("PRECISION = 1\n")
            file.write("TRACE_LEVEL = 1\n")
            file.write("TIME_LIMIT = {}\n".format(time_limit))

        with open("./GLKH/{}.gtsp".format(name), mode='w') as file:
            file.write("NAME: {}\n".format(name))
            file.write("TYPE: GTSP\n")
            file.write("COMMENT: Nothing\n")

            n = len(ik_solutions) + 1
            k = len(pose_i_to_ik_i) + 1
            file.write("DIMENSION: {}\n".format(n))
            file.write("GTSP_SETS: {}\n".format(k))
            file.write("EDGE_DATA_FORMAT: EDGE_LIST\n")
            file.write("EDGE_DATA_SECTION\n")
            for i in range(n):
                for j in range(i+1, n):
                    if i != j and distance_matrix[i,j]!=LARGE_WEIGHT:
                        file.write("{} {} {}\n".format(i+1, j+1, distance_matrix[i, j]))

            file.write("GTSP_SET_SECTION\n")
            for i in range(len(pose_i_to_ik_i)):
                row = "{} ".format(i+1)
                for j in pose_i_to_ik_i[i]:
                    row += "{} ".format(j+1)
                row += "-1\n"
                file.write(row)
            
            file.write("{} {} -1\n".format(k, n))
            file.write("EOF\n")
