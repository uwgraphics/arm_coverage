from datetime import datetime
import numpy as np
import transformations as T

from static_variables import *

def get_curr_time():
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def check_pose(pose1, pose2, tolerances):
    distance = 0
    if tolerances[0] < 1e-6 and tolerances[1] < 1e-6 and tolerances[2] < 1e-6:
        distance = np.linalg.norm(pose1[0:3] - pose2[0:3])
    elif tolerances[0] > 1e-6 and tolerances[1] > 1e-6 and tolerances[2] < 1e-6:
        d_in_w = np.array(pose2[0:3]) - np.array(pose1[0:3])
        # distance in goal pose frame
        d_in_g = np.dot(d_in_w, T.quaternion_matrix([pose2[6], pose2[3], pose2[4], pose2[5]])[:3, :3])

        x = np.abs(d_in_g[0])
        y = np.abs(d_in_g[1])
        
        if x > tolerances[0] or y > tolerances[1]:
            return False
        
        z = np.abs(d_in_g[2])
        distance = z
    else:
        raise ValueError("Not implemented yet", tolerances)

    angle = 0

    if tolerances[3] < 1e-6 and tolerances[4] < 1e-6 and tolerances[5] < 1e-6:
        angle = T.quaternion_difference([pose1[6], pose1[3], pose1[4], pose1[5]], [pose2[6], pose2[3], pose2[4], pose2[5]])
    elif tolerances[5] > 100 and tolerances[3] < 1e-6 and tolerances[4] < 1e-6:
        angle = T.quaternion_to_angle_about_xy([pose1[6], pose1[3], pose1[4], pose1[5]], [pose2[6], pose2[3], pose2[4], pose2[5]])
    elif tolerances[3] > 1e-6 and tolerances[4] > 1e-6 and tolerances[5] > 100:
        q1 = [pose1[6], pose1[3], pose1[4], pose1[5]]
        q2 = [pose2[6], pose2[3], pose2[4], pose2[5]]
        scaled_axis = T.quaternion_to_scaledAxis(T.quaternion_dispQ(q1, q2))
        rx = scaled_axis[0]
        ry = scaled_axis[1]

        if np.abs(rx) > tolerances[3] or np.abs(ry) > tolerances[4]:
            return False
        
        angle = 0.0

    else:
        raise ValueError("Not implemented yet", tolerances)

    if  distance > 0.001:
        return False

    if angle > 0.01:
        return False
    
    return True

class Motions: 
    def __init__(self):
        self.n = 0
        self.times = []
        self.motions = []
        self.c_paths = []

        self.ds = []
        self.num_reconfigs = []
        self.ds_plot = []
        self.num_reconfigs_plot = []
        self.times_plot = []

        self.best_idx = None
        self.best_motion = None
        self.best_c_path = None

        self.sampling_ik_time = 0
        self.start_time = 0