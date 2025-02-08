from datetime import datetime
import numpy as np
import transformations as T
import ctypes
import yaml
import os
import random
import math
import sys

from static_variables import *
from robot_kinematics import RobotKinematics

pwd = os.getcwd()
sys.path.append(pwd)
print("sys.path: ", sys.path)
from wrappers.python_wrapper import RelaxedIKRust
from utils import *

class Robot:

    def __init__(self, robot_name):
        pwd = os.getcwd()
        
        setting_file_path =  pwd +'/configs/' + SETTING_FILE_PATHS[robot_name]
        # print("setting_file_path: ", setting_file_path)
        setting_file = open(setting_file_path, 'r')
        settings = yaml.load(setting_file, Loader=yaml.FullLoader)
       
        urdf_file = open( pwd +'/configs/urdfs/' + settings["urdf"], 'rb')
        urdf_string = urdf_file.read()

        self.rk = RobotKinematics(robot_name, urdf_string, BASE_LINKS[robot_name], TIP_LINKS[robot_name])

        self.n_joints = self.rk._num_jnts

        self.joint_names = [self.rk._robot_name + "-" + jn for jn in self.rk._joint_names]

        self.joint_type = self.rk._joint_type
        self.theta_lower_limits = np.array(self.rk._joint_lower_limits) 
        self.theta_upper_limits = np.array(self.rk._joint_upper_limits)  
        self.max_joint_velocity = np.array(self.rk._joint_vel_limits)

        curr_path = os.getcwd()   
        self.relaxed_ik = RelaxedIKRust(setting_file_path)

        os.chdir(curr_path)


    def reach_with_trac_ik(self, target_pose, seed_state = []):
        if len(seed_state) == 0:
            seed_state = self.generate_random_angle()
        trac_ik_sol = self.trac_ik.get_ik( seed_state,
                            target_pose[0],
                            target_pose[1],
                            target_pose[2],
                            target_pose[3],
                            target_pose[4],
                            target_pose[5],
                            target_pose[6])
        if trac_ik_sol is None:
            return []

        ik_solutions = list(trac_ik_sol)

        curr_pose = self.rk.forward_position_kinematics(ik_solutions)
        if check_pose(curr_pose, target_pose):
            return ik_solutions
        else:
            return []
        
    def reach_with_spot_ik(self, target_pose, start_config = []):

        res = self.spot_ik.solve(target_pose[0], target_pose[1], target_pose[2], target_pose[3], target_pose[4], target_pose[5], target_pose[6], start_config)

        return res

    def reset_relaxed_ik(self, start_config):
        x = (ctypes.c_double * len(start_config))()
        for i in range(len(start_config)):
            x[i] = start_config[i]
        self.relaxed_ik.reset(x)

    def reach_with_relaxed_ik(self, target_pose, constrain_velocity, tolerances, start_config = [],  test = False, max_iter = 10 ):

        if len(start_config) > 0:
            x = (ctypes.c_double * len(start_config))()
            for i in range(len(start_config)):
                x[i] = start_config[i]
            self.relaxed_ik.reset(x)
        else:
            x = (ctypes.c_double * self.n_joints)()
            random_config = self.generate_random_angle()
            assert len(random_config) == self.n_joints
            for i in range(self.n_joints):
                x[i] = random_config[i]
            self.relaxed_ik.reset(x)

        if test:
            print("start_config: {}".format(start_config))
            print("target_pose: {}".format(target_pose))

        counter = 0
        while counter < max_iter:
            ik_solutions = self.relaxed_ik.solve_position(target_pose[0:3], target_pose[3:7], tolerances, constrain_velocity)

            curr_pose = self.rk.forward_position_kinematics(ik_solutions)

            if test:
                print("ik_solutions: {} curr_pose: {}".format(ik_solutions, curr_pose))
            
            if check_pose(curr_pose, target_pose, tolerances):
                return ik_solutions
            
            counter += 1

        return []

    def generate_random_placement(self):
        return [random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-math.pi, math.pi)]

    def generate_random_angle(self):
        theta = np.zeros(self.n_joints)
        for i in range(self.n_joints):
            # theta[i] =  self.theta_lower_limits[i] + random.uniform(0, self.theta_upper_limits[i] - self.theta_lower_limits[i])
            if self.joint_type[i] == 'revolute':
                theta[i] =  self.theta_lower_limits[i] + random.uniform(0, self.theta_upper_limits[i] - self.theta_lower_limits[i])
            elif self.joint_type[i] == 'continuous':
                theta[i] = random.uniform(-4 * math.pi, 4 * math.pi)
            elif self.joint_type[i] == 'prismatic':
                theta[i] =  self.theta_lower_limits[i] + random.uniform(0, self.theta_upper_limits[i] - self.theta_lower_limits[i])
            else:
                print("Unknown joint type: {}".format(self.joint_type[i]))
        return theta
