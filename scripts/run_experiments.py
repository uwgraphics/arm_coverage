import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from datetime import datetime
import csv
import glob 
from utils import *
import time
import static_variables

from joint_gtsp import JointGTSP
from cartesian_tsp import CartesianTSP
from hrchy_gtsp import HrchyGTSP

def write_performance(writer, task, solver, motions):
    best_idx = motions.best_idx
    writer.writerow({'task': task, 'solver': solver, 
                    'timestamp': datetime.now().strftime("%Y%m%d-%H%M%S"),
                    'joint_movement': motions.ds[best_idx], 
                    'num_reconfigs': motions.num_reconfigs[best_idx], 
                    'num_points_on_surface': len(motions.motions[best_idx])})

if __name__ == '__main__':
    
    # glob all files under ee_poses
    files = glob.glob("./surfaces/*.csv")

    performance_file_name = "performance/performance_{}.csv".format(datetime.now().strftime("%Y%m%d-%H%M%S"))

    performance_file = open(performance_file_name, mode='w')

    fields = ['task', 'solver', 'timestamp',  'joint_movement', 'num_reconfigs', 'num_points_on_surface']

    with performance_file:
        writer = csv.DictWriter(performance_file, fieldnames=fields)
        writer.writeheader()

        for file_name in files:

            num_samples = 100
            time_limit = 100
            print("file_name: ", file_name)
            task = file_name.split('/')[-1].split('.')[0].split('_')[-1]
            print("task: ", task)
            if task == "wok":
                tolerances = [0.05, 0.05, 0., 0., 0., 999.]
            if task == "wokreal":
                tolerances = [0., 0., 0., 0., 0., 999.]
            elif task == "floor":
                tolerances = [0., 0., 0., 0., 0., 0.]
            else:
                print("task not found")
            
            print("=====>  Running Cartesian-TSP")
            cartesian_tsp_solver = CartesianTSP(file_name, tolerances)
            cartesian_tsp_solver.solve(time_limit, num_samples)
            write_performance(writer, task, "cartesian_tsp", cartesian_tsp_solver.motions)
            performance_file.flush()

            print("=====>  Running Hrchy-Joint-GTSP")
            hrchy_gtsp_solver = HrchyGTSP(file_name, tolerances)
            hrchy_gtsp_solver.solve(time_limit, num_samples)
            write_performance(writer, task, "hrchy_gtsp", hrchy_gtsp_solver.motions)
            performance_file.flush()

            print("=====>  Running Joint-GTSP")
            joint_gtsp_solver = JointGTSP(file_name, tolerances)
            joint_gtsp_solver.solve(time_limit, num_samples)
            write_performance(writer, task, "joint_gtsp", joint_gtsp_solver.motions)
            performance_file.flush()

    performance_file.close()


        
