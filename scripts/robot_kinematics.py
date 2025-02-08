
import numpy as np
import transformations as T
import math

from urdf_parser_py.urdf import URDF
from robot_kdl.kdl_parser import kdl_tree_from_urdf_model
import PyKDL

class RobotKinematics(object):
    def __init__(self, robot_name, urdfString, base_link, tip_link):
        
        self._robot = URDF.from_xml_string(urdfString)
        self._kdl_tree = kdl_tree_from_urdf_model(self._robot)
        self._base_link = base_link

        self._robot_name = robot_name
        
        self._tip_link = tip_link
        self._tip_frame = PyKDL.Frame()
        self._arm_chain = self._kdl_tree.getChain(self._base_link,
                                                  self._tip_link)

        self._joint_lower_limits = []
        self._joint_upper_limits = []
        self._joint_vel_limits = []
        self._joint_type = []
            
        self._joint_names = self._robot.get_chain(self._base_link , self._tip_link , links=False, fixed=False)
        self._num_jnts = self._arm_chain.getNrOfJoints()

        for jn in self._joint_names:
            for j in self._robot.joints:
                if j.name == jn:
                    if j.type == 'fixed':
                        continue
                    elif j.type == 'continuous':
                        self._joint_lower_limits.append(-4 * np.pi)
                        self._joint_upper_limits.append(4 * np.pi)
                        self._joint_vel_limits.append(j.limit.velocity)
                        self._joint_type.append('continuous')
                        # self._joint_names.append(j.name)
                    elif j.type == 'prismatic' or j.type == 'revolute':
                        self._joint_lower_limits.append(j.limit.lower)
                        self._joint_upper_limits.append(j.limit.upper)
                        self._joint_vel_limits.append(j.limit.velocity)
                        self._joint_type.append(j.type )
                    else:
                        raise ValueError('Joint type ' + j.type + ' is not supported.')
                    break;

        assert len(self._joint_names) == self._num_jnts, "len(self._joint_names): {}, self._num_jnts: {}".format(len(self._joint_names), self._num_jnts)

        assert len(self._joint_type) == self._num_jnts, "len(self._joint_type): {}, self._num_jnts: {}".format(len(self._joint_type), self._num_jnts)

        # KDL Solvers
        self._fk_p_kdl = PyKDL.ChainFkSolverPos_recursive(self._arm_chain)
        self._fk_v_kdl = PyKDL.ChainFkSolverVel_recursive(self._arm_chain)
        self._ik_v_kdl = PyKDL.ChainIkSolverVel_pinv(self._arm_chain)
        self._ik_p_kdl = PyKDL.ChainIkSolverPos_NR(self._arm_chain,
                                                   self._fk_p_kdl,
                                                   self._ik_v_kdl)
        self._jac_kdl = PyKDL.ChainJntToJacSolver(self._arm_chain)
        self._dyn_kdl = PyKDL.ChainDynParam(self._arm_chain,
                                            PyKDL.Vector.Zero())

    def print_robot_description(self):
        nf_joints = 0
        for j in self._robot.joints:
            if j.type != 'fixed':
                nf_joints += 1
        print ("URDF non-fixed joints: %d;" % nf_joints)
        print ("URDF total joints: %d" % len(self._robot.joints))
        print ("URDF links: %d" % len(self._robot.links))
        print ("KDL joints: %d" % self._kdl_tree.getNrOfJoints())
        print ("KDL segments: %d" % self._kdl_tree.getNrOfSegments())

    def print_kdl_chain(self):
        for idx in xrange(self._arm_chain.getNrOfSegments()):
            print ('* ' + self._arm_chain.getSegment(idx).getName())

    def joints_to_kdl(self, type, values=None):
        kdl_array = PyKDL.JntArray(self._num_jnts)

        cur_type_values = values
        
        for idx in range(self._num_jnts):
            kdl_array[idx] = cur_type_values[idx]
        if type == 'velocities':
            kdl_array = PyKDL.JntArrayVel(kdl_array)
        return kdl_array

    def kdl_to_mat(self, data):
        mat =  np.mat(np.zeros((data.rows(), data.columns())))
        for i in range(data.rows()):
            for j in range(data.columns()):
                mat[i,j] = data[i,j]
        return mat

    def forward_position_kinematics(self,joint_values=None, base_placement = None):
        end_frame = PyKDL.Frame()
        self._fk_p_kdl.JntToCart(self.joints_to_kdl('positions',joint_values),
                                 end_frame)
        pos_wrt_base = end_frame.p
        rot = PyKDL.Rotation(end_frame.M)
        rot_wrt_base = rot.GetQuaternion()

        if base_placement is not None:
            x, y, theta = base_placement
            q_wxyz = T.quaternion_multiply(T.quaternion_about_axis(theta, [0, 0, 1]), [rot_wrt_base[3], rot_wrt_base[0], rot_wrt_base[1], rot_wrt_base[2]])
            rot_wrt_world = np.array([q_wxyz[1], q_wxyz[2], q_wxyz[3], q_wxyz[0]])
            pos_wrt_world = [x + math.cos(theta) * pos_wrt_base[0] - math.sin(theta) * pos_wrt_base[1],
                             y + math.sin(theta) * pos_wrt_base[0] + math.cos(theta) * pos_wrt_base[1],
                                pos_wrt_base[2]] 
            
            return np.array([pos_wrt_world[0], pos_wrt_world[1], pos_wrt_world[2],
                            rot_wrt_world[0], rot_wrt_world[1], rot_wrt_world[2], rot_wrt_world[3]])
        else:
            return np.array([pos_wrt_base[0], pos_wrt_base[1], pos_wrt_base[2],
                            rot_wrt_base[0], rot_wrt_base[1], rot_wrt_base[2], rot_wrt_base[3]])

    def forward_velocity_kinematics(self,joint_velocities=None):
        end_frame = PyKDL.FrameVel()
        self._fk_v_kdl.JntToCart(self.joints_to_kdl('velocities',joint_velocities),
                                 end_frame)
        return end_frame.GetTwist()

    def inverse_kinematics(self, position, orientation, seed):
        ik = PyKDL.ChainIkSolverVel_pinv(self._arm_chain)
        pos = PyKDL.Vector(position[0], position[1], position[2])
        
        rot = PyKDL.Rotation()
        rot = rot.Quaternion(orientation[0], orientation[1],
                                orientation[2], orientation[3])
        
        # Populate seed with current angles if not provided
        seed_array = PyKDL.JntArray(self._num_jnts)
        
        seed_array.resize(len(seed))
        for idx, jnt in enumerate(seed):
            seed_array[idx] = jnt

        goal_pose = PyKDL.Frame(rot, pos)
            
        result_angles = PyKDL.JntArray(self._num_jnts)

        if self._ik_p_kdl.CartToJnt(seed_array, goal_pose, result_angles) >= 0:
            result = np.array(list(result_angles))
            return result
        else:
            return None

    def manipulability(self,joint_values=None):
        jacobian = self.jacobian(joint_values)
        return np.sqrt(np.linalg.det(jacobian * jacobian.T))

    def jacobian(self,joint_values=None):
        jacobian = PyKDL.Jacobian(self._num_jnts)
        self._jac_kdl.JntToJac(self.joints_to_kdl('positions',joint_values), jacobian)
        return self.kdl_to_mat(jacobian)

    def jacobian_transpose(self,joint_values=None):
        return self.jacobian(joint_values).T

    def jacobian_pseudo_inverse(self,joint_values=None):
        return np.linalg.pinv(self.jacobian(joint_values))

    def inertia(self,joint_values=None):
        inertia = PyKDL.JntSpaceInertiaMatrix(self._num_jnts)
        self._dyn_kdl.JntToMass(self.joints_to_kdl('positions',joint_values), inertia)
        return self.kdl_to_mat(inertia)

    def cart_inertia(self,joint_values=None):
        js_inertia = self.inertia(joint_values)
        jacobian = self.jacobian(joint_values)
        return np.linalg.inv(jacobian * np.linalg.inv(js_inertia) * jacobian.T)

