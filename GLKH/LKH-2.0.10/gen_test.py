import numpy as np
import os
from subprocess import check_call
import time

# a graph with nxn nodes 
class Graph:
    def __init__(self, m):
        self.n = m * m
        self.m = m
        self.cost = np.ones((self.n+1, self.n+1)) * 536870911
        for i in range(self.n+1):
            self.cost[i][i] = 0
            self.cost[i][self.n] = 0
            self.cost[self.n][i] = 0

    def generate(self):
        for x in range(self.m):
            for y in range(self.m):

                i = x * self.m + y

                if x > 0:
                    j = (x - 1) * self.m + y
                    self.cost[i][j] = np.random.randint(1, 100)
                
                if x < self.m - 1:
                    j = (x + 1) * self.m + y
                    self.cost[i][j] = np.random.randint(1, 100)
                
                if y > 0:
                    j = x * self.m + y - 1
                    self.cost[i][j] = np.random.randint(1, 100)
                
                if y < self.m - 1:
                    j = x * self.m + y + 1
                    self.cost[i][j] = np.random.randint(1, 100)

    def print(self):
        # print as matrix
        for i in range(self.n):
            for j in range(self.n):
                print(self.cost[i][j], end=' ')
            print()

    def write_full_matrix(self):
        with open("test_dense.tsp", "w") as f:
            f.write("NAME: test_dense\n")
            f.write("TYPE: ATSP\n")
            f.write("DIMENSION: {}\n".format(self.n+1))
            f.write("EDGE_WEIGHT_TYPE: EXPLICIT\n")
            f.write("EDGE_WEIGHT_FORMAT: FULL_MATRIX\n")
            f.write("EDGE_WEIGHT_SECTION\n")
            for i in range(self.n+1):
                for j in range(self.n+1):
                    f.write("{} ".format(int(self.cost[i][j])))
                f.write("\n")
            f.write("EOF\n")

    def solve_full_matrix(self):
        self.write_full_matrix()
        # run LKH
        start = time.time()
        with open("test_dense.log", "w") as f:
            check_call(["./LKH", "test_dense.par"], stdout=f)
        print(time.time() - start)

    def write_sparse_matrix(self):
        with open("test_sparse.tsp", "w") as f:
            f.write("NAME: test_sparse\n")
            f.write("TYPE: ATSP\n")
            f.write("DIMENSION: {}\n".format(self.n+1))
            f.write("EDGE_DATA_FORMAT: EDGE_LIST\n")
            f.write("EDGE_DATA_SECTION:\n")
            for i in range(self.n+1):
                for j in range(self.n+1):
                    if self.cost[i][j] != 536870911 and i != j:
                        f.write("{} {} {}\n".format(i+1, j+1, int(self.cost[i][j])))
            f.write("EOF\n")
    
    def solve_sparse_matrix(self):
        self.write_sparse_matrix()
        # run LKH
        start = time.time()
        with open("test_sparse.log", "w") as f:
            check_call(["./LKH", "test_sparse.par"], stdout=f)
        print(time.time() - start)
        # os.system("./LKH test_sparse.par")
    
if __name__ == '__main__':
    g = Graph(30)
    g.generate()
    g.solve_full_matrix()
    g.solve_sparse_matrix()
    # g.print()

