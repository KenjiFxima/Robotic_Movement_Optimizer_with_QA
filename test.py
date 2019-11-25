import math
<<<<<<< HEAD
import sys
import numpy as np

=======
>>>>>>> ec550e7829356d401baee9f63536ed13a23081f9

def dist_cal(start, end):
    dist = math.sqrt((start[0] - end[0]) ** 2 + (start[1] - end[1]) ** 2)
    return dist
<<<<<<< HEAD

def read_distances(filename):
    task = []
    set = []
    velocity = []
    with open(filename, 'r') as f:
        for line in f:
            # Skip comments
            if line[0] == '#':
                continue
            else:
                l = list(map(int, (line.strip()).split(',')))
                set.append(list(l[0:2]))
                set.append(list(l[2:4]))
                velocity.append(l[4])
    n = int(len(set) / 2)
    for i in range(n):
        task.append(set)
    task = np.array(task)

    return task, velocity

arg = sys.argv[1]
task, velocity = read_distances(arg)
print(task)
n = len(velocity)
start =[0,0]
w0 = [[] for i in range(n)]
w = [[] for i in range(n*2)]
for i in range(n * 2):
    for j in range(n):
        w0[j].append(dist_cal(start, task[0][i]))
    for j in range(n*2):
        w[j].append(dist_cal(task[0][i],task[0][j]))
w0 = np.array(w0)
w = np.array(w)
print(w0)
print(w)
=======
>>>>>>> ec550e7829356d401baee9f63536ed13a23081f9
'''
ans = dist_cal((0,0),(17,479))
ans += dist_cal((287,619),((43,565)))
ans += dist_cal((488,353),(450,142))
ans += dist_cal((688,580),(595,627))
ans += dist_cal((592,52),(630,170))
ans += dist_cal((136,403),(0,0))
<<<<<<< HEAD

=======
'''
>>>>>>> ec550e7829356d401baee9f63536ed13a23081f9
ans = dist_cal((0,0),(136,403))
ans += dist_cal((630,170),((43,565)))
ans += dist_cal((488,353),(287,619))
ans += dist_cal((17,479),(595,627))
ans += dist_cal((592,52),(688,580))
ans += dist_cal((450,142),(0,0))

print(dist_cal((287,619),(17,479)))

<<<<<<< HEAD
print(ans)
'''
=======
print(ans)
>>>>>>> ec550e7829356d401baee9f63536ed13a23081f9
