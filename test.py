import math

def dist_cal(start, end):
    dist = math.sqrt((start[0] - end[0]) ** 2 + (start[1] - end[1]) ** 2)
    return dist
'''
ans = dist_cal((0,0),(17,479))
ans += dist_cal((287,619),((43,565)))
ans += dist_cal((488,353),(450,142))
ans += dist_cal((688,580),(595,627))
ans += dist_cal((592,52),(630,170))
ans += dist_cal((136,403),(0,0))
'''
ans = dist_cal((0,0),(136,403))
ans += dist_cal((630,170),((43,565)))
ans += dist_cal((488,353),(287,619))
ans += dist_cal((17,479),(595,627))
ans += dist_cal((592,52),(688,580))
ans += dist_cal((450,142),(0,0))

print(dist_cal((287,619),(17,479)))

print(ans)