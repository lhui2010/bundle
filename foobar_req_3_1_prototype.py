# coding: utf-8
def answers(m):
    denominator = sum(m[0]) * sum(m[1]) - m[0][1] * m[1][0]
    norminator = []
    for i in range(2, len(m)):
        norminator.append(m[0][i] * sum(m[1])+ m[1][i]*m[0][1])
    norminator.append(denominator)
    
    return norminator
