# coding: utf-8
def answer(s):
    odd_list=[2,3,5,7,11,13,17,19]
    odd_list.reverse()
    s_iter=s
    sum=1

    if(s_iter.count(s_iter[0:1]) == len(s_iter)):
        sum=len(s_iter)
    else:
        while (len(odd_list)>0):
            for odd in odd_list:
                if(len(s_iter)%odd > 0):
                    odd_list.pop(0)
                else:
                    cake_size=len(s_iter)//odd
                    sub_a=s_iter[0:cake_size]
                    sub_b=s_iter[cake_size:cake_size+cake_size]
                    if(sub_a == sub_b):
                        sum = sum * odd
                        s_iter = sub_a
                    else:
                        odd_list.pop(0)
    return sum

import fileinput

for line in fileinput.input():
    print(answer(line.strip()))

# coding: utf-8
def answer(s):
    odd_list=[19, 17, 13, 11, 7, 5, 3, 2]
    s_iter=s
    sum=1
    while (len(odd_list)>0):
        for odd in odd_list:
            slen=len(s_iter)
            if(slen%odd != 0):
                odd_list.pop(0)
            else:
                cake_size=slen//odd
                Flag=True
                sub_a=s_iter[0:cake_size]
                for check in range(cake_size, slen, cake_size):
                    sub_b=s_iter[check:check+cake_size]
                    if(sub_a != sub_b):
                        Flag=False
                        break
                if(Flag):
                    sum = sum * odd
                    s_iter = sub_a
                else:
                    odd_list.pop(0)
    return sum
