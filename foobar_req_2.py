# coding: utf-8
def answer(s):
    sum=0
    for check in range(0, len(s)):
        if(s[check:check+1] != '>'):
            continue
        else:
            subs=s[check:len(s)]
            sum=subs.count('<')+sum
    sum=sum+sum
    return sum
