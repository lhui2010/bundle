# coding: utf-8

def answer(pegs):
    # your code here
    strand = -1

    #descending order
    indices = list(range(len(pegs)-1))
    indices.reverse()

    #use formular get x
    sum = pegs[-1]
    for index in indices:
        sum = sum + pegs[index] * 2 * strand
        strand = strand * -1
    sum = sum + pegs[0] * strand
    numerator = sum * 2 * strand
    
    #use formular get y
    if(strand >0):
        denominator = 3
        #if modulus by 3
        if(numerator % denominator == 0):
            numerator = numerator // denominator
            denominator = 1
    else:
        denominator = 1
        
    #and radius
    radius = numerator / denominator

    #make sure each peg's radius do not exceed the neighboring distance
    #Descending order
    BadFlag = False
    this_radius  = radius/2

    for index in indices:
        this_radius = pegs[index + 1] - pegs[index] - this_radius
        threshold = pegs[index + 1] - pegs[index]
        if(this_radius >= threshold or this_radius <= 1 ):
            BadFlag=True

    #Final judgement
    if(radius < 1 or BadFlag ):
        numerator = -1
        denominator = -1

    return(numerator, denominator)




print("4,20")
print(answer([4,20]))
print("4,30, 40, 50")
print(answer([4,30, 40, 50]))
print ("4,40")
print(answer([4,40]))
print(answer([4,40, 80, 120, 160, 200, 240, 250]))
