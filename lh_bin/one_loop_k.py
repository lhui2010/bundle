
for i in range(-4,4+1):
    O = " "*abs(i) + "*" + " "*2*(4-abs(i)) + "*" +  " "*abs(i)
    K = "*" + " "*abs(i) + "*"
    print(O + "  " +  K)
