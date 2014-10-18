using brml

# problem 3.17

A,B,C=1,2,3

pA = PotArray(A, [1//2, 1//2])

pBgA = PotArray([B A], [1//2 1//5; 1//2 1//5; 0 3//5])

pCgB = PotArray([C B], [1//2 1//4 3//8; 1//2 3//4 5//8])

pABC = pCgB * pBgA * pA
pAC = sum(pABC, B)

pA=sum(pAC,C)
pC=sum(pAC,A)

println("pAC-pA*pC=")
println(pAC-pA*pC)



# problem 3.20


w,h,inc=1,2,3
low,high=1,2,3

pInc = PotArray(inc, [0.8 0.2])
pWgI = PotArray([w inc], [0.7 0.2; 0.3 0.1; 0 0.4; 0 0.3])
pHgI = PotArray([h inc], [0.2 0; 0.8 0; 0 0.3; 0 0.7])

pWHI = pWgI * pHgI * pInc
pWH = sum(pWHI, inc)

pW = sum(pWH, h)
pH = sum(pWH, w)

println("pWH - pW*pH= should not be zero")
pWH - pW*pH

pWHgI = pWHI / pInc

println("pWHgI - pWgI*pHgI= should be zero")
pWHgI - pWgI*pHgI


# problem 3.22


c12,c13,c23,c32,c31,c21 = 1,2,3,4,5,6
c122, c232 = 1,2
pC12 = PotArray(c12, [0.9 0.1])
pC13 = PotArray(c13, [0.9 0.1])
pC23 = PotArray(c23, [0.9 0.1])
pC32 = PotArray(c32, [0.9 0.1])
pC31 = PotArray(c31, [0.9 0.1])
pC21 = PotArray(c21, [0.9 0.1])

pC232gC23C21C13 = PotArray([c232 c23 c21 c12])


c122, c232, c121323, c232113 = 1,2,3,4

indexC12 = PotArray(c121323, [0 0 0 0 1 1 1 1])

pC121323 = PotArray(c121323, [0.9^3 0.9^2*0.1 0.9^2*0.1 0.9*0.1^2 0.9^2*0.1 0.9*0.1^2 0.9*0.1^2 0.1^3])
pC232113 = PotArray(c232113, [0.9^3 0.9^2*0.1 0.9^2*0.1 0.9*0.1^2 0.9^2*0.1 0.9*0.1^2 0.9*0.1^2 0.1^3])

pC122gC121323 = PotArray([c122 c121323], [1 1 1 0 0 0 0 0; 0 0 0 1 1 1 1 1])
pC232gC232113 = PotArray([c232 c232113], [1 1 1 0 0 0 0 0; 0 0 0 1 1 1 1 1])

pC122 = PotArray(c122, [0; 1])
pC232 = PotArray(c232, [1; 0])

lamC122 = sum(pC122 * pC122gC121323, c122)
lamC232 = sum(pC232 * pC232gC232113, c232)

p1C121323 = lamC122 * pC121323
p1C232113 = lamC232 * pC232113

p1C12 = p1C121323



























