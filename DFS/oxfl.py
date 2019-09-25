import numpy as np

crow = np.array([9,11,5,8,2,13,14,1,4,3,7,10,6,12])
food = np.array([4,10,7,2,13,8,11,5,14,12,6,1,3,9])
sol = np.zeros(14)

#C1 y C2 van de 0 a N-1
C1 = 2
C2 = 11

sol[C1:C2] = crow[C1:C2]

temp = np.setdiff1d(food,crow[C1:C2]) 
sol[0:C1] = temp[0:C1]
#sol[C2+1:crow.shape[0]] = 


print(sol)
print(temp) 