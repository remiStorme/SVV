from math import sqrt

class VonMises:

    def __init__(self,tensor,yieldCrit,elastLimit,k):
        self.tensor = tensor
        self.yieldCrit = yieldCrit
        self.elastLimit = elastLimit
        self.k = k
        self.tv = 3*self.k**3

    def StressCriterion(self):
        if self.tv**2 >= self.elastLimit:
            print(f"Material will yield")
        x = self.evaluation()
        if self.tv >= x:
            print(f"Material will yield")

    def __str__(self):
        return f"Tv = {self.tv}\nk = {self.k}\nElastic Limit = {self.elastLimit}"

    def evaluation(self):
        return sqrt((1/2)*((self.tensor[0,0]-self.tensor[1,1])**2 +
                       (self.tensor[1,1]-self.tensor[2,2])**2 +
                       (self.tensor[2,2]-self.tensor[0,0])**2 +
                     6*(self.tensor[0,1]**2 + self.tensor[1,2]**2 + self.tensor[0,2]**2)))


