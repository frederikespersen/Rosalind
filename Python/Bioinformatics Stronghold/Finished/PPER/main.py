input = '85 8'
input = input.split()

n = int(input[0])
k = int(input[1])

from math import factorial as f
P = f(n)/f(n-k)

print(int(P) % 1000000)