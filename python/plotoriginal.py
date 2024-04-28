import matplotlib.pyplot as plt
import math

def classify_point(c):
    if c == "A":
        return '$A$',"red"
    elif c == "G":
        return '$G$',"green"
    elif c == "C":
        return '$C$',"blue"
    else: 
        return '$T$',"purple"

with open("files/input.txt","r") as f:
    seq = f.readline().strip()
    dot_bracket = f.readline().strip()

pairs = []
stack = []
for i, char in enumerate(dot_bracket):
    if char == '(':
        stack.append(i)
    elif char == ')':
        pairs.append((stack.pop(), i))

interval = math.radians(360/len(seq))
x = []
y = []

plt.figure(figsize=(12, 12))
plt.axis('off')
for index in range(len(seq)):
    x1 = len(seq)//20 * math.cos((index)*interval)
    y1 = len(seq)//20 * math.sin((index)*interval)
    x.append(x1)
    y.append(y1)
    m,c = classify_point(seq[index])
    plt.scatter(x1, y1, marker=m, c=c, zorder=2)

plt.plot(x, y, 'red', zorder=1, mfc="none", linewidth=1)

for (p1,p2) in pairs:
    h = [x[p1], x[p2]]
    k = [y[p1], y[p2]]
    plt.plot(h, k, 'blue', zorder=1, linestyle="dotted")


plt.savefig("files/originalstructure.png")
# plt.show()
