import matplotlib.pyplot as plt
from mpl_interactions import ioff, panhandler, zoom_factory
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

f = open("files/output.txt","r")
pairs = []
seq = f.readline()
for line in f:
    pairs.append((int(line.split(" ")[0]),int(line.split(" ")[1])))

interval = math.radians(360/len(seq))
x = []
y = []
with plt.ioff():
    figure,axis = plt.subplots()
plt.axis('off')
for index in range(len(seq)):
    x1 = len(seq)//20* math.cos((index)*interval)
    y1 = len(seq)//20 *math.sin((index)*interval)
    x.append(x1)
    y.append(y1)
    m,c = classify_point(seq[index])
    plt.scatter(x1,y1,marker=m,c=c,zorder =2)

plt.plot(x,y,'red',zorder=1,mfc="none",linewidth = 1)
for (p1,p2) in pairs:
    h = [x[p1],x[p2]]
    k = [y[p1],y[p2]]
    plt.plot(h,k,'blue',zorder=1,linestyle="dotted")
figure.set_size_inches((12,12))
disconnect_zoom = zoom_factory(axis)
pan_handler = panhandler(figure)
plt.savefig("files/structure.png")
plt.show()