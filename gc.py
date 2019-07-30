# Calculates GC content

gene = open("BCRA1_BAP1.txt", "r")
print(gene)

# set the values to 0
g=0
a=0
c=0
t=0

# skip de the first line of header
gene.readline()

for line in gene:
    line = line.lower() # to lower case
    for char in line:
        if char == 'g':
            g += 1
        if char == 'a':
            a += 1
        if char == 'c':
            c += 1
        if char == 't':
            t += 1

print "number of g's: " + str(g)
print "number of c's: " + str(c)
print "number of a's: " + str(a)
print "number of t's: " + str(t)

# converting to float
gc = (g + c + 0.)/(a + t + c + g + 0.)
print ("gc content ") + str(gc)
