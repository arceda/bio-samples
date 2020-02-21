import os
import subprocess
import kameris_formats
import numpy as np
from Bio import SeqIO
import re

current_dir = os.path.dirname(os.path.abspath(__file__))

def read_fasta(fa_file): 
    data = []

    sequences = SeqIO.parse(fa_file, "fasta")

    for record in sequences:
        data.append([record.id, record.seq.upper()])
        #print(record.id)
        #print(record.seq.upper())     

    return data

def cgr(seq, k, order='ACGT'):
    result = np.zeros(pow(4, k))
    #print(result.shape, result)

    x = pow(2, k-1)
    y = pow(2, k-1)

    print("len seq:", len(seq))
    for i, c  in enumerate(seq):
    #for i, c  in enumerate(range(20)):
        #print(i, c)
        if seq[i] != order[0] and seq[i] != order[1] and seq[i] != order[2] and seq[i] != order[3]:
            continue

        x >>= 1  #x /= 2
        if seq[i] == order[2] or seq[i] == order[3]:
            x += pow(2, k-1)

        y >>= 1  #y /= 2
        if seq[i] == order[0] or seq[i] == order[3]:
            y += pow(2, k-1)
        
        
        if i >= (k - 1):
            #print((pow(2, k))*y + x)
            pos = int( pow(2, k)*y + x )
            #print(pos)
            result[pos] += 1

    return result


data = read_fasta(current_dir + "/hiv1-genomes/A1.fasta")
result = []
for d in data:
    result.append(cgr(d[1], 5)) # the second value save the sequence, the first one is the id





"""
#call implementation in c++


output_file = current_dir + "/results/cgr-k=5"
input_files = current_dir + "/hiv1-genomes/"
k = 5
path_cmd = current_dir + "/../kameris-backend/build/generation_cgr_original " 
cmd = path_cmd + " cgr '" + input_files + "' '" + output_file + "' " + str(k) + " 16 "
#cmd = "/home/vicente/projects/BIOINFORMATICS/bio-samples/viral/kameris-backend/build/generation_cgr cgr '/home/vicente/Descargas/hiv1-genomes/'  '/home/vicente/Descargas/output/cgr-k=9.mm-repr' 6 16"
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
(output, err) = p.communicate()
 
## Wait for date to terminate. Get return returncode ##
p_status = p.wait()

output = str(output).replace('\\n', '\n')
print("Command output : ", output)
print("Command exit status/return code : ", p_status)

cgrs = []
reader = kameris_formats.repr_reader(output_file)
for i in range(reader.count):                
    cgr = reader.read_matrix(i)
    #print("\ni:", i)
    #print("len(cgr)", len(cgr[0]))
    #print("cgr", cgr)

    cgrs.append(cgr[0])
    #cgrs.append(cgr[0]/cgr[0].sum())


#############################################################################
#############################################################################
    # the both comands return the same
    
path_cmd = current_dir + "/../kameris-backend/build/generation_cgr_linux_avx2 "
cmd = path_cmd + " cgr '" + input_files + "' '" + output_file + "' " + str(k) + " 16 "
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
(output, err) = p.communicate()
p_status = p.wait()
output = str(output).replace('\\n', '\n')
print("Command output : ", output)
print("Command exit status/return code : ", p_status)

cgrs2 = []
reader2 = kameris_formats.repr_reader(output_file)
for i in range(reader2.count):                
    cgr = reader2.read_matrix(i)
    #print("\ni:", i)
    #print("len(cgr)", len(cgr[0]))
    #print("cgr", cgr)

    cgrs2.append(cgr[0]/cgr[0].sum())

"""