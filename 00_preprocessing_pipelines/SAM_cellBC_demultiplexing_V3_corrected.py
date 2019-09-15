import os
import sys

big_sam_file = open(sys.argv[1], "r")
bc_file = open(sys.argv[2], "r")
protocol = sys.argv[3]
sample = sys.argv[4]
output_path = sys.argv[5]


bc_list = []
dict_demux = {}
bc_line = bc_file.readline()
while bc_line:
    bc = bc_line.split()[0]
    if bc != 'XC':
        bc_list.append(bc)
        dict_demux[bc] = []
    bc_line = bc_file.readline()
bc_file.close()


header_list = []
read_line = big_sam_file.readline()
while read_line:
    if read_line.startswith("@"):
        header_list.append(read_line)
    else:
        line_words = read_line.split()
        BC_tag = [i for i in line_words if i.startswith('BC:Z:')][0]
        cellBC = BC_tag.split(":")
        if len(cellBC) != 3:
            print(cellBC)
        cellBC = cellBC[2]
        if cellBC in dict_demux.keys():
            line_partial = '\t'.join(line_words[0:4])
            line_partial = line_partial + "\n"
            dict_demux[cellBC].append(line_partial)
    read_line = big_sam_file.readline()
big_sam_file.close()
            

for cBC in bc_list:
    filenames = output_path + "/" + protocol + "." + "mixed." + sample + "." + cBC + ".sam"
    fout = open(filenames, 'w')
    for hline in header_list:
        fout.write(hline)
    for mapped_line in dict_demux[cBC]:
        fout.write(mapped_line)
    fout.close()
        
        