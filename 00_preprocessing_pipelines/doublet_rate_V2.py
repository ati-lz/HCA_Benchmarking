import os
import sys

allcells_hsap = []
allcells_mmus = []
allcells_hsap_mmus = []

samfiles_path = sys.argv[1] #path to cell demultiplexe sam files
#doublet_table_path = sys.argv[2]
#doublet_table_reads_path = sys.argv[3]

fout = open(sys.argv[2], 'w') # the hsap, mmus percentages
fout2 = open(sys.argv[3],'w') # the hsam, mmus percentages with number
fout3 = open(sys.argv[4],'w') # only human cells
fout4 = open(sys.argv[5],'w') # only mouse cells
fout5 = open(sys.argv[6],'w') # only dog cells
fout6 = open(sys.argv[7],'w') # only human cells bam files ID

for root,dirs,files in os.walk(samfiles_path):
    species = ''
    for file in files:
        path = os.path.join(root,file)
        fin = open(path,'r')
        cell_id = file.split(".")[3] #change this accordingly
        cell_name = ".".join(file.split(".")[0:4]) #change this accordingly
        cell_name_bam = cell_name + ".sam.sorted.bam"
        

        dict_mmus = {}
        dict_hsap = {}
        dict_can = {}
        total_reads = 0
        uniq_hsap = 0
        uniq_mmus = 0
        uniq_can = 0
        for line in fin:
            if not line.startswith("@"):
                total_reads = total_reads +1
                words = line.split()
                read_id = words[0]
                chrom = words[2]
                if "hsap" in chrom:
                    if read_id not in dict_hsap:
                        dict_hsap[read_id] = 1
                    else:
                        dict_hsap[read_id] += 1
                elif "mmus" in chrom:
                    if read_id not in dict_mmus:
                        dict_mmus[read_id] = 1
                    else:
                        dict_mmus[read_id] += 1
                elif "can" in chrom:
                    if read_id not in dict_can:
                        dict_can[read_id] = 1
                    else:
                        dict_can[read_id] += 1
                        
        for read in dict_hsap:
            if read not in dict_mmus and read not in dict_can:
                uniq_hsap += dict_hsap[read]
                
        for read in dict_mmus:
            if read not in dict_hsap and read not in dict_can:
                uniq_mmus += dict_mmus[read]  

        for read in dict_can:
            if read not in dict_hsap and read not in dict_mmus:
                uniq_can += dict_can[read]   

        total = float(uniq_hsap + uniq_mmus + uniq_can)
        if uniq_hsap != 0:
            hsap_per = round((uniq_hsap / total),2) * 100
        else:
            hsap_per = 0
        if uniq_mmus != 0:
            mmus_per = round((uniq_mmus / total),2) * 100
        else:
            mmus_per = 0
        if uniq_can != 0:
            can_per = round((uniq_can / total),2) * 100
        else:
            can_per = 0
        #allcells_hsap.append(hsap_per)
        #allcells_mmus.append(mmus_per)
        #allcells_hsap_mmus.append((hsap_per, mmus_per))
        if hsap_per >= 70:
            fout3.write(cell_name + '\n')
            fout6.write(cell_name_bam + '\n')
            species = "Human"
        elif mmus_per >= 70:
            fout4.write(cell_name + '\n')
            species = "Mouse"
        elif can_per >= 70:
            fout5.write(cell_name + '\n')
            species = "Dog"
        else:
            species = "Doublet"
        fstr_line1 = "cell_id\tspecies\thsap_percentage\tmmus_percentage\n"
        fstr_line2 = "cell_id\tspecies\thsap_uniq_reads\tmmus_uniq_reads\thsap_percentage\tmmus_percentage\n"
        fstr = cell_name + "\t" + species + "\t" + str(hsap_per) + "\t" + str(mmus_per) + "\t" + str(can_per)+ "\n"
        fstr2 = cell_name + "\t" + species + "\t"+ str(uniq_hsap) + "\t" + str(uniq_mmus) + "\t" + str(uniq_can) + "\t" + str(hsap_per) + "\t" + str(mmus_per) + "\t" + str(can_per) + "\n"
        fout2.write(fstr2)
        fout.write(fstr)

fout.close()
fout2.close()
fout3.close()
fout4.close()
fout5.close()
fout6.close()


'''
allcells_hsap = [10,25,90,100,40,98,15]
allcells_mmus = [90,75,10,0,60,2,85]

labels = ["NA"]*len(allcells_hsap)
colors = ["red","green","blue"]
for i in range(0,len(allcells_hsap)):
    if allcells_hsap[i] >= 90:
        labels[i] = "Human"
        #labels[i] = 1
    elif allcells_hsap[i] <= 10:
        labels[i] = "Mouse"
        #labels[i] = 2
    else:
        labels[i] = "Mixed"
        #labels[i] = 3

df = pd.DataFrame(dict(Human=allcells_hsap, Mouse=allcells_mmus, color=labels))
fig, ax = plt.subplots()
colors = {'Human':'red', 'Mouse':'green', 'Mixed':'gray'}
ax.scatter(df['Human'], df['Mouse'], c=df['color'].apply(lambda x: colors[x]))

plt.show()


#matplotlib.pyplot.scatter(allcells_hsap,allcells_mmus, c=labels,cmap=matplotlib.colors.ListedColormap(colors))
#matplotlib.pyplot.show()
'''