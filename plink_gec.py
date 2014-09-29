import subprocess, shlex, fileinput, sys, os, glob, shutil

# make hidden directory
os.makedirs(".output", exist_ok=True)

# ask user for filename from map/ped files
filename=input("\nFile name?\n")

# ask user which SNP chip to set SNP boundaries from
snpchip=input("\nWhich SNP chip? (50K or 70K)\n")

# set SNP boundaries in PLINK commands
plink1="./plink --file "
plink2=" --silent --horse --nonfounders --allow-no-sex --snps "
plink3=" --recode --out "
if snpchip=="50K":
    snp1_22="chr1.72460-chr22.49944830"
    snp23_32="chr23.41804-chrX.124108865"
elif snpchip=="70K":
    snp1_22="chr1.15656-chr22.49944830"
    snp23_32="chr23.4510-chrX.124108865"
print("\nPLINKing...\n")

# string PLINK commands together & run as subprocess for GEC input
plink=plink1+filename+plink2+snp1_22+plink3+".output/chr1_22"
plink=shlex.split(plink)
p=subprocess.Popen(plink).wait()
plink=plink1+filename+plink2+snp23_32+plink3+".output/chr23_32"
plink=shlex.split(plink)
p=subprocess.Popen(plink).wait()

os.chdir(".output")

# rename chromosomes greater than 22 in map file for GEC recognition
for line in fileinput.input("chr23_32.map", inplace=1):
    if line.startswith("23"):
        line=line.replace("23","1",1)
        print(line, end='')
    elif line.startswith("24"):
        line=line.replace("24", "2",1)
        print(line, end='')
    elif line.startswith("25"):
        line=line.replace("25", "3",1)
        print(line, end='')
    elif line.startswith("26"):
        line=line.replace("26", "4",1)
        print(line, end='')
    elif line.startswith("27"):
        line=line.replace("27", "5",1)
        print(line, end='')
    elif line.startswith("28"):
        line=line.replace("28", "6",1)
        print(line, end='')
    elif line.startswith("29"):
        line=line.replace("29", "7",1)
        print(line, end='')
    elif line.startswith("30"):
        line=line.replace("30", "8",1)
        print(line, end='')
    elif line.startswith("31"):
        line=line.replace("31", "9",1)
        print(line, end='')
    elif line.startswith("32"):
        line=line.replace("32", "10",1)
        print(line, end='')

os.chdir("..")

print("Doing some GEC stuff...\n")

# string GEC commands together and run as subprocess
gec1="java -jar -Xmx1g gec.jar --no-web --effect-number --linkage-file "
gec2=" --genome --out "
gec=gec1+".output/chr1_22"+gec2+".output/chr1_22"
gec=shlex.split(gec)
p1=subprocess.Popen(gec, stdout=subprocess.PIPE)
chr1_22_out=p1.communicate()[0]
chr1_22_gec_out=open('.output/gec_out', 'w')
chr1_22_gec_out.write(chr1_22_out.decode('utf-8'))
gec=gec1+".output/chr23_32"+gec2+".output/chr23_32"
gec=shlex.split(gec)
p2=subprocess.Popen(gec, stdout=subprocess.PIPE)
chr23_32_out=p2.communicate()[0]
chr23_32_gec_out=open('.output/gec_out2', 'w')
chr23_32_gec_out.write(chr23_32_out.decode('utf-8'))
chr23_32_gec_out.close()

os.chdir(".output")

# get GEC summary stats
with open("chr1_22.sum", "r") as chr1_22:
    for i in range(1):
        line=chr1_22.readline().strip()
    for line in chr1_22:
        out1=line.strip().split('\t')
with open("chr23_32.sum", "r") as chr23_32:
    for i in range(1):
        line=chr23_32.readline().strip()
    for line in chr23_32:
        out2=line.strip().split('\t')

# calculate total observed and effective markers
obs=eval(out1[0])+eval(out2[0])
eff=eval(out1[1])+eval(out2[1])

# calculate corresponding p-value thresholds
suggestive=0.1/eff
significant=0.05/eff
hsignificant=0.001/eff

# rename chromosomes >22 again
for line in fileinput.input("chr23_32.block.txt", inplace=1):
        out3=line.strip().split('\t')
        if out3[0]=="1":
            print("23", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="2":
            print("24", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="3":
            print("25", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="4":
            print("26", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="5":
            print("27", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="6":
            print("28", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="7":
            print("29", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="8":
            print("30", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="9":
            print("31", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="10":
            print("32", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')

# combine block files
blocks=open("gec_blocks.txt", 'w')
for line in fileinput.input("chr1_22.block.txt"):
    blocks.write(line)
blocks2=open("gec_blocks.txt", 'a')
for line in fileinput.input("chr23_32.block.txt"):
    blocks2.write(line)

# move combined block file to working directory
os.rename("gec_blocks.txt", "../gec_blocks.txt")

# format GEC output and remove extraneous information
gec_out3=open("gec_out3.txt", 'a')
for line in fileinput.input("gec_out"):
    if line.startswith("The number"):
        line=line.replace("chr1_22.map ", "")
        if "X" in line or "Y" in line or "MT" in line:
            line=""
        else:
            gec_out3.write(line)
    elif "MAF" in line:
        gec_out3.write(line)
    elif line.startswith("The estimated"):
        gec_out3.write(line)
for line in fileinput.input("gec_out2"):
    if line.startswith("The number"):
        line=line.replace("chr23_32.map ", "")
        if "chromosome 1 " in line:
            line=line.replace("chromosome 1 ", "chromosome 23 ")
            gec_out3.write(line)
        elif "chromosome 2 " in line:
            line=line.replace("chromosome 2 ", "chromosome 24 ")
            gec_out3.write(line)
        elif "chromosome 3 " in line:
            line=line.replace("chromosome 3 ", "chromosome 25 ")
            gec_out3.write(line)
        elif "chromosome 4 " in line:
            line=line.replace("chromosome 4 ", "chromosome 26 ")
            gec_out3.write(line)
        elif "chromosome 5 " in line:
            line=line.replace("chromosome 5 ", "chromosome 27 ")
            gec_out3.write(line)
        elif "chromosome 6 " in line:
            line=line.replace("chromosome 6 ", "chromosome 28 ")
            gec_out3.write(line)
        elif "chromosome 7 " in line:
            line=line.replace("chromosome 7 ", "chromosome 29 ")
            gec_out3.write(line)
        elif "chromosome 8 " in line:
            line=line.replace("chromosome 8 ", "chromosome 30 ")
            gec_out3.write(line)
        elif "chromosome 9 " in line:
            line=line.replace("chromosome 9 ", "chromosome 31 ")
            gec_out3.write(line)
        elif "chromosome 10 " in line:
            line=line.replace("chromosome 10 ", "chromosome 32 ")
            gec_out3.write(line)
    elif "MAF" in line:
        gec_out3.write(line)
    elif line.startswith("The estimated"):
        if "chromosome 1\n" in line:
            line=line.replace("chromosome 1", "chromosome 23")
            gec_out3.write(line)
        elif "chromosome 2" in line:
            line=line.replace("chromosome 2", "chromosome 24")
            gec_out3.write(line)
        elif "chromosome 3" in line:
            line=line.replace("chromosome 3", "chromosome 25")
            gec_out3.write(line)
        elif "chromosome 4" in line:
            line=line.replace("chromosome 4", "chromosome 26")
            gec_out3.write(line)
        elif "chromosome 5" in line:
            line=line.replace("chromosome 5", "chromosome 27")
            gec_out3.write(line)
        elif "chromosome 6" in line:
            line=line.replace("chromosome 6", "chromosome 28")
            gec_out3.write(line)
        elif "chromosome 7" in line:
            line=line.replace("chromosome 7", "chromosome 29")
            gec_out3.write(line)
        elif "chromosome 8" in line:
            line=line.replace("chromosome 8", "chromosome 30")
            gec_out3.write(line)
        elif "chromosome 9" in line:
            line=line.replace("chromosome 9", "chromosome 31")
            gec_out3.write(line)
        elif "chromosome 10" in line:
            line=line.replace("chromosome 10", "chromosome 32")
            gec_out3.write(line)

# finish formatting GEC output
output=[]
output.append(['CHR', 'SNPs', 'MAF>0.05', 'EFF', 'RATIO'])
with open("gec_out3.txt", 'r') as gec_out:
    lines=[line for line in gec_out]
    for i in range(0,len(lines),3):
        line1=lines[i].strip().split(' ')
        line2=lines[i+1].strip().split(' ')
        line3=lines[i+2].strip().split(' ')
        output.append([line1[6], line1[12].replace('.', ''), line2[0],
                       format(float(line3[8]), '.2f'),
                       format(float(line3[14]), '.2f')])
with open("gec_output", 'w') as gec_output:
    for item in output:
        gec_output.write("{}\n".format(item).replace(',','').replace("'","").replace('[','').replace(']','').replace(' ','\t'))

# rename output to include original filename
os.rename("gec_output", "../"+filename+"_gec_out.txt")
os.chdir("..")
os.rename("gec_blocks.txt", filename+"_gec_blocks.txt")
shutil.rmtree(".output")

# print calculated values
print("Observed markers: ", obs, sep='\t')
print("Effective markers: ", eff, sep='\t')
print("Suggestive p-value: ", format(suggestive, '.3e'), sep='\t')
print("Significant p-value: ", format(significant, '.3e'), sep='\t')
print("Highly significant p-value: <", format(hsignificant, '.3e'), "\n", sep='\t')
print("Block-wise summary can be found in ", filename, ".block.txt", sep='')
print("Effective marker summary by chromosome can be found in ", filename, "_gec_out.txt\n", sep='')
