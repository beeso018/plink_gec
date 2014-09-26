import subprocess, shlex, fileinput, sys, os, glob

os.makedirs(".output", exist_ok=True)

filename=input("\nFile name?\n")

plink1="./plink --file "
plink2=" --silent --dog --nonfounders --allow-no-sex --snps "
plink3=" --recode --out "
snp1_22="BICF2G630707759-G1212f41S99"
snp23_39="BICF2P653617-TIGRP2P134917_rs8450925"
print("\nPLINKing...\n")

plink=plink1+filename+plink2+snp1_22+plink3+".output/chr1_22"
plink=shlex.split(plink)
p=subprocess.Popen(plink).wait()
plink=plink1+filename+plink2+snp23_39+plink3+".output/chr23_39"
plink=shlex.split(plink)
p=subprocess.Popen(plink).wait()

os.chdir(".output")

for line in fileinput.input("chr23_39.map", inplace=1):
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
    elif line.startswith("33"):
        line=line.replace("33", "11",1)
        print(line, end='')
    elif line.startswith("34"):
        line=line.replace("34", "12",1)
        print(line, end='')
    elif line.startswith("35"):
        line=line.replace("35", "13",1)
        print(line, end='')
    elif line.startswith("36"):
        line=line.replace("36", "14",1)
        print(line, end='')
    elif line.startswith("37"):
        line=line.replace("37", "15",1)
        print(line, end='')
    elif line.startswith("38"):
        line=line.replace("38", "16",1)
        print(line, end='')
    elif line.startswith("X"):
        line=line.replace("X", "17",1)
        print(line, end='')

os.chdir("..")

print("Doing some GEC stuff...\n")

gec1="java -jar -Xmx1g gec.jar --no-web --effect-number --linkage-file "
gec2=" --genome --out "

gec=gec1+".output/chr1_22"+gec2+".output/chr1_22"
gec=shlex.split(gec)
p1=subprocess.Popen(gec, stdout=subprocess.PIPE)
chr1_22_out=p1.communicate()[0]
chr1_22_gec_out=open('.output/gec_out', 'w')
chr1_22_gec_out.write(chr1_22_out.decode('utf-8'))

gec=gec1+".output/chr23_39"+gec2+".output/chr23_39"
gec=shlex.split(gec)
p2=subprocess.Popen(gec, stdout=subprocess.PIPE)
chr23_39_out=p2.communicate()[0]
chr23_39_gec_out=open('.output/gec_out2', 'w')
chr23_39_gec_out.write(chr23_39_out.decode('utf-8'))
chr23_39_gec_out.close()

os.chdir(".output")

with open("chr1_22.sum", "r") as chr1_22:
    for i in range(1):
        line=chr1_22.readline().strip()
    for line in chr1_22:
        out1=line.strip().split('\t')
with open("chr23_39.sum", "r") as chr23_39:
    for i in range(1):
        line=chr23_39.readline().strip()
    for line in chr23_39:
        out2=line.strip().split('\t')

obs=eval(out1[0])+eval(out2[0])
eff=eval(out1[1])+eval(out2[1])
suggestive=0.1/eff
significant=0.05/eff
hsignificant=0.001/eff

for line in fileinput.input("chr23_39.block.txt", inplace=1):
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
        elif out3[0]=="11":
            print("33", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="12":
            print("34", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="13":
            print("35", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="14":
            print("36", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="15":
            print("37", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="16":
            print("38", out3[1], out3[2], out3[3], out3[4], sep='\t', end='\n')
        elif out3[0]=="17":
            print("39", out3[1], out3[2], out3[3], out3[4], sep='\t', end='')

blocks=open("gec_blocks.txt", 'w')
for line in fileinput.input("chr1_22.block.txt"):
    blocks.write(line)
blocks2=open("gec_blocks.txt", 'a')
for line in fileinput.input("chr23_39.block.txt"):
    blocks2.write(line)

os.rename("gec_blocks.txt", "../gec_blocks.txt")
for i in glob.glob(u'chr*'):
    os.unlink(i)

gec_out3=open("gec_out3.txt", 'a')
for line in fileinput.input("gec_out"):
    if line.startswith("The number"):
        line=line.replace("chr1_22.map", "")
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
        line=line.replace("chr23_39.map", "")
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
        elif "chromosome 11 " in line:
            line=line.replace("chromosome 11 ", "chromosome 33 ")
            gec_out3.write(line)
        elif "chromosome 12 " in line:
            line=line.replace("chromosome 12 ", "chromosome 34 ")
            gec_out3.write(line)
        elif "chromosome 13 " in line:
            line=line.replace("chromosome 13 ", "chromosome 35 ")
            gec_out3.write(line)
        elif "chromosome 14 " in line:
            line=line.replace("chromosome 14 ", "chromosome 36 ")
            gec_out3.write(line)
        elif "chromosome 15 " in line:
            line=line.replace("chromosome 15 ", "chromosome 37 ")
            gec_out3.write(line)
        elif "chromosome 16 " in line:
            line=line.replace("chromosome 16 ", "chromosome 38 ")
            gec_out3.write(line)
        elif "chromosome 17 " in line:
            line=line.replace("chromosome 17 ", "chromosome 39 ")
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
        elif "chromosome 11" in line:
            line=line.replace("chromosome 11", "chromosome 33")
            gec_out3.write(line)
        elif "chromosome 12" in line:
            line=line.replace("chromosome 12", "chromosome 34")
            gec_out3.write(line)
        elif "chromosome 13" in line:
            line=line.replace("chromosome 13", "chromosome 35")
            gec_out3.write(line)
        elif "chromosome 14" in line:
            line=line.replace("chromosome 14", "chromosome 36")
            gec_out3.write(line)
        elif "chromosome 15" in line:
            line=line.replace("chromosome 15", "chromosome 37")
            gec_out3.write(line)
        elif "chromosome 16" in line:
            line=line.replace("chromosome 16", "chromosome 38")
            gec_out3.write(line)
        elif "chromosome 17" in line:
            line=line.replace("chromosome 17", "chromosome 39")
            gec_out3.write(line)

os.rename("gec_out3.txt", "../gec_out3.txt")
os.chdir("..")
os.rename("gec_out3.txt", filename+"_gec_out.txt")
os.rename("gec_blocks.txt", filename+"_gec_blocks.txt")

print("Observed markers: ", obs)
print("Effective markers: ", eff)
print("Suggestive p-value: <", format(suggestive, '.3e'))
print("Significant p-value: <", format(significant, '.3e'))
print("Highly significant p-value: <", format(hsignificant, '.3e'), "\n")
print("Block-wise summary can be found in ", filename, ".block.txt", sep='')
print("Effective marker summary by chromosome can be found in ", filename, "_gec_out.txt\n", sep='')
