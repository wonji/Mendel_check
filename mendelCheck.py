# mendelCheck.py

print "\n"
print "####################################################################################"
print "####################################################################################"
print "################ mendelCheck.py - dnjswlzz11@gmail.com (2017.02.27) ################"
print "####################################################################################"
print "####################################################################################"
print "\n"


print "usage: python mendelCheck.py plink_input output variant_rate individual_rate missing_code\n"

import sys
import datetime

s1 = datetime.datetime.now()
print "\nThe anlaysis is started at "+str(s1)+".\n"

print "Input parameter : "+str(sys.argv[1])
print "Output parameter : "+str(sys.argv[2])
print "Individuals with more than "+str(float(sys.argv[4])*100)+"% Mendel errors will be listed."
print "SNPs with more than "+str(float(sys.argv[3])*100)+"% Mendel errors will be listed."
print "Missing code : "+str(sys.argv[5])+"\n"

ped = open(str(sys.argv[1])+".ped")
map = open(str(sys.argv[1])+".map")
varErr = float(sys.argv[3])
indErr = float(sys.argv[4])
missing_code = str(sys.argv[5])

mend = open(str(sys.argv[2])+'.mendel',"w")
lmendel = open(str(sys.argv[2])+'.lmendel',"w")
fmendel = open(str(sys.argv[2])+'.fmendel',"w")
imendel = open(str(sys.argv[2])+'.imendel',"w")
lfilt = open(str(sys.argv[2])+'.lfilt',"w")
ifilt = open(str(sys.argv[2])+'.ifilt',"w")
log = open(str(sys.argv[2])+'.log',"w")

fam=[]
id=[]
non_founder=[]
geno=[]

snp = []
for line in map:
 ele=line.split()
 if int(ele[0]) <24: snp.append(ele)

print "There are "+str(len(snp))+" SNPs. (Only autosomes and X chr are included)"

i=0
for line in ped:
 ele = line.split()
 fam.append(ele[:6])
 geno.append(ele[6:len(snp)*2+6])
 id.append(ele[1])
 if not (ele[2]=='0' and ele[3]=='0'):
  non_founder.append(i)
 i+=1


for ii in non_founder:
 if not fam[ii][2] in id: fam[ii][2] = '0'
 if not fam[ii][3] in id: fam[ii][3] = '0'

i=0
non_founder=[]
for ele in fam:
 if not (ele[2]=='0' and ele[3]=='0'):
  non_founder.append(i)
 i+=1

print "Ungenotyped parents will be set to missing."
print str(len(fam))+" individuals are detected."
print str(len(fam)-len(non_founder))+" founders and "+str(len(non_founder))+" non-founders are detected."


imenDict = {}
for i in range(len(fam)):
 imenDict[fam[i][1]]=0



fmenDict = {}
for i in range(len(fam)):
 if not fam[i][0] in fmenDict.keys():
  fmenDict[fam[i][0]]=0

print str(len(fmenDict))+" families are detected."

lmenDict = {}

mendel = ["FID\tKID\tCHR\tSNP\tCODE\tERROR"]

def isConsistent(a1,a2,c,missing_code):
  return c == str(missing_code) or a1 == str(missing_code) or a2 == str(missing_code) or c == a1 or c == a2

print "Checking mendelian inconsistency rate..."

for j in non_founder:
 parent = fam[j][2:4]
 if parent[0]=='0':
  mo = id.index(parent[1])
  geno_mo = geno[mo]
  geno_off = geno[j]
  for jj in range(0,len(geno_mo),2):
   if not geno_mo[jj] == str(missing_code) and geno_mo[jj]==geno_mo[jj+1]:
    if not geno_off[jj] == str(missing_code) and geno_off[jj]==geno_off[jj+1]:
     if not geno_mo[jj]==geno_off[jj]:
      ## .mendel
      out1 = []
      msg = ['3',"*/* x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
      out1.extend(fam[j][0:2])
      out1.extend(snp[jj/2][0:2])
      out1.extend(msg)
      ## .imendel
      imenDict[parent[1]]=imenDict.get(parent[1])+1
      imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
      ## .fmendel
      fmenDict[fam[j][0]]=fmenDict.get(fam[j][0])+1
      ## .lmendel
      lmenDict['\t'.join(snp[jj/2])] = lmenDict.get('\t'.join(snp[jj/2]),0)+1
      mendel.append('\t'.join(out1))
 elif parent[1]=='0':
  fa = id.index(parent[0])
  geno_fa = geno[fa]
  geno_off = geno[j]
  for jj in range(0,len(geno_fa),2):
   if not geno_fa[jj] == str(missing_code) and geno_fa[jj]==geno_fa[jj+1]:
    if not geno_off[jj] == str(missing_code) and geno_off[jj]==geno_off[jj+1]:
     if not geno_fa[jj]==geno_off[jj]:
      ## .mendel
      out1 = []
      msg = ['2',geno_fa[jj]+"/"+geno_fa[jj+1]+" x */* -> "+geno_off[jj]+"/"+geno_off[jj+1]]
      out1.extend(fam[j][0:2])
      out1.extend(snp[jj/2][0:2])
      out1.extend(msg)
      ## .imendel
      imenDict[parent[0]]=imenDict.get(parent[0])+1
      imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
      ## .fmendel
      fmenDict[fam[j][0]]=fmenDict.get(fam[j][0])+1
      ## .lmendel
      lmenDict['\t'.join(snp[jj/2])] = lmenDict.get('\t'.join(snp[jj/2]),0)+1
      mendel.append('\t'.join(out1))
 else :
  fa = id.index(parent[0])
  mo = id.index(parent[1])
  geno_fa = geno[fa]
  geno_mo = geno[mo]
  geno_off = geno[j]
  for jj in range(0,len(geno_fa),2):
   b1 = isConsistent(geno_fa[jj],geno_fa[jj+1],geno_off[jj],missing_code) and isConsistent(geno_mo[jj],geno_mo[jj+1],geno_off[jj+1],missing_code)
   b2 = isConsistent(geno_fa[jj],geno_fa[jj+1],geno_off[jj+1],missing_code) and isConsistent(geno_mo[jj],geno_mo[jj+1],geno_off[jj],missing_code)
   if not (b1 or b2):
    if int(snp[jj/2][0]) < 23:
     ## .mendel & .imendel 
     if geno_fa[jj:jj+2] == geno_mo[jj:jj+2]:
      if geno_off[jj] == geno_off[jj+1] :
       imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
       out1 = []
       msg = ['4',geno_fa[jj]+"/"+geno_fa[jj+1]+" x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
       out1.extend(fam[j][0:2])
       out1.extend(snp[jj/2][0:2])
       out1.extend(msg)
      else :
       imenDict[parent[0]]=imenDict.get(parent[0])+1
       imenDict[parent[1]]=imenDict.get(parent[1])+1
       imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
       out1 = []
       msg = ['1',geno_fa[jj]+"/"+geno_fa[jj+1]+" x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
       out1.extend(fam[j][0:2])
       out1.extend(snp[jj/2][0:2])
       out1.extend(msg)
     elif not geno_fa[jj] == missing_code and geno_fa[jj]==geno_fa[jj+1]:
      imenDict[parent[0]]=imenDict.get(parent[0])+1
      imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
      out1 = []
      msg = ['2',geno_fa[jj]+"/"+geno_fa[jj+1]+" x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
      out1.extend(fam[j][0:2])
      out1.extend(snp[jj/2][0:2])
      out1.extend(msg)
     else :
      imenDict[parent[1]]=imenDict.get(parent[1])+1
      imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
      out1 = []
      msg = ['3',geno_fa[jj]+"/"+geno_fa[jj+1]+" x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
      out1.extend(fam[j][0:2])
      out1.extend(snp[jj/2][0:2])
      out1.extend(msg)
     ## .fmendel
     fmenDict[fam[j][0]]=fmenDict.get(fam[j][0])+1
     ## .lmendel
     lmenDict['\t'.join(snp[jj/2])] = lmenDict.get('\t'.join(snp[jj/2]),0)+1
     mendel.append('\t'.join(out1))
    elif int(snp[jj/2][0]) == 23 and int(fam[j][4])==2:
     ## .mendel & .imendel 
     if geno_fa[jj:jj+2] == geno_mo[jj:jj+2]:
      if geno_off[jj] == geno_off[jj+1] :
       imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
       out1 = []
       msg = ['4',geno_fa[jj]+"/"+geno_fa[jj+1]+" x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
       out1.extend(fam[j][0:2])
       out1.extend(snp[jj/2][0:2])
       out1.extend(msg)
      else :
       imenDict[parent[0]]=imenDict.get(parent[0])+1
       imenDict[parent[1]]=imenDict.get(parent[1])+1
       imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
       out1 = []
       msg = ['1',geno_fa[jj]+"/"+geno_fa[jj+1]+" x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
       out1.extend(fam[j][0:2])
       out1.extend(snp[jj/2][0:2])
       out1.extend(msg)
     elif geno_fa[jj]==geno_fa[jj+1]:
      imenDict[parent[0]]=imenDict.get(parent[0])+1
      imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
      out1 = []
      msg = ['2',geno_fa[jj]+"/"+geno_fa[jj+1]+" x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
      out1.extend(fam[j][0:2])
      out1.extend(snp[jj/2][0:2])
      out1.extend(msg)
     else :
      imenDict[parent[1]]=imenDict.get(parent[1])+1
      imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
      out1 = []
      msg = ['3',geno_fa[jj]+"/"+geno_fa[jj+1]+" x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
      out1.extend(fam[j][0:2])
      out1.extend(snp[jj/2][0:2])
      out1.extend(msg)
     ## .fmendel
     fmenDict[fam[j][0]]=fmenDict.get(fam[j][0])+1
     ## .lmendel
     lmenDict['\t'.join(snp[jj/2])] = lmenDict.get('\t'.join(snp[jj/2]),0)+1
     mendel.append('\t'.join(out1))
    elif int(snp[jj/2][0]) == 23 and int(fam[j][4])==1:
     if not geno_mo[jj] == str(missing_code) and geno_mo[jj]==geno_mo[jj+1]:
      if not geno_off[jj] == str(missing_code) and geno_off[jj]==geno_off[jj+1]:
       if not geno_mo[jj]==geno_off[jj]:
        ## .mendel
        out1 = []
        msg = ['5',"*/* x "+geno_mo[jj]+"/"+geno_mo[jj+1]+" -> "+geno_off[jj]+"/"+geno_off[jj+1]]
        out1.extend(fam[j][0:2])
        out1.extend(snp[jj/2][0:2])
        out1.extend(msg)
        ## .imendel
        imenDict[parent[1]]=imenDict.get(parent[1])+1
        imenDict[fam[j][1]]=imenDict.get(fam[j][1])+1
        ## .fmendel
        fmenDict[fam[j][0]]=fmenDict.get(fam[j][0])+1
        ## .lmendel
        lmenDict['\t'.join(snp[jj/2])] = lmenDict.get('\t'.join(snp[jj/2]),0)+1
        mendel.append('\t'.join(out1))



### Output files ###

fmenRES = ["FID\tN"]
for key, value in fmenDict.iteritems():
 fmenRES.append(str(key)+"\t"+str(value))

fmenRESI = '\n'.join(fmenRES)

lmenRES = ["CHR\tSNP\tgenetic_pos\tphysical_pos\tN"]
varFilt = []
aa=lmenDict.keys()
aa.sort()
for key in aa:
 lmenRES.append(key+"\t"+str(lmenDict[key]))
 rate = float(lmenDict[key])/len(non_founder)
 if rate >= varErr : varFilt.append(str(key))

print str(len(varFilt))+" SNPs are listed due to Mendel errors"

lmenRESI = '\n'.join(lmenRES)
varFiltI = '\n'.join(varFilt)

bb=imenDict.keys()
bb.sort()
imenRES = ["FID\tIID\tN"]
indFilt=[]
for key in bb:
 nn = id.index(key)
 imenRES.append('\t'.join(fam[nn][0:2])+"\t"+str(imenDict[key]))
 if not fam[nn][2:4]==['0','0']:
  rate = float(imenDict[key])/len(snp)
  if rate >= indErr : indFilt.append('\t'.join(fam[nn][0:2]))

print str(len(indFilt))+" subjects are listed due to Mendel errors"

imenRESI = '\n'.join(imenRES)
indFiltI = '\n'.join(indFilt)

menRES = '\n'.join(mendel)

denom = len(non_founder)*len(snp)
numer = len(mendel)-1
print "The overall mendelian error rate is "+str(float(numer)/float(denom))

print "Output files are being written..."
mend.write(str(menRES))
mend.close()
lmendel.write(str(lmenRESI))
lmendel.close()
fmendel.write(str(fmenRESI))
fmendel.close()
imendel.write(str(imenRESI))
imendel.close()
lfilt.write(str(varFiltI))
lfilt.close()
ifilt.write(str(indFiltI))
ifilt.close()

s2 = datetime.datetime.now()
print "\n["+str(s2-s1)+"] The anlaysis is successfully finished at "+str(s2)+".\n"


log1 = "####################################################################################\n####################################################################################\n################ mendelCheck.py - dnjswlzz11@gmail.com (2015.07.10) ################\n####################################################################################\n####################################################################################\n\n"
log2 = "usage: python mendelCheck.py plink_input output variant_rate individual_rate missing_code\n\n"
log3 = "Input parameter : "+str(sys.argv[1])+"\nOutput parameter : "+str(sys.argv[2])+"\nIndividuals with more than "+str(float(sys.argv[4])*100)+"% Mendel errors will be listed.\nSNPs with more than "+str(float(sys.argv[3])*100)+"% Mendel errors will be listed.\nMissing code : "+str(sys.argv[5])+"\n\n"
log4 = "There are "+str(len(snp))+" SNPs. (Only autosomes and X chr are included)\nUngenotyped parents will be set to missing.\n"+str(len(fam))+" individuals are detected.\n"+str(len(fam)-len(non_founder))+" founders and "+str(len(non_founder))+" non-founders are detected.\n"+str(len(fmenDict))+" families are detected.\n"
log5 = "Checking mendelian inconsistency rate...\n"+str(len(varFilt))+" SNPs are listed due to Mendel errors\n"+str(len(indFilt))+" subjects are listed due to Mendel errors\nThe overall mendelian error rate is "+str(float(numer)/float(denom))+"\nOutput files are being written...\n\n["+str(s2-s1)+"] The anlaysis is successfully finished at "+str(s2)+".\n"

log.write(log1+log2+log3+log4+log5)
log.close()
