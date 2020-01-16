import os
import requests
from pathos.pools import ProcessPool
import re
import subprocess
#############################
encoding = 'utf-8'
#############################################https://www.cancergenomeinterpreter.org/mutations
cancer="catalog_of_validated_oncogenic_mutations.tsv"
##############################################
gdna,site={},{}
if not os.path.exists("hotspot.tsv"):
    outfile=open("hotspot.tsv","w")
    outfile.write("#Chr\tPos\tRef\tAlt\tGene\tProtein\tContext\tgdna\tInfo\n")
    outfile.close()
else:
    infile=open("hotspot.tsv","r")
    for line in infile:
        line=line.strip()
        array=line.split("\t")
        gdna[array[7]]=1
        site[array[0]+"\t"+array[1]+"\t"+array[2]+"\t"+array[3]]=1
    infile.close()
####################################
infile = open(cancer, "r")
num,var,info=0,[],[]
for line in infile:
    line=line.strip()
    array=line.split("\t")
    num+=1
    if num!=1 and re.search('somatic',line):
            string=array[0]+"\t"+array[2]+"\t"+array[5]+"\t"
            tmp=array[1].split("__")
            for key in tmp:
                if not key in gdna:
                    var.append(key)
                    info.append(string+key+"\t.")
infile.close()
def run(id,otherinfo):
    server = "https://grch37.rest.ensembl.org/vep/human/hgvs/" + id
    res = requests.get(server, headers={"Content-Type": "application/json"})
    if not res.ok:
        print(id)
    else:
        outfile = open("hotspot.tsv", "a+")
        for item in res.json():
            if item['allele_string']:
                Ref = item['allele_string'].split("/")[0]
                Alt = item['allele_string'].split("/")[1]
                Pos = item['start']
                Chr = item['seq_region_name']
                outfile.write("chr%s\t%s\t%s\t%s\t%s\n" % (Chr, Pos, Ref, Alt,otherinfo))
        outfile.close()
def civic():
    vcf=open("civic_variants.vcf","r")#https://civicdb.org/home
    outfile = open("hotspot.tsv", "a+")
    for line in vcf:
        line=line.strip()
        if not line.startswith("#") and re.search('Germline',line,re.I):
            array=line.split("\t")
            ALT = array[4].split("/")
            for key in ALT:
                tmp="chr"+array[0]+"\t"+array[1]+"\t"+array[3]+"\t"+key
                if not tmp in site:
                    string="chr"+array[0]+"\t"+array[1]+"\t"+array[3]+"\t"+key+"\t"
                    p=re.compile(r'GN=(\S+)')
                    string+=p.findall(line.split(";")[0])[0]+"\t.\t"+line.split("|")[-2]+"\t.\t.\n"
                    outfile.write(string)
    vcf.close()
    outfile.close()
def Docm():#http://www.docm.info,Curation of the literature to produce a high quality set of pathogenic somatic mutations is not trival.
    outfile = open("hotspot.tsv", "a+")
    chr=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]
    for key in chr:
        html="http://www.docm.info/api/v1/variants.json?chromosomes="+str(key)
        res = requests.get(html)
        for item in res.json():
            string="chr"+item['chr']+"\t"+str(item['start'])+"\t"+item['read']+"\t"+item['variant']
            if not string in site:
                string = "chr" + item['chr'] + "\t" + str(item['start']) + "\t" + item['read'] + "\t" + item['variant']+"\t"
                string+=item['gene']+"\t"+item['amino_acid']+"\t"+"."+"\t"+"."+"\t.\n"
                outfile.write(string)
    outfile.close()
def clinvar_cosmic():#添加clinvar和comic数据库中共有的位点，而且这些位点CNT>=50(出现的样本数)
    clinvar="/data/Database/clinvar/clinvar.vcf"
    infile=open(clinvar,"r")
    clinvar_site={}
    for line in infile:
        if not line.startswith("#"):
            line=line.strip()
            array=line.split("\t")
            tmp="chr"+array[0]+"\t"+array[1]+"\t"+array[3]+"\t"+array[4]
            if not tmp in site:
                clinvar_site[tmp]=1
    infile.close()
    cosmic="/data/Database/COSMIC/release_v88/CosmicCodingMuts.vcf"
    infile=open(cosmic,"r")
    cosmic_clinvar,CNT,gene,trans={},{},{},{}
    for line in infile:
        if not line.startswith("#") and not re.search('SNP',line):
            line=line.strip()
            array=line.split("\t")
            tmp="chr"+array[0]+"\t"+array[1]+"\t"+array[3]+"\t"+array[4]
            p=re.compile(r'CNT=(\d+)')
            if not tmp in site and tmp in clinvar_site:
                cosmic_clinvar[tmp]="chr"+array[0]+"\t"+array[1]+"\t"+array[3]+"\t"+array[4]
                if tmp in CNT:
                    CNT[tmp]+=int(p.findall(line)[0])
                else:
                    CNT[tmp] = int(p.findall(line)[0])
    infile.close()
    outfile = open("hotspot.tsv", "a+")
    for key in cosmic_clinvar:
        if CNT[key]>=50:
            outfile.write(cosmic_clinvar[key]+"\t.\t.\t.\t.\t.\n")
    outfile.close()
def PharmGKB():
    ###https://www.pharmgkb.org
    PharmGKB_var=open("/home/fanyucai/test/hotspot/clinical_ann_metadata.tsv","r")
    rs=open("rsID.txt","w")
    for line in PharmGKB_var:
        line=line.strip()
        array=line.split("\t")
        if array[1].startswith("rs"):
            rs.write("%s\n"%array[1])
    rs.close()
    PharmGKB_var.close()
    if not os.path.exists("output.avinput"):
        cmd="/software/perl/perl-v5.28.1/bin/perl /software/docker_tumor_base/Resource/Annovar/convert2annovar.pl"
        cmd+=" --format rsid rsID.txt --avsnpfile /software/docker_tumor_base/Resource/Annovar/humandb/hg19_avsnp147.txt >output.avinput"
        subprocess.check_call(cmd,shell=True)
    elif not os.path.exists("PharmGKB.hg19_multianno.txt"):
        par = " -protocol refGene,cytoBand,exac03,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,gnomad211_exome,gnomad211_genome,cosmic88_coding,clinvar_20190305,ljb26_all,intervar_20180118"
        par += ",1000g2015aug_sas,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur "
        par += " -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f "
        par += " -nastring . -polish "
        subprocess.check_call(
            "/software/perl/perl-v5.28.1/bin/perl /software/docker_tumor_base/Resource/Annovar/table_annovar.pl output.avinput /software/docker_tumor_base/Resource/Annovar/humandb -buildver hg19 -out PharmGKB -remove %s " % (
                par), shell=True)
    else:
        #################
        cnt = {}
        file1 = open("/data/Database/COSMIC/release_v88/CosmicCodingMuts.vcf", "r")
        for line in file1:
            if not line.startswith("#"):
                line = line.strip()
                array = line.split()
                pattern = re.compile(r'CNT=(\d+)')
                a = pattern.findall(line)
                cnt[array[2]] = int(a[0])
        file1.close()
        ##########################
        infile=open("PharmGKB.hg19_multianno.txt","r")
        outfile = open("hotspot.tsv", "a+")
        database = ['1000g2015aug_all', '1000g2015aug_eas', 'ExAC_ALL', 'esp6500siv2_all', 'ExAC_EAS','genome_AF','exome_AF','exome_AF_eas','genome_AF_eas']
        dict = {}
        for line in infile:
            line = line.strip()
            array = line.split("\t")
            name = []
            if line.startswith("Chr"):
                for i in range(len(array)):
                    name.append(array[i])
                    dict[array[i]] = i
            else:
                freq = 0
                freq_counts = 0
                result = ""
                dot = 0
                for i in database:
                    if array[dict[i]] == ".":
                        dot += 1
                    elif array[dict[i]] != "." and float(array[dict[i]]) <= 0.01:
                        freq += 1
                    else:
                        freq_counts += 1
                if freq_counts==0:  # not common snp
                    if array[dict['CLNSIG']].startswith("Pathogenic") or array[dict['CLNSIG']].startswith(
                            "Likely_pathogenic") or array[dict['CLNSIG']].startswith("drug_response"):
                        result = "true"
                    if array[dict['InterVar_automated']].startswith("Pathogenic") or array[
                        dict['InterVar_automated']].startswith("Likely pathogenic"):
                        result = "true"
                    if freq >= 1:  # at least 1<MAF
                        result = "true"
                    if dot == len(database):
                        result = "true"
                if result == "true":
                    tmp = "chr" + array[0] + "\t" + array[1] + "\t" + array[3] + "\t" + array[4]
                    if not tmp in site:
                        outfile.write(tmp+"\t%s\t.\t.\t.\t.\n"%(array[6]))
        outfile.close()
if __name__=="__main__":
    pool = ProcessPool(nodes=100)
    pool.map(run, var, info)
    civic()
    Docm()
    clinvar_cosmic()
    PharmGKB()
