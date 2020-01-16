## [1.CIVIC数据](https://civicdb.org/home)
###  1.1设置环境文件

    下载https://civicdb.org/downloads/nightly/nightly-civicpy_cache.pkl
    然后在~/.bash_profile中设置：
    export CIVICPY_CACHE_FILE=/path/to/nightly-civicpy_cache.pkl

###  1.2加载安装python3模块

    pip3 install civicpy==1.0.0rc2

###  1.3生成CIViC所有条目的vcf格式文件

    from civicpy import civic, exports
    with open('civic_variants.vcf', 'w', newline='') as file:
        w = exports.VCFWriter(file)
        all_variants = civic.get_all_variants()
        w.addrecords(all_variants)
        w.writerecords()
        
**备注:在vcf文件中有对每个位点的具体分类，只提取somatic位点**

## [2.cancergenomeinterpreter数据库](https://www.cancergenomeinterpreter.org/home)

### 2.1下载变异位点
[catalog_of_validated_oncogenic_mutations_latest.zip](https://www.cancergenomeinterpreter.org/data/catalog_of_validated_oncogenic_mutations_latest.zip?ts=20180216)
### 2.2获得位点信息
	主要是API(https://grch37.rest.ensembl.org/vep/human/hgvs/)获得Chr/Pos/Ref/Alt的信息

**备注:原始文件中对每个位点的具体分类，只提取somatic位点**


## [3.Docm数据库](http://www.docm.info)

    主要是借助该数据库自己的API，按照染色体获得所有变异位点信息,好像该数据库都是人工整理的somatic位点。

##	4.[cosmic](https://cancer.sanger.ac.uk/cosmic)与[clinvar](https://www.ncbi.nlm.nih.gov/clinvar/)
	这两个数据库取交集并限制变异位点在cosmic数据库的CNT>=50，此外在cosmic数据库也就是VCF中注释为SNP的去掉

##	5.[PharmGKB数据库](https://www.pharmgkb.org)

    提取数据库中的rs号（dbsnp_v147),使用annovar的convert2annovar.pl程序转化成染色体位置以及变异信息,并对人群频率数据库进行过滤

以上代码以及相关文件：
[https://github.com/fanyucai1/hotspot](https://github.com/fanyucai1/hotspot)

##  备注

综上所述得到的变异位点取并集构成目前的第一版本的热点数据库，最后更新时间2020-1-16，欢迎批评指正！！！