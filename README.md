##   1.civic数据库的处理
###  1.1设置环境文件

    下载https://civicdb.org/downloads/nightly/nightly-civicpy_cache.pkl
    然后在~/.bash_profile中设置：
    export CIVICPY_CACHE_FILE=/path/to/nightly-civicpy_cache.pkl

###  1.2加载安装python3模块

    pip3 install civicpy==1.0.0rc2

###  1.3生成CIViC数据库中所有条目的vcf格式文件

    from civicpy import civic, exports
    with open('civic_variants.vcf', 'w', newline='') as file:
        w = exports.VCFWriter(file)
        all_variants = civic.get_all_variants()
        w.addrecords(all_variants)
        w.writerecords()

##   2.cancergenomeinterpreter 数据库对处理

    主要是API(https://grch37.rest.ensembl.org/vep/human/hgvs/)获得Chr/Pos/Ref/Alt的信息

##   3.Docm数据库的处理

    主要是借助该数据库自己的API，按照染色体获得所有变异位点信息

##   4.去cosmic编码区与clinvar的交集，并限制变异位点在cosmic数据库的CNT>=50

    这一部分主要是借助cosmic和clinvar数据库的vcf文件

##  5.考虑PharmGKB数据库

    提取数据库中的rs号（dbsnp_v147),使用annovar的convert2annovar.pl程序转化成染色体位置以及变异信息

综上所述得到的变异位点取并集构成目前的第一版本的热点数据库，最后更新时间2019-12-30，欢迎批评指正！！！