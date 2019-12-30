#   设置环境文件

    下载https://civicdb.org/downloads/nightly/nightly-civicpy_cache.pkl
    然后在~/.bash_profile中设置：
    export CIVICPY_CACHE_FILE=/path/to/nightly-civicpy_cache.pkl

#   加载安装python3模块

    pip3 install civicpy==1.0.0rc2

#   生成CIViC数据库中所有条目的vcf格式文件

    from civicpy import civic, exports
    with open('civic_variants.vcf', 'w', newline='') as file:
        w = exports.VCFWriter(file)
        all_variants = civic.get_all_variants()
        w.addrecords(all_variants)
        w.writerecords()

#   