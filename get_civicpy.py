from civicpy import civic, exports

with open('civic_variants.vcf', 'w', newline='') as file:
    w = exports.VCFWriter(file)
    all_variants = civic.get_all_variants()
    w.addrecords(all_variants)
    w.writerecords()