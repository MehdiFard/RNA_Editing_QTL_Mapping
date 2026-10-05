#!/usr/bin/env python3

import sys
import numpy as np

from optparse import OptionParser
    
from sklearn.decomposition import PCA
from sklearn import preprocessing

from scipy.stats import rankdata
from scipy.stats import norm

def qqnorm(x):
    n=len(x)
    a=3.0/8.0 if n<=10 else 0.5
    return(norm.ppf( (rankdata(x)-a)/(n+1.0-2.0*a) ))

# Quantile normalization + rank-based inverse normal transform (--norm quantile_int), copied from Pantry
# phenotyping/scripts/normalize_phenotypes.py
# (functions by Francois Aguet; https://github.com/PejLab/Pantry, MIT License)
def normalize_quantiles(M):
    """Quantile normalization of columns (samples) to the average empirical distribution
    (replicates R preprocessCore::normalize.quantiles; Bolstad et al. 2003)"""
    M = M.copy()
    Q = M.argsort(axis=0)
    m,n = M.shape
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n
    for i in range(n):
        dupes = np.zeros(m, dtype=int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1
        M[Q[:,i],i] = quantiles
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1
    return M

def inverse_normal_transform(M):
    """Transform rows (sites) to a standard normal distribution; ties are averaged"""
    R = np.apply_along_axis(rankdata, 1, M)
    return norm.ppf(R/(M.shape[1]+1))

def stream_table(f, ss = ''):
    fc = '#'
    while fc[0] == "#":
        fc = f.readline().strip()
        head = fc.split(ss)

    for ln in f:
        ln = ln.strip().split(ss)
        attr = {}

        for i in range(len(head)):
            try: attr[head[i]] = ln[i]
            except: break
        yield attr

def main(ratio_file, pcs=50, norm_method="quantile_int"):
    
    dic_pop, fout = {}, {}
    try: open(ratio_file)
    except: 
        sys.stderr.write("Can't find %s..exiting\n"%(ratio_file))
        sys.exit(1)
    if not ratio_file.endswith(".noXYM.txt"): # Added by Me
        out_name = ratio_file.replace(".txt", ".noXYM.txt")
    else:
        out_name = ratio_file


    sys.stderr.write("Starting...\n")
    for i in range(1,23):
        fout[i] = open(f"{out_name}.phen_chr{i}", 'w') # Added by Me
    fout_ave = open(f"{out_name}.ave", 'w') # Added by Me
    valRows, valRowsnn, geneRows = [], [], []
    finished = False
    # Open the input file properly depending on .gz or not
    header = open(ratio_file).readline().split()[1:]

    for i in fout:
        fout[i].write("\t".join(["#Chr","start", "end", "ID"]+header)+'\n')

    for dic in stream_table(open(ratio_file),' '):

        chrom = dic['chrom'].replace("chr",'')
        chr_ = chrom.split(":")[0]
        if chr_ in 'XY': continue
        NA_indices, valRow, aveReads = [], [], []
        tmpvalRow = []

        i = 0
        for sample in header:

            if sample not in dic:
                sys.stderr.write("Error: row %s has %d fields, sample %s is missing\n"%(chrom, len(dic), sample))
                sys.exit(1)
            count = dic[sample]
            num, denom = count.split('/')
            if float(denom) < 1:
                count = "NA"
                tmpvalRow.append("NA")
                NA_indices.append(i)
            else:
                # add a 0.5 pseudocount
                count = (float(num)+0.5)/((float(denom))+0.5)
                tmpvalRow.append(count) 
                aveReads.append(count)


        # If ratio is missing for over 40% of the samples, skip
        if tmpvalRow.count("NA") > len(tmpvalRow)*0.4:
            continue

        ave = np.mean(aveReads)

        # Set missing values as the mean of all values
        for c in tmpvalRow:
            if c == "NA": valRow.append(ave)
            else: valRow.append(c)

        # If there is too little variation, skip (there is a bug in fastqtl which doesn't handle cases with no variation)
        if np.std(valRow) < 0.005: continue

        chr_, s, e = chrom.split(":")
        if len(valRow) > 0:
            if chr_ in ['M', 'X', 'Y']: # Added by Me
                continue  # Skip to the next iteration # Added by Me
            chrom_int = int(chr_)
            fout[chrom_int].write("\t".join([chr_,s,e,chrom]+[str(x) for x in valRow])+'\n')
            fout_ave.write(" ".join(["%s"%chrom]+[str(min(aveReads)), str(max(aveReads)), str(np.mean(aveReads))])+'\n')

            # scale normalize
            valRowsnn.append(valRow)                
            # valRow = preprocessing.scale(valRow)
            # Rounded: mean-imputed samples should get z = 0 exactly, but float rounding gives tiny +/- noise
            # (~1e-14) that depends on sample order, so qqnorm broke these ties differently in every run
            valRow = np.round(preprocessing.scale(valRow), 10)

            valRows.append(valRow)
            geneRows.append("\t".join([chr_,s,e,chrom]))
            if len(geneRows) % 1000 == 0:
                sys.stderr.write("Parsed %s introns...\n"%len(geneRows))

    for i in fout:
        fout[i].close()
    fout_ave.close()

    if len(valRows) == 0:
        sys.stderr.write("Error: no sites passed the filters in %s (too few samples or too little variation)\n"%(ratio_file))
        sys.exit(1)

    if norm_method == "quantile_int":
        # Quantile normalization of samples, then rank-based inverse normal transform of each site across samples
        # (as in Pantry and the GTEx eQTL pipeline)
        # (starts from the same filtered, mean-imputed ratios; the per-site z-score is not used)
        sys.stderr.write("Normalization: quantile_int (quantile normalization + per-site rank-based inverse normal transform)\n")
        matrix = inverse_normal_transform(normalize_quantiles(np.array(valRowsnn, dtype=float)))
    else:
        # zscore_qqnorm: per-site z-score (above), then qqnorms on the columns (LeafCutter; Li et al. GTEx_edQTL)
        sys.stderr.write("Normalization: zscore_qqnorm (per-site z-score + per-sample quantile-quantile normalization)\n")
        matrix = np.array(valRows)
        for i in range(len(matrix[0,:])):
            matrix[:,i] = qqnorm(matrix[:,i])
        
    # write the corrected tables
    fout = {}
    for i in range(1,23):
        # fn="%s.qqnorm_chr%d"%(ratio_file,i)
        fn = "%s.qqnorm_chr%d" % (out_name, i) # Added by Me
        print(("Outputting: " + fn))
        fout[i] = open(fn,'w')
        fout[i].write("\t".join(['#Chr','start','end','ID'] + header)+'\n')
    lst = []
    for i in range(len(matrix)):
        chrom, s = geneRows[i].split()[:2]
        
        lst.append((int(chrom.replace("chr","")), int(s), "\t".join([geneRows[i]] + [str(x) for x in  matrix[i]])+'\n'))

    lst.sort()
    for ln in lst:
        fout[ln[0]].write(ln[2])
        
    # fout_run = open("%s.prepare.sh"%ratio_file,'w')
    fout_run = open("%s.prepare.sh"%out_name,'w') # Added by Me

    for i in fout:
        fout[i].close()
        fout_run.write("bgzip -f %s.qqnorm_chr%d\n"%(out_name, i)) # Added by Me
        fout_run.write("tabix -p bed %s.qqnorm_chr%d.gz\n"%(out_name, i)) # Added by Me
    fout_run.close()

    # sys.stdout.write("Use `sh %s.prepare.sh' to create index for fastQTL (requires tabix and bgzip).\n"%ratio_file)
    sys.stdout.write("Use `sh %s.prepare.sh' to create index for fastQTL (requires tabix and bgzip).\n"%out_name) # Added by Me

    if pcs>0:
        pca = PCA(n_components=pcs)                                                                                                                                                                            
        pca.fit(matrix)  
        # pca_fn=ratio_file+".PCs"
        pca_fn=out_name+".PCs" # Added by Me
        print(("Outputting PCs: " + pca_fn))
        pcafile = open(pca_fn,'w')  
        pcafile.write("\t".join(['id']+header)+'\n')
        pcacomp = list(pca.components_)
        for i in range(len(pcacomp)):
            pcafile.write("\t".join([str(i+1)]+[str(x) for x in pcacomp[i]])+'\n')

        pcafile.close()

if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [-p num_PCs] input_perind.counts.gz")
    parser.add_option("-p", "--pcs", dest="npcs", default = 50, help="number of PCs output")
    parser.add_option("-n", "--norm", dest="norm", default = "quantile_int", choices=["quantile_int", "zscore_qqnorm"],
                      help="quantile_int (default): quantile normalization + per-site rank-based inverse normal transform; "
                           "zscore_qqnorm: per-site z-score + per-sample quantile-quantile normalization")
    (options, args) = parser.parse_args()
    if len(args)==0:
        sys.stderr.write("Error: no ratio file provided...\n")
        sys.exit(1)
    main(args[0], int(options.npcs), options.norm)
    
### END