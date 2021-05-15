#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, re, sys, time
import pandas as pd
import numpy as np

def disp(text):
    print('[{}] {}'.format(time.asctime(), text), file=sys.stderr)

def Load_Reference(ifile,chroms):
    ref, cr, seq = {}, '', ''
    disp('Loading Reference Genome: {}'.format(ifile))
    for line in open(ifile):
        if line[0] == '>':
            if len(cr) > 0:
                if len(chroms)==0 or cr in chroms: ref[cr] = seq.upper()
            cr, seq = line[1:-1].split()[0], ''
        else: seq += line.strip()
    if len(chroms) == 0 or cr in chroms: ref[cr] = seq.upper()
    del seq
    return ref

def Mark_Reference(ref,refmark,CG, CHG, CHH):
    disp('Marking Reference Genome')
    for cr in ref:
        refcr, refmarkcr = ref[cr], refmark[cr]
        index = refcr.find('C', 0, len(refcr)-2)
        while index >= 0:
            if refcr[index+1] == 'G': refmarkcr[index] = CG
            elif refcr[index+2] == 'G': refmarkcr[index] = CHG
            else: refmarkcr[index] = CHH
            index = refcr.find('C', index+1, len(refcr)-2)
        index = refcr.find('G', 2, len(refcr))
        while index >= 0:
            if refcr[index-1] == 'C': refmarkcr[index] = CG
            elif refcr[index-2] == 'C': refmarkcr[index] = CHG
            else: refmarkcr[index] = CHH
            index = refcr.find('G', index+1, len(refcr))
    return refmark

def Load_One_Read(line,ref,coverage,sam_format,unique,pair,rm_dup,trim_fillin,chroms):
    col = line.split('\t')
    if sam_format:
        if line[0] == '@': return []
        flag = col[1]
        if 'u' in flag: return []
        if unique and 's' in flag: return []
        if pair and 'P' not in flag: return []
        cr, pos, cigar, seq, strand, insert = col[2], int(col[3])-1, col[5], col[9], '', int(col[8])
        if cr not in chroms: return []
        strand_index = line.find('ZS:Z:')
        if strand_index == -1: return []
        #assert strand_index >= 0, 'Missing strand information "ZS:Z:xx". Please use BSMAP output.'
        strand = line[strand_index+5:strand_index+7]
        gap_pos, gap_size = 0, 0
        while 'I' in cigar or 'D' in cigar:
            for sep in 'MID':
                try: gap_size = int(cigar.split(sep, 1)[0])
                except ValueError: continue
                break
            if sep == 'M': gap_pos += gap_size
            elif sep == 'I': seq = seq[:gap_pos] + seq[gap_pos+gap_size:]
            elif sep == 'D':
                seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
                gap_pos += gap_size
            cigar = cigar[cigar.index(sep)+1:]
    else:
        flag = col[3][:2]
        if flag == 'NM' or flag == 'QC': return []
        if unique and flag != 'UM': return []
        if pair and col[7] == '0': return []
        seq, strand, cr, pos, insert, mm = col[1], col[6], col[4], int(col[5])-1, int(col[7]), col[9]
        if cr not in chroms: return []
        if ':' in mm:
            tmp = mm.split(':')
            gap_pos, gap_size = int(tmp[1]), int(tmp[2])
            if gap_size < 0: seq = seq[:gap_pos] + seq[gap_pos-gap_size:]  # insertion on reference
            else: seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
    if pos + len(seq) >= len(ref[cr]): return []
    if rm_dup:  # remove duplicate hits
        if strand == '+-' or strand == '-+': frag_end, direction = pos+len(seq), 2
        else: frag_end, direction = pos, 1
        if coverage[cr][frag_end] & direction: return []
        coverage[cr][frag_end] |= direction
    if trim_fillin > 0: # trim fill in nucleotides
        if strand == '+-' or strand == '-+': seq = seq[:-trim_fillin]
        elif strand == '++' or strand == '--': seq, pos = seq[trim_fillin:], pos+trim_fillin
    if sam_format and insert > 0: seq = seq[:int(col[7])-1-pos] # remove overlapped regions in paired hits, SAM format only
    return (seq, strand[0], cr, pos)

def Load_Alignment(ifiles,ref,refmark,coverage,\
                   meth0,meth1,depth,meth_ct,depth_ct,\
                   sam_path,unique,pair,rm_dup,trim_fillin,chroms,\
                   seq_context,CT_SNP):
    pipes = []
    for ifile in ifiles:
        if ifile[-4:].upper() == '.SAM': sam_format, fin = True, os.popen('%ssamtools view -XS %s' % (sam_path, ifile))
        elif ifile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view -X %s' % (sam_path, ifile))
        else: sam_format, fin = False, open(ifile)
        pipes.append((ifile,sam_format,fin))
    BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}
    nmap = 0
    for ifile, sam_format, fin in pipes:
        disp('Load Alignment: {}'.format(ifile))
        nline = 0
        for line in fin:
            nline += 1
            #if nline % 10000000 == 0: disp('Read {} lines'.format(nline))
            map_info = Load_One_Read(line,ref,coverage,sam_format,unique,pair,rm_dup,trim_fillin,chroms)
            if len(map_info) == 0: continue
            seq, strand, cr, pos = map_info
            pos2 = pos + len(seq)
            nmap += 1
            match, convert, rc_match, rc_convert = BS_conversion[strand]
            refseq, refmarkcr = ref[cr], refmark[cr]
            indexes=[];methylated=False
            for index in re.finditer(match,refseq[pos:pos2]):
                index0=index.span()[0]
                indexes.append(index0)
                if seq[index0] == match and refmarkcr[index0+pos] in seq_context:
                    methylated=True
            depth_cr = depth[cr];
            meth0_cr = meth0[cr];meth1_cr = meth1[cr];
            if len(indexes)>0:
                for index in indexes:
                    if refmarkcr[index+pos] in seq_context and depth_cr[index+pos] < (2**16-1):
                        if seq[index] == convert:
                            depth_cr[index+pos] += 1
                            if methylated==True:
                                meth1_cr[index+pos] += 1
                        elif seq[index] == match:
                            depth_cr[index+pos] += 1
                            meth0_cr[index+pos] += 1
                            meth1_cr[index+pos] += 1
            if CT_SNP == 0: continue
            indexes=[]
            for index in re.finditer(rc_match,refseq[pos:pos2]):
                index0=index.span()[0]
                indexes.append(index0)
            depth_ct_cr = depth_ct[cr];
            meth_ct_cr = meth_ct[cr];
            if len(indexes)>0:
                for index in indexes:
                    if refmarkcr[index+pos] in seq_context and depth_ct_cr[index+pos] < (2**16-1):
                        if seq[index] == rc_convert:
                            depth_ct_cr[index+pos] += 1
                        elif seq[index] == rc_match:
                            depth_ct_cr[index+pos] += 1
                            meth_ct_cr[index+pos] += 1
        fin.close()
        disp('Read {} lines'.format(nline))
    if CT_SNP > 0:
        return meth0,meth1,depth,meth_ct,depth_ct,nmap
    elif CT_SNP==0:
        return meth0,meth1,depth,nmap

def Combine_Methylation_Both_Strands(ref,uncombined):
    for cr in uncombined:
        uncombined_cr, ref_cr = uncombined[cr], ref[cr]
        pos = ref_cr.find('CG')
        while pos >= 0:
            try:
                uncombined_cr[pos] += uncombined_cr[pos+1]
            except OverflowError:
                uncombined_cr[pos] = int((uncombined_cr[pos] + uncombined_cr[pos+1]) / 2)
            uncombined_cr[pos+1] = 0
            pos = ref_cr.find('CG', pos+2)
    return uncombined

def Out_CAMDA(out_prefix,wig_prefix,wig_bin,min_depth,ref,refmark,CT_SNP,seq_context,meth0,meth1,depth,meth_ct,depth_ct,nmap):
    header=['chr','pos','strand','context','ratio','eff_CT_count','C_count','CT_count','rev_G_count','rev_GA_count','CI_lower','CI_upper']
    fo_mr = open(out_prefix+"_CpG_MethRatio.tsv", 'w')
    fo_mr.write('\t'.join(header)+'\n')
    fo_camda = open(out_prefix+"_CpG_CAMDA.tsv", 'w')
    fo_camda.write('\t'.join(header)+'\n')
    if wig_prefix!=None:
        fo_mr_wig = open(wig_prefix+"_CpG_MethRatio.wig", 'w');fo_mr_wig.write('track type=wiggle_0 name='+wig_prefix+'_MethRatio\n')
        fo_camda_wig = open(wig_prefix+"_CpG_CAMDA.wig", 'w');fo_camda_wig.write('track type=wiggle_0 name='+wig_prefix+'_CAMDA\n')
        disp('Output MethRatio files and wiggle files')
    else:
        disp('Output MethRatio files')

    z95, z95sq = 1.96, 1.96 * 1.96
    nc, nd, dep0 = 0, 0, min_depth
    seq_context_str = ['CG','CHG','CHH']
    for cr in sorted(depth.keys()):
        depth_cr, meth0_cr, meth1_cr, refcr, refmarkcr = depth[cr], meth0[cr], meth1[cr], ref[cr], refmark[cr]
        if CT_SNP > 0: depth_ct_cr, meth_ct_cr = depth_ct[cr], meth_ct[cr]
        if wig_prefix!=None:
            fo_mr_wig.write('variableStep chrom={} span={}\n'.format(cr, wig_bin))
            fo_camda_wig.write('variableStep chrom={} span={}\n'.format(cr, wig_bin))
            bin = wigd = wigm0 = wigm1 = 0
        for i, dd in enumerate(depth_cr):
            if dd < dep0: continue
            if CT_SNP > 0:
                m1, d1 = meth_ct_cr[i], depth_ct_cr[i]
                if m1 != d1:
                    if CT_SNP == 2: continue
                    d = float(dd) * m1 / d1
                else: d = float(dd)
            else: d = float(dd)
            if refmarkcr[i] not in seq_context: continue
            else: seq = seq_context_str[refmarkcr[i]-1]
            if refcr[i] == 'C': strand = '+'
            else: strand = '-'
            m_mr = meth0_cr[i];m_full=meth1_cr[i];m_camda=m_full-m_mr;
            try: ratio_mr = min(m_mr, d) / d
            except ZeroDivisionError: continue
            try: ratio_full = min(m_full, d) / d
            except ZeroDivisionError: continue
            ratio_camda = ratio_full - ratio_mr
            nc += 1
            nd += d
            if wig_prefix!=None:
                if i // wig_bin != bin:
                    if wigd > 0:
                        fo_mr_wig.write('{:.0f}\t{:.3f}\n'.format(bin*wig_bin+1, min(wigm0/wigd,1)))
                        fo_camda_wig.write('{:.0f}\t{:.3f}\n'.format(bin*wig_bin+1, min((wigm1-wigm0)/wigd,1)))
                    bin = i // wig_bin#use integer division
                    wigd = wigm0 = wigm1 = 0.0
                wigd += d
                wigm0 += m_mr
                wigm1 += m_full
            denorminator = 1 + z95sq / d
            pmid_mr = ratio_mr + z95sq / (2 * d)
            sd_mr = z95 * ((ratio_mr*(1-ratio_mr)/d + z95sq/(4*d*d)) ** 0.5)
            CIl_mr, CIu_mr = (pmid_mr - sd_mr) / denorminator, (pmid_mr + sd_mr) / denorminator
            pmid_camda = ratio_camda + z95sq / (2 * d)
            sd_camda = z95 * ((ratio_camda*(1-ratio_camda)/d + z95sq/(4*d*d)) ** 0.5)
            CTl_camda, CTu_camda = (pmid_camda - sd_camda) / denorminator, (pmid_camda + sd_camda) / denorminator
            if CT_SNP>0:
                fo_mr.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\n'.format(cr, i+1, strand, seq, ratio_mr, d, m_mr, dd, m1, d1, CIl_mr, CIu_mr))
                fo_camda.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\n'.format(cr, i+1, strand, seq, ratio_camda, d, m_camda, dd, m1, d1, CTl_camda, CTu_camda))
            else:
                fo_mr.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\tNA\tNA\t{:.3f}\t{:.3f}\n'.format(cr, i+1, strand, seq, ratio_mr, d, m_mr, dd, CIl_mr, CIu_mr))
                fo_camda.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.2f}\t{}\t{}\tNA\tNA\t{:.3f}\t{:.3f}\n'.format(cr, i+1, strand, seq, ratio_camda, d, m_camda, dd, CTl_camda, CTu_camda))
    fo_mr.close();fo_camda.close();
    if wig_prefix!=None:fo_mr_wig.close();fo_camda_wig.close();
    disp('Total {} valid mappings, {} covered cytosines, average coverage: {:.2f} fold.'.format(nmap, nc, float(nd)/nc))

def read_methy_files(ifile, cols=[0,1,2,6,7]):
    names = ['chr', 'pos', 'strand', 'methy', 'total']
    disp('Loading MethRatio: {}'.format(ifile))
    meth_file = pd.read_csv(ifile, sep='\t', header=0, usecols=cols, names=names, compression='infer',low_memory=False)
    meth_file.index = meth_file['pos']
    meth_file.drop(['pos'], axis=1, inplace=True)
    return meth_file

def merge_strand_each_chr(df):
    df_p = df[df['strand']=='+']
    df_n = df[df['strand']=='-']
    df_n.index =  df_n.index.values - 1
    merge_index = np.sort(np.unique(np.append(df_n.index.values, df_p.index.values)))
    df_merge = pd.DataFrame(np.zeros(shape=[len(merge_index),2]), index=merge_index)
    df_merge.loc[df_p.index,:] = df_merge.loc[df_p.index,:] + df_p.loc[:,['methy','total']].values
    df_merge.loc[df_n.index,:] = df_merge.loc[df_n.index,:] + df_n.loc[:,['methy','total']].values
    df_merge.columns = ['methy','total']
    df_merge = df_merge.loc[0:,:] # remove the minus index pos -1
    return df_merge

def merge_strand(df):
    chs = df["chr"].unique().tolist()
    #chs = ['chr' + str(i) for i in range(1, 23)]  + ['chrX','chrY']
    df_merge = pd.DataFrame()
    for ch in chs:
        chr_sub = df[df["chr"] == ch]
        if chr_sub.shape[0] > 0:
            #disp('Merging MethRatio on {}'.format(ch))
            chr_sub = merge_strand_each_chr(chr_sub)
            chr_sub['chr']=pd.Series([ch] * chr_sub.shape[0], index=chr_sub.index)
            df_merge=df_merge.append(chr_sub)
    return df_merge

def Region_Meth_Ratio(ratio_sub,start=0,end=0):
    #ratio_sub: methratio df of one chrom
    region_meth=ratio_sub.loc[start:end,:]
    count_C=region_meth.shape[0]
    if count_C > 0:
        region_meth=merge_strand(df=region_meth)
        methy_C=region_meth["methy"].sum()
        total_C=region_meth["total"].sum()
        region_methratio=methy_C*1.0/total_C
    else:
        region_methratio=np.nan
    return [count_C,region_methratio]

def Bam2ReadCT(ifiles,ref,refmark,coverage,sam_path,unique,pair,rm_dup,trim_fillin,seq_context,chroms,output):
    pipes = []
    for ifile in ifiles:
        if ifile[-4:].upper() == '.SAM': sam_format, fin = True, os.popen('%ssamtools view -XS %s' % (sam_path, ifile))
        elif ifile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view -X %s' % (sam_path, ifile))
        else: sam_format, fin = False, open(ifile)
        pipes.append((ifile,sam_format,fin))
    BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}
    complement={"A":"T","T":"A","C":"G","G":"C"}
    nmap = 0
    fo_readct = open(output, 'w')
    for ifile, sam_format, fin in pipes:
        disp('Load Alignment: {}'.format(ifile))
        nline = 0; n_NotMatch = 0;
        for line in fin:
            nline += 1
            map_info = Load_One_Read(line,ref,coverage,sam_format,unique,pair,rm_dup,trim_fillin,chroms)
            if len(map_info) == 0: continue
            seq, strand, cr, pos = map_info
            pos2 = pos + len(seq)
            nmap += 1
            match, convert, rc_match, rc_convert = BS_conversion[strand]
            refseq, refmarkcr = ref[cr], refmark[cr]
            indexes=[]
            for index in re.finditer(match,refseq[pos:pos2]):
                index0=index.span()[0]
                indexes.append(index0)
            CT_pos_seq={};
            if len(indexes)>0:
                for index in indexes:
                    if refmarkcr[index+pos] in seq_context:
                        if seq[index] in [match, convert]:
                            if strand=="+":
                                CT_pos_seq[index+pos]=seq[index]
                            else:
                                CT_pos_seq[index+pos]=complement[seq[index]]
            if CT_pos_seq=={}:
                #fo_readct.write("\t".join([cr,strand,"NA","NA","0","NA","NA"])+"\n")
                n_NotMatch += 1
            else:
                CT_pos=sorted(CT_pos_seq);CT_seq=[];
                for key in CT_pos:
                    CT_seq.append(CT_pos_seq[key])
                CT_pos=[x+1 for x in CT_pos]
                fo_readct.write("\t".join([cr,strand,str(CT_pos[0]),str(CT_pos[-1]),str(len(CT_pos)),",".join(map(str,CT_pos)),",".join(CT_seq)])+"\n")
        fin.close()
    fo_readct.close()
    disp('Read {} lines'.format(nline))
    disp('Total {} valid mappings, {} of them donot cover cytosine in requested context(-x option)'.format(nmap,n_NotMatch))
