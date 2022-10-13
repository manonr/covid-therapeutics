# owned by Richard Myers

#!/usr/bin/env python
'''
module load python/3.6
Don't run on the head nodes!!!!

rem - 2020-04-22
'''
import sys
import os
import re
import argparse
from read_seq import *
###############################################################################
class CodeFrame(object):
    def __init__(self, ffile, dfile, ref, incx):
        self.seqs = self.__proc_seqs(ffile)
        self.cfrs = self.__proc_def(dfile)
        self.ref = ref
        self.incx = incx
        self.tabd = {}
        self.misd = {}

    def __proc_seqs(self, ffile):
        if not os.path.isfile(ffile):
            print("Cannot open file: {}".format(ffile))
            exit()
        seqs = read_seqs(ffile)
        for fseq in seqs:
            regex = re.compile(r'[\r]')
            fseq.sid = regex.sub('', fseq.sid)
            fseq.seq = regex.sub('', fseq.seq)
            fseq.seq = fseq.seq.upper()
        return seqs

    def __proc_def(self, dfile):
        if not os.path.isfile(dfile):
            print("Cannot open file: {}".format(dfile))
            exit()
        fhl = open(dfile, 'r')
        lines = fhl.read().split('\n')
        cfrs = []
        for line in lines:
            tmp = line.split('\t')
            if len(tmp) >= 3:
                cfrs.append(line)
        fhl.close()
        return cfrs

    def __translate(self, dseq):
        aam = {'TTT' : 'F', 'TTC' : 'F', 'TTA' : 'L', 'TTG' : 'L',
               'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L',
               'ATT' : 'I', 'ATC' : 'I', 'ATA' : 'I', 'ATG' : 'M',
               'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V',
               'TCT' : 'S', 'TCC' : 'S', 'TCA' : 'S', 'TCG' : 'S',
               'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P',
               'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T',
               'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A',
               'TAT' : 'Y', 'TAC' : 'Y', 'TAA' : '*', 'TAG' : '*',
               'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q',
               'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG' : 'K',
               'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E',
               'TGT' : 'C', 'TGC' : 'C', 'TGA' : '*', 'TGG' : 'W',
               'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R',
               'AGT' : 'S', 'AGC' : 'S', 'AGA' : 'R', 'AGG' : 'R',
               'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G'}
        aseq = ""
        for i in range(0, len(dseq), 3):
            cod = "{}{}{}".format(dseq[i], dseq[i+1], dseq[i+2])
            aac = 'X'
            if cod in aam:
                aac = aam[cod]
            aseq = "{}{}".format(aseq, aac)
        return aseq

    def __comp_ref(self, raas, cfr):
        #print(raas)
        if self.ref not in raas:
            print("cannot find {} in sequence set".format(self.ref))
        else:
            refa = raas[self.ref]
            fout = "var_{}.tsv".format(cfr)
            hhdl = open(fout, 'w')
            for aas in raas:
                if aas != self.ref:
                    if aas not in self.tabd:
                        self.tabd[aas] = {}
                        self.misd[aas] = {}
                    if cfr not in self.tabd[aas]:
                        self.tabd[aas][cfr] = ""
                        self.misd[aas][cfr] = ""
                    hhdl.write("{}\t".format(aas))
                    for i in range(0, len(refa)):
                        if raas[aas][i] != refa[i]:
                            if self.incx == 'yes' or raas[aas][i] != 'X':
                                hhdl.write("{}{}{},".format(refa[i], i + 1, raas[aas][i]))
                                self.tabd[aas][cfr] = "{},{}{}{}".format(self.tabd[aas][cfr],
                                                                         refa[i],
                                                                         i + 1,
                                                                         raas[aas][i])
                            else:
                                self.misd[aas][cfr] = "{},{}{}{}".format(self.misd[aas][cfr],
                                                                         refa[i],
                                                                         i + 1,
                                                                         raas[aas][i])
                    hhdl.write('\n')
            hhdl.close()

    def __print_tabd(self):
        ohl = open("summary_vars.tsv", 'w')
        ohl.write("COG ID\t")
        for cfr in self.cfrs:
            ohl.write("{}\t".format(cfr.split('\t')[0]))
        ohl.write('\n')
        for cid, cfrs in self.tabd.items():
            ohl.write("{}\t".format(cid))
            for cfr in self.cfrs:
                cfrl = cfr.split('\t')[0]
                ohl.write("{}\t".format(cfrs[cfrl][1:]))
            ohl.write('\n')
        ohl.close()

    def __proc_missing(self):
        self.cfrs.append('missing')
        for cid, cfrs in self.misd.items():
            mstr = ""
            for cfr, missing in cfrs.items():
                locs = []
                #print("{}\t{}\t{}".format(cid, cfr, missing))
                tmp = missing.split(',')
                for miss in tmp:
                    pos = miss[1:-1]
                    if pos != "":
                        locs.append(int(pos))
                    #print("{}\t{}\t{}".format(cid, cfr, pos))
                ostr = ""
                if len(locs) > 0:
                    #print(locs)
                    start = locs[0]
                    stop = locs[0]
                    curr = locs[0]
                    for pos in locs:
                        if pos > curr + 1:
                            ostr = "{},{}-{}".format(ostr, start, stop)
                            start = pos
                            curr = pos
                        if pos == curr + 1:
                            stop = pos
                            curr = pos
                    ostr = "{},{}-{}".format(ostr, start, stop)
                    ostr = ostr[1:]
                if ostr != "":
                    ostr = "{}:{}".format(cfr, ostr)
                    #print("{}\t{}\t{}".format(cid, cfr, ostr))
                    mstr = "{},{}".format(mstr, ostr)
            #print("{} {}".format(cid, mstr[1:]))
            self.tabd[cid]['missing'] = mstr


    def gen_nuc(self):
        for cfr in self.cfrs:
            tmp = cfr.split('\t')
            if len(tmp) >= 3:
                raas = {}
                fout = "cds_{}.fas".format(tmp[0])
                gout = "aa_{}.fas".format(tmp[0])
                fhl = open(fout, 'w')
                ghl = open(gout, 'w')
                st1 = int(tmp[1]) - 1
                ed1 = int(tmp[2])
                for cid in self.seqs:
                    oseq = ""
                    if len(tmp) == 3:
                        oseq = cid.seq[st1:ed1]
                    elif len(tmp) == 5:
                        st2 = int(tmp[3]) - 1
                        ed2 = int(tmp[4])
                        oseq = "{}{}".format(cid.seq[st1:ed1], cid.seq[st2:ed2])
                    aseq = self.__translate(oseq)
                    raas[cid.sid] = aseq
                    fhl.write(">{}\n{}\n".format(cid.sid, oseq))
                    ghl.write(">{}\n{}\n".format(cid.sid, aseq))
                if self.ref != "none":
                    self.__comp_ref(raas, tmp[0])
                ghl.close()
                fhl.close()
        if self.ref != "none":
            self.__proc_missing()
            self.__print_tabd()
###############################################################################
def get_args():
    """
    parse command line arguments
    input: sys.args
    output: args hash
    """
    args = argparse.ArgumentParser(description="get_proteins.py: process climb aln into coding seq requires python/3.6")
    args.add_argument("--fasta", "-f", required = True, help = "path to fastaseq")
    args.add_argument("--def", "-d", required = True, help = "path coding defs")
    args.add_argument("--ref", "-r", default = "none", help = "name of refseq in sequences")
    args.add_argument("--incx", "-x", default = "yes", help = "include X - yes/no default yes")
    return args

#Main##########################################################################
def main():
    """main function"""
    args = vars(get_args().parse_args())
    cfr = CodeFrame(args['fasta'], args['def'], args['ref'], args['incx'])
    cfr.gen_nuc()
    
###############################################################################
if __name__ == "__main__":
    exit(main()) 
