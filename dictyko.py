#!/usr/bin/env python
"""
dictyko.py
Last update: March 3, 2015
Author: Rafael D. Rosengarten

The goal of this program is to automate the design of
Dictyostelium knock-out vector constructs.

"""
import argparse
import os
import sys

from Bio import Entrez
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqFeature import SeqFeature, FeatureLocation

from primer3 import *
from validate_email import validate_email


parser = argparse.ArgumentParser(description='Automate pLPBLP knockout vectors for Dictyostelium')
parser.add_argument('genes', help='Enter single gene name or .txt file with list of names')
parser.add_argument('--email', '-e', required=True, help='Required to query NCBI databases.')
args = parser.parse_args()


Entrez.email = args.email
is_valid = validate_email(Entrez.email)
if is_valid is False:
    sys.exit("Invalid email. Please try again.")

amb = IUPACAmbiguousDNA()

def fetch_gene_coordinates(gene):
    """Get the genomic coordinates from NCBI for the gene of interest.

    http://www.biopython.org/pipermail/biopython/2010-December/006922.html
    https://gist.github.com/chapmanb/727625

    """
    fetch_log = "Fetching_gene_coordinates for gene " + gene + "\n"
    log.write(fetch_log)
    search_term = gene + " AX4"
    handle = Entrez.esearch(db="gene", term=search_term)
    rec = Entrez.read(handle)
    gene_id = rec["IdList"][0] # assuming best match works
    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
    rec = Entrez.read(handle)[0]
    gene_locus = rec["Entrezgene_locus"][4]  ###had to modify this line to find appropriate records
    region = gene_locus["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]
    ORF_start = int(region["Seq-interval_from"]) + 1
    ORF_end = int(region["Seq-interval_to"]) + 1
    gi_id = region["Seq-interval_id"]["Seq-id"]["Seq-id_gi"]
    strand = region["Seq-interval_strand"]["Na-strand"].attributes["value"]
    return gi_id, ORF_start, ORF_end, strand

###subroutine to actually contact NCBI and get genbank
def get_gb_seq(gi_id, ORF_start, ORF_end, strand):
    get_gb_log = "get_gb_seq for gene " + gene + "\n"
    log.write(get_gb_log)
    strand = 2 if strand.lower() == "minus" else 1
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id=gi_id,
            seq_start=ORF_start, seq_stop=ORF_end, strand=strand)
    return handle.read()

###subroutine to define locus encompassing gene of interest
def locus_maps(gene, flank):
    locus_log = "locus_maps " + gene + "\n"
    log.write(locus_log)
    gi_id, ORF_start, ORF_end, strand = fetch_gene_coordinates(gene)
    ORF = get_gb_seq(gi_id, ORF_start, ORF_end, strand)    # ORF
    if (ORF_start - flank) < 0:
        locus_start = 0
    else:
        locus_start = (ORF_start - flank)
    locus_end = ORF_end + flank   ### will need conditional regarding  end of chromosome
    locusgb = get_gb_seq(gi_id, locus_start, locus_end, strand) # ORF +/- 5kb

    # write entire gb file as temp, to be parsed into simpler version
    with open("temp.gb", "w") as temp:
        temp.write(locusgb)
    outfile = gene + "_locus.gb" ## set unique outfile name

    # trim genbank feature table to include only genes (for the gene name) and CDS's for the exon/intron positions
    for record in SeqIO.parse("temp.gb","genbank"):
        record.name = gene
        check = 0
        record.features = [f for f in record.features if (f.type == "CDS" or f.type == "gene")]
        for feat in record.features:                # we need to deal with instances for predicted genes
            if (feat.type == "gene"):                # that do not have a "gene" feature qualifier,
                if 'gene' not in feat.qualifiers:    # and thus remain unlabeled on the resulting map
                    feat.qualifiers["gene"] = feat.qualifiers["locus_tag"]

    # error checking in case the NCBI search went wrong
    # this would ideally be done in the fetch_gene_coordinates module,
    # but it wasn't clear to me how to tell if the search went awry.
    # here I simply check the genbank feature table and if the gene of interest is missing
    # we report and log an error.
                if (gene in feat.qualifiers["gene"]) or (gene in feat.qualifiers["locus_tag"]):
                    check = check + 1
        if (check>0):
            SeqIO.write(record, newdir+"/"+outfile, "genbank")
            confirm = "Writing genbank file for: " + gene + "\n"
            log.write(confirm)
            print confirm
            return outfile
        else:
            error = "Error: Could not accurately identify gene " + gene + "\n"
            print error
            log.write(error)
            return "error"

###subroutine to run primer3
def runp3(geneHA, sequence, HApos, HA):
    p3_log = "runp3 for " + geneHA + "\n"
    log.write(p3_log)
    if (HA == 5):
        pos1 = HApos - 1500
    elif (HA == 3):
        pos1 = HApos - 100

    #modified from:     http://benpruitt.github.io/primer3-py/quickstart.html#primer-design
    p3output = bindings.designPrimers(
        {
            'SEQUENCE_ID': geneHA,
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': [pos1,1600],
            'SEQUENCE_TARGET': [HApos,1]
        },
        {
            'PRIMER_OPT_SIZE': 25,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 30,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 45.0,
            'PRIMER_MAX_TM': 70.0,
            'PRIMER_MIN_GC': 18.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [800,1500],
        })
    ## will return Error if p3 fails
    if (p3output['PRIMER_PAIR_NUM_RETURNED'] < 1):
        print "Failed to find primers for ", geneHA, "\n"
        log.write("Failed to find primers for "+geneHA+"\n")
        return "error"
    else:
        return p3output

###subroutine to parse primer3 output
def parsep3(p3output, geneHA):
    parse_log = "parsep3 for " + geneHA + "\n"
    log.write(parse_log)
##parse relevant info from output dictionary
    Lprimer = str(p3output['PRIMER_LEFT_0_SEQUENCE'])
    Ltm = round(p3output['PRIMER_LEFT_0_TM'])
    Ltuple = p3output['PRIMER_LEFT_0']
    Lstart, Lleng = str(Ltuple[0]), str(Ltuple[1])
    Rprimer = str(p3output['PRIMER_RIGHT_0_SEQUENCE'])
    Rtm = round(p3output['PRIMER_RIGHT_0_TM'])
    Rtuple = p3output['PRIMER_RIGHT_0']
    Rstart, Rleng = str(Rtuple[0]), str(Rtuple[1])
    HAleng = str(int(Rstart) - int(Lstart))
##append output table
    with open('gene_specific_primers.csv', 'a') as pf:
        #print "writing primers for ", geneHA
        log.write("writing primers for "+geneHA+"\n")
        pf.write(geneHA+",")
        pf.write(Lprimer+",")
        pf.write(str(Ltm)+",")
        pf.write(Lstart+",")
        pf.write(Lleng+",")
        pf.write(Rprimer+",")
        pf.write(str(Rtm)+",")
        pf.write(Rstart+",")
        pf.write(Rleng+",")
        pf.write(HAleng+"\n")

    return (Lstart, Lleng, Rstart, Rleng, Ltm, Rtm, Lprimer, Rprimer)


###subroutine to search Homology Arm amplicons for RE sites and to identify usable sites in pLPBLP
def HA_REs(HA5str, HA3str):
    HARE_log = "HA_REs for gene " + gene + "\n"
    log.write(HARE_log)
    HA5 = Seq(HA5str, amb)
    HA3 = Seq(HA3str, amb)
##assign multiple cloning site enzymes #this would have to be modified for a different vector
    mcsA = [KpnI,SalI,HindIII]
    mcsB = [SpeI,BamHI,PstI]
    typeIIs = [BbsI, BsaI]
    absentA5, absentA3, absentB5, absentB3, absentTypeII = [], [], [], [], []
#initialize markers
    A5_B3 = 0
    B5_A3 = 0
    error = 0
#look for typeIIs enzymes and keep those that are absent
    for enz in typeIIs:
        if (len(enz.search(HA5)) == 0) and (len(enz.search(HA3)) == 0):
            absentTypeII.append(enz)
    if len(absentTypeII) > 0:
        typeII = absentTypeII[0]
    else:
        typeII = "no"
#look for cloning site enzymes, as above
    for enz in mcsA:
        if (len(enz.search(HA5)) == 0):
            absentA5.append(enz)
        if (len(enz.search(HA3)) == 0):
            absentA3.append(enz)
    for enz in mcsB:
        if (len(enz.search(HA5)) == 0):
            absentB5.append(enz)
        if (len(enz.search(HA3)) == 0):
            absentB3.append(enz)
#determine which homology arm goes into which multiple cloning site (A or B)
#default to putting HA5' into MCSA
    if ((len(absentA5) > 0) and (len(absentB3) > 0)):
        A5_B3 = 1
        if (len(absentA5) > 1):
            addRE55 = absentA5[0]
            addRE53 = absentA5[1]
        elif (len(absentA5) == 1):
            addRE55, addRE53 = absentA5[0], absentA5[0]
        if (len(absentB3) > 1):
            addRE33 = absentB3[0]
            addRE35 = absentB3[1]
        elif (len(absentB3) == 1):
            addRE35, addRE33 = absentB3[0], absentB3[0]
    ##do I need to escape conditional once this is satisfied?
#otherwise put HA5' into MCSB
    elif ((len(absentB5) > 0) and (len(absentA3) > 0)):
        B5_A3 = 1
        if (len(absentB5) > 1):
            addRE53 = absentB5[0]
            addRE55 = absentB5[1]
        elif (len(absentB5) == 1):
            addRE55, addRE53 = absentB5[0], absentB5[0]
        if (len(absentA3) > 1):
            addRE35 = absentA3[0]
            addRE33 = absentA3[1]
        elif (len(absentA3) == 1):
            addRE35, addRE33 = absentA3[0], absentA3[0]
    else:
        print "Could not find appropriate combination of unique restriction sites. Inspect manually.\n"
        error = 1
#set markers for future modules
    if A5_B3 == 1:
        HAorder = "A5_B3"
    elif B5_A3 == 1:
        HAorder = "B5_A3"

    return error,addRE55,addRE53,addRE35,addRE33,typeII, HAorder

###subroutine to take subsection of GB file and maintain partial annotations (that would be automatically lost)
def GBslice(Record, HAseq):
    slice_log = "GBslice for gene " + gene + "\n"
    log.write(slice_log)
    fullseq = Record.seq
    sl_st = fullseq.find(HAseq)
    sl_end = sl_st + len(HAseq)
    slicerecord = Record[sl_st:sl_end]

#define and add annotations for features that would be lost due to partial overlap with slice ends
    for feature in Record.features:
        if (sl_st in feature) and (feature.type == "gene"):
            exact = feature.location.nofuzzy_end
            slicefeature1_end = exact - sl_st
            slicefeature1_loc = FeatureLocation(0, slicefeature1_end)
            slicefeature1_type = feature.type
            slicefeature1_qualifier = {'gene':feature.qualifiers["gene"]}
            slicefeature1_strand = feature.strand
            slicefeature1 = SeqFeature(slicefeature1_loc,type=slicefeature1_type,strand=slicefeature1_strand,qualifiers=slicefeature1_qualifier)
            slicerecord.features.append(slicefeature1)
#define and add annotations for features that would be lost due to partial overlap with slice ends
        if (sl_end in feature) and (feature.type == "gene"):
            exact = feature.location.nofuzzy_start
            slicefeature2_start = len(slicerecord) - (sl_end - exact)
            slicefeature2_loc = FeatureLocation(slicefeature2_start, len(slicerecord))
            slicefeature2_type = feature.type
            slicefeature2_qualifier = {'gene':feature.qualifiers["gene"]}
            slicefeature2_strand = feature.strand
            slicefeature2 = SeqFeature(slicefeature2_loc,type=slicefeature2_type,strand=slicefeature2_strand,qualifiers=slicefeature2_qualifier)
            slicerecord.features.append(slicefeature2)

    return slicerecord

###subroutine to add restriction enzymes to appropriate ends of HAs, including typeIIS
def addREprimers(frag, RE5, RE3, Lleng, Rleng, typeII, HAorder, HAnum):
    addREprime_log = "addREprimers for gene " + gene + "\n"
    log.write(addREprime_log)

##concatenate the RE sites plus 'aaa' triplet to ends of contig
    if typeII == "no":
        typeII_seq = ''
        revII = ''
    else:
        typeII_seq = typeII.site
        revII = Seq(typeII.site, amb).reverse_complement()

#use the HA-MCS order markers from before to add the appropriate enzyme sites to the appropriate sides of the HAs
    if HAorder == "A5_B3":
        if HAnum == 5:
            add5 = Seq(('aaa' + RE5.site + typeII_seq).lower(), amb)
            add3 = Seq((RE3.site + 'aaa').lower(), amb)
            tII_strand = 1
        elif HAnum == 3:
            add5 = Seq(('aaa' + RE5.site).lower(), amb)
            add3 = Seq((str(revII) + RE3.site + 'aaa').lower(), amb)
            tII_strand = -1
    elif HAorder == "B5_A3":
        frag = frag.reverse_complement()    ##flips the HAs to accommodate locations on opposite sites of BsR marker
        Ltmp = Lleng    #switch primer length values
        Rtmp = Rleng
        Lleng = Rtmp
        Rleng = Ltmp

        if HAnum == 3:
            add5 = Seq(('aaa' + RE5.site + typeII_seq).lower(), amb)
            add3 = Seq((RE3.site + 'aaa').lower(), amb)
            tII_strand = 1
        elif HAnum == 5:
            add5 = Seq(('aaa' + RE5.site).lower(), amb)
            add3 = Seq((str(revII) + RE3.site + 'aaa').lower(), amb)
            tII_strand = -1

    concat = (add5 + frag[:len(frag)] + add3)
    catseq = concat.seq
    if tII_strand == 1:
        typeII_start = 9
    else:
        typeII_start = len(concat) - len(add3)
    typeII_end = typeII_start + len(typeII_seq)

##make features for the gene-specific primer (gsp)
    gspF_location = FeatureLocation(len(add5), (len(add5)+int(Lleng)))
    gspR_location = FeatureLocation((len(concat)-(len(add3)+int(Lleng))), (len(concat)-len(add3)))
    labelF = "HA"+str(HAnum)+"_gspFwd"
    labelR = "HA"+str(HAnum)+"_gspRev"
    gspF_qual = {'label':labelF}
    gspR_qual = {'label':labelR}
    gspF_feature = SeqFeature(gspF_location,type='primer_bind',strand=1,qualifiers=gspF_qual)
    gspR_feature = SeqFeature(gspR_location,type='primer_bind',strand=-1,qualifiers=gspR_qual)
    concat.features.append(gspF_feature)
    concat.features.append(gspR_feature)

###make features for primers to order
    Fprime_loc = FeatureLocation(0, (len(add5)+int(Lleng)))
    Rprime_loc = FeatureLocation((len(concat)-(len(add3)+int(Lleng))), len(concat))
    Fprime_qual = {'label':'primer_Fwd'}
    Rprime_qual = {'label':'primer_Rev'}
    Fprime_feature = SeqFeature(Fprime_loc,type='primer',strand=1,qualifiers=Fprime_qual)
    Rprime_feature = SeqFeature(Rprime_loc,type='primer',strand=-1,qualifiers=Rprime_qual)
    concat.features.append(Fprime_feature)
    concat.features.append(Rprime_feature)

##make features for typeII
    if typeII != "no":
        typeII_loc = FeatureLocation(typeII_start, typeII_end)
        typeII_qual = {'label':'typeIIs'}
        typeII_feature = SeqFeature(typeII_loc,type='misc_bind',strand=tII_strand,qualifiers=typeII_qual)
        concat.features.append(typeII_feature)

    fwd = concat.seq[0:(len(add5)+int(Lleng))]
    rev = (concat.seq[(len(concat)-(len(add3)+int(Lleng))):len(concat)]).reverse_complement()
    return concat, fwd, rev

###subroutine to write GB file with HA's inserted into pLPBLP vector
##requires pLPBLP vector present in working directory
def pLPBLP_merge(HA5, HA3, RE55, RE53, RE35, RE33, HAorder):
    merge_log = "pLPBLP_merge for gene " + gene + "\n"
    log.write(merge_log)
    pLPBLP = SeqIO.read("pLPBLP.gb", "genbank")

#determine the insertional distances based on the orientation of the HAs relative to the MCSs
    if HAorder == "A5_B3":
        A_HA = HA5
        B_HA = HA3
        #use the string.find notation rather than RE.search to get start position of motif rather than cut site
        LPA5 = (str(pLPBLP.seq)).find(RE55.site)
        LPA3 = (str(pLPBLP.seq)).find(RE53.site)
        dist35 = (str(pLPBLP.seq)).find(RE33.site)
        mcaA_remove = (6 + LPA3) - LPA5
#same as above, except use the RE3 sites instead of RE5 sites to set the mcsA digestion
    elif HAorder == "B5_A3":
        A_HA = HA3
        B_HA = HA5
        LPA5 = (str(pLPBLP.seq)).find(RE35.site)
        LPA3 = (str(pLPBLP.seq)).find(RE33.site)
        dist35 = (str(pLPBLP.seq)).find(RE53.site)
        mcaA_remove = (6 + LPA3) - LPA5

#slice up pLPBLP GB
#shift position 0 to first RE in mcsA
    pLPBLP_shiftA = pLPBLP[LPA5:]+pLPBLP[:LPA5]
    mcsA_digest = pLPBLP_shiftA[mcaA_remove:]
#add arm into mcsA

    concatA = A_HA[3:(len(A_HA)-3)] + mcsA_digest

#add arm into mcsB
    if HAorder == "A5_B3":
        LPB5 = (str(concatA.seq)).find(RE35.site)
        LPB3 = (str(concatA.seq)).find(RE33.site)
    if HAorder == "B5_A3":
        LPB5 = (str(concatA.seq)).find(RE55.site)
        LPB3 = (str(concatA.seq)).find(RE53.site)

    mcaB_remove = (6 + LPB3) - LPB5
    concat_shiftB = concatA[LPB5:]+concatA[:LPB5]
    mcsB_digest = concat_shiftB[mcaB_remove:]

#I would like to add a bit to write this as 'circular' DNA
    merged_map = B_HA[3:(len(B_HA)-3)] + mcsB_digest

    return merged_map

###subroutine to add HA annotations to final GB insertion map

def build_final_map (locus, KO, Fseq5, Fleng5, Rseq3, Rleng3, HAorder):
##reverse complement right primers
    Fseq5rc = str((Seq(Fseq5, amb)).reverse_complement())
    Rseq3rc = str((Seq(Rseq3, amb)).reverse_complement())

#get ends of locus positions
    locLeft = (locus.seq).find(Fseq5) - 1
    locRight = (locus.seq).find(Rseq3rc) + len(Rseq3rc)

    leftRec = GBslice(locus, (locus[:locLeft]).seq)
    rightRec = GBslice(locus, (locus[locRight:]).seq)

##get annotated insert from KO plasmid
#shift position of map so position 0 outside of insert
    if HAorder == "A5_B3":
        shiftpos = (KO.seq).find(Rseq3rc) + int(Rleng3)
        KOshift = KO[shiftpos:]+KO[:shiftpos]
        inStart = (KOshift.seq).find(Fseq5)
        inEnd = (KOshift.seq).find(Rseq3rc) + int(Rleng3)

    elif HAorder == "A3_B5":
        shiftpos = (KO.seq).find(Fseq5rc) + int(Fleng5)
        KOshift = KO[shiftpos:]+KO[:shiftpos]
        inStart = (KOshift.seq).find(Rseq3)
        inEnd = (KOshift.seq).find(Fseq5rc) + int(Fleng5)

    finalmap = leftRec + KOshift[inStart:inEnd] + rightRec

##it would be nice to merge fragmented features at record junctions

    finalmap.name = gene + "_BsR"
    return finalmap

##############
#main program#
##############
genefile = args.genes

## perhaps give user option to set flank
flank = 5000 #int(raw_input("How much flanking sequence would you like to grab?: \n"))

log = open("log.txt", "w")

genelist = []

if "txt" in genefile:
    with open(genefile, "r") as f:
        genelist = f.readlines()
else:
    genelist.append(genefile)

##initiate output file
with open('gene_specific_primers.csv', 'w') as pf:
    pf.write('GeneHA,Left_seq,Left_tm,Left_start,Left_length,Right_seq,Right_tm,Right_start,Right_len,HA_len\n')

p2o = open('primers_to_order.csv', 'w')
p2o.write('Primer,Seq,GSP_tm,Subclone_RE,TypeII_RE,amplicon_length\n')

##cycle through genes in input
##and call subroutines to get genomic info from NCBI
for gene in genelist:
    gene = gene.translate(None,"\t\n\r") #remove end of line crap

    newdir = gene + "_outfiles"
    os.mkdir(newdir)

    gbrecord = locus_maps(gene, flank) ###returns outfile name to use in primer stuff
    if gbrecord == "error":
        os.rmdir(newdir)
        continue

##get sequence record from GB file and do stuff with it
    seqpath = newdir+"/"+gbrecord
    seqRecord = SeqIO.read(seqpath, "genbank")            #read in seqRecord for gene +- 5kb
    sequence = str(seqRecord.seq)                        #get actual sequence
    HA5pos = 5000                                        #set HA anchor positions
    HA3pos = len(sequence) - 5000
    geneHA5 = gene + "_HA5"
    geneHA3 = gene + "_HA3"

##design primers for homology arm amplification
##this module (runp3) has failed for 1 gene (gtaC), with the primer3 program getting stuck
#    need to figure out what went wrong in that case, and how to escape the error
    HA5_primers = runp3(geneHA5, sequence, HA5pos, 5)    #get primers for HA5
    if HA5_primers == "error":
        continue

##parse primer3 output
    HA5_Lstart, HA5_Lleng, HA5_Rstart, HA5_Rleng, fwd5TM, rev5TM, Fseq5, Rseq5 = parsep3(HA5_primers, geneHA5) #parse relevant info

    HA3_primers = runp3(geneHA3, sequence, HA3pos, 3)    #get primers for HA3
    if HA3_primers == "error":
        continue

    HA3_Lstart, HA3_Lleng, HA3_Rstart, HA3_Rleng, fwd3TM, rev3TM, Fseq3, Rseq3 = parsep3(HA3_primers, geneHA3)    #parse relevant info

##define homology arms based on primer locations calculated above
    HA5_seq = sequence[int(HA5_Lstart):(int(HA5_Rstart)+1)]
    HA3_seq = sequence[int(HA3_Lstart):(int(HA3_Rstart)+1)]

##check availability of RE sites in HA's for the REs used in the pLPBLP KO vector (Dicty field standard)
    REerror,addRE55,addRE53,addRE35,addRE33,typeII, HAorder = HA_REs(HA5_seq, HA3_seq)     #check for available RE sites in HAs
    if REerror == 1:
        continue

##Build GB files for homology arms
    HA5_fragment = GBslice(seqRecord, HA5_seq)        #get annotated slice of seqRecord
    HA5_record, fwd5, rev5 = addREprimers(HA5_fragment, addRE55, addRE53, HA5_Lleng, HA5_Rleng, typeII, HAorder, 5) #build HA amplicon GB
    HA5_record.name = gene + "_HA5"
    SeqIO.write(HA5_record, newdir+"/"+gene+"_HA5.gb", "genbank")
    HA3_fragment = GBslice(seqRecord, HA3_seq)        #get annotated slice of seqRecord
    HA3_record, fwd3, rev3 = addREprimers(HA3_fragment, addRE35, addRE33, HA3_Lleng, HA3_Rleng, typeII, HAorder, 3) #build HA amplicon GB
    HA3_record.name = gene + "_HA3"
    SeqIO.write(HA3_record, newdir+"/"+gene+"_HA3.gb", "genbank")

##Merge pLPBLP vector map with HA maps
    pLPBLP_HAs = pLPBLP_merge(HA5_record, HA3_record, addRE55, addRE53 ,addRE35 ,addRE33, HAorder)
    pNAME = gene + "_pLPBLP.gb"
    pLPBLP_HAs.name = gene + "_KO"
    SeqIO.write(pLPBLP_HAs,  newdir+"/"+pNAME, "genbank")

##Construct and write final insert map
    #insert_map = seqRecord[:(int(HA5_Rstart)+1)] + BsR_insert + seqRecord[int(HA3_Lstart):]
    final_map = build_final_map(seqRecord, pLPBLP_HAs, Fseq5, HA5_Lleng, Rseq3, HA3_Rleng, HAorder)
    mapname = gene + "_locus_BsR_insertion.gb"
    SeqIO.write(final_map,  newdir+"/"+mapname, "genbank")

##print primers_to_order
    ha5f = str(gene+'_HA5_fwd,'+fwd5+','+str(fwd5TM)+','+str(addRE55)+','+str(typeII)+','+str(len(HA5_record))+' bp\n')
    ha5r = str(gene+'_HA5_rev,'+rev5+','+str(rev5TM)+','+str(addRE53)+',,\n')
    ha3f = str(gene+'_HA3_fwd,'+fwd3+','+str(fwd3TM)+','+str(addRE35)+',,'+str(len(HA3_record))+' bp\n')
    ha3r = str(gene+'_HA3_rev,'+rev3+','+str(rev3TM)+','+str(addRE33)+','+str(typeII)+'\n')
    p2o.write(ha5f)
    p2o.write(ha5r)
    p2o.write(ha3f)
    p2o.write(ha3r)

log.close()
pf.close()
p2o.close()








